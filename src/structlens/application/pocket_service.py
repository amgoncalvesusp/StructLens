"""Application boundary for deterministic blind pocket detection.

The application layer owns input selection, progress/cancellation and
provenance.  Alpha-sphere construction, solvent exposure and candidate
ranking remain in :mod:`structlens.core.pockets`; this keeps GUI and file
system concerns out of the geometric implementation.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from importlib.metadata import version as package_version
from threading import Event
from types import MappingProxyType
from typing import cast

from structlens.core.errors import AnalysisCancelledError
from structlens.core.evidence import Availability, Diagnostic, DiagnosticSeverity
from structlens.core.models import AtomRecord, ProteinChain, ResidueId
from structlens.core.parsing import ParsedStructure
from structlens.core.pockets import (
    POCKET_RADII_VERSION,
    PocketCandidate,
    PocketDetectionSettings,
    vdw_radius_angstrom,
)
from structlens.core.pockets.clustering import cluster_alpha_spheres
from structlens.core.pockets.delaunay import PocketAtom, detect_alpha_spheres
from structlens.core.pockets.ranking import rank_pocket_candidates
from structlens.core.provenance import FrozenJSON, MethodProvenance

ProgressCallback = Callable[[str], None]

_METHOD_ID = "structlens.pocket.detect"
_METHOD_VERSION = "0.4.0"
_HYDROGEN_ELEMENTS = frozenset({"H", "D", "T"})
_NUMERICAL_FAILURE_CODES = frozenset(
    {
        "pocket.detect.resource_limit",
        "pocket.detect.qhull_failure",
    }
)


@dataclass(frozen=True, slots=True)
class PocketDetectionReport:
    """Immutable result of one blind pocket-detection run.

    ``NOT_DETECTED`` is intentionally distinct from a numerical or invalid
    input failure.  Empty candidates therefore never masquerade as a zero
    volume or a successful measurement.
    """

    availability: Availability
    candidates: tuple[PocketCandidate, ...] = ()
    diagnostics: tuple[Diagnostic, ...] = ()
    counts: Mapping[str, int] = field(default_factory=dict)
    provenance: MethodProvenance | None = None

    def __post_init__(self) -> None:
        availability = self.availability
        if not isinstance(availability, Availability):
            try:
                availability = Availability(availability)
            except (TypeError, ValueError) as exc:
                raise ValueError(f"unknown pocket availability: {self.availability!r}") from exc
            object.__setattr__(self, "availability", availability)
        candidates = tuple(self.candidates)
        if any(not isinstance(item, PocketCandidate) for item in candidates):
            raise TypeError("candidates must contain PocketCandidate values")
        diagnostics = tuple(self.diagnostics)
        if any(not isinstance(item, Diagnostic) for item in diagnostics):
            raise TypeError("diagnostics must contain Diagnostic values")
        counts: dict[str, int] = {}
        if not isinstance(self.counts, Mapping):
            raise TypeError("counts must be a mapping")
        for key, value in self.counts.items():
            if not isinstance(key, str) or not key.strip():
                raise ValueError("count names must be non-empty strings")
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ValueError("counts must contain non-negative integers")
            counts[key] = value
        if "candidates" not in counts:
            counts["candidates"] = len(candidates)
        elif counts["candidates"] != len(candidates):
            raise ValueError("counts['candidates'] must match the candidate collection")
        object.__setattr__(self, "candidates", candidates)
        object.__setattr__(self, "diagnostics", diagnostics)
        object.__setattr__(self, "counts", MappingProxyType(counts))
        if self.provenance is not None and not isinstance(self.provenance, MethodProvenance):
            raise TypeError("provenance must be a MethodProvenance or None")

    @property
    def candidate_count(self) -> int:
        """Number of displayed, ranked candidates."""

        return len(self.candidates)

    def to_json(self) -> dict[str, object]:
        """Return a fresh JSON-ready representation for export adapters."""

        return {
            "availability": self.availability.value,
            "candidates": [candidate.to_json() for candidate in self.candidates],
            "diagnostics": [diagnostic.to_json() for diagnostic in self.diagnostics],
            "counts": dict(self.counts),
            "provenance": self.provenance.to_json() if self.provenance is not None else None,
        }


class StructurePocketService:
    """Detect blind pockets from the selected polymer heavy atoms."""

    def __init__(self, settings: PocketDetectionSettings | None = None) -> None:
        self._settings = settings or PocketDetectionSettings()

    @property
    def settings(self) -> PocketDetectionSettings:
        return self._settings

    def analyze(
        self,
        parsed: ParsedStructure,
        *,
        progress_callback: ProgressCallback | None = None,
        cancel_event: Event | None = None,
    ) -> PocketDetectionReport:
        """Run the four auditable detection stages on one parsed structure."""

        if not isinstance(parsed, ParsedStructure):
            raise TypeError("parsed must be a ParsedStructure")
        _check_cancel(cancel_event)
        _progress(progress_callback, "prepare")
        _check_cancel(cancel_event)
        atoms, preparation_diagnostics, excluded_counts = _selected_polymer_atoms(parsed)
        counts = {"input_atoms": len(atoms), **excluded_counts}
        provenance = _build_provenance(parsed, self._settings, input_atoms=len(atoms))
        if len(atoms) < 4:
            diagnostics = preparation_diagnostics + (
                _diagnostic(
                    "pocket.detect.insufficient_atoms",
                    "Fewer than four selected polymer heavy atoms with validated radii were available.",
                    source_id=parsed.selection.display_name,
                    remediation="Select a protein or nucleic-acid polymer with at least four known heavy atoms.",
                ),
            )
            return PocketDetectionReport(
                availability=Availability.INVALID_INPUT,
                diagnostics=diagnostics,
                counts=counts,
                provenance=provenance,
            )

        _check_cancel(cancel_event)
        _progress(progress_callback, "tessellate")
        _check_cancel(cancel_event)
        alpha_result = detect_alpha_spheres(
            atoms,
            self._settings,
            cancel_check=lambda: cancel_event is not None and cancel_event.is_set(),
        )
        _check_cancel(cancel_event)
        alpha_spheres = tuple(alpha_result.spheres)
        counts["alpha_spheres"] = len(alpha_spheres)
        diagnostics = preparation_diagnostics + tuple(alpha_result.diagnostics)
        if _has_numerical_failure(alpha_result.diagnostics):
            return PocketDetectionReport(
                availability=Availability.NUMERICAL_FAILURE,
                diagnostics=diagnostics,
                counts=counts,
                provenance=provenance,
            )

        _check_cancel(cancel_event)
        _progress(progress_callback, "cluster")
        _check_cancel(cancel_event)
        cluster_result = cluster_alpha_spheres(alpha_spheres, atoms, self._settings)
        diagnostics += tuple(cluster_result.diagnostics)
        counts["clusters"] = len(cluster_result.candidates)
        if _has_numerical_failure(cluster_result.diagnostics):
            return PocketDetectionReport(
                availability=Availability.NUMERICAL_FAILURE,
                diagnostics=diagnostics,
                counts=counts,
                provenance=provenance,
            )

        _check_cancel(cancel_event)
        _progress(progress_callback, "rank")
        _check_cancel(cancel_event)
        ranked = rank_pocket_candidates(
            tuple(cluster_result.candidates),
            maximum_candidates=self._settings.maximum_candidates,
        )
        counts["candidates"] = len(ranked)
        if not ranked:
            return PocketDetectionReport(
                availability=Availability.NOT_DETECTED,
                diagnostics=diagnostics,
                counts=counts,
                provenance=provenance,
            )
        return PocketDetectionReport(
            availability=Availability.AVAILABLE,
            candidates=ranked,
            diagnostics=diagnostics,
            counts=counts,
            provenance=provenance,
        )


def detect_blind_pockets(
    parsed: ParsedStructure,
    settings: PocketDetectionSettings | None = None,
    *,
    progress_callback: ProgressCallback | None = None,
    cancel_event: Event | None = None,
) -> PocketDetectionReport:
    """Convenience facade for one blind pocket-detection run."""

    return StructurePocketService(settings).analyze(
        parsed,
        progress_callback=progress_callback,
        cancel_event=cancel_event,
    )


def _selected_polymer_atoms(
    parsed: ParsedStructure,
) -> tuple[tuple[PocketAtom, ...], tuple[Diagnostic, ...], dict[str, int]]:
    selected_records = sorted(
        (
            record
            for chain in parsed.protein_structure.chains
            if _chain_is_selected(chain, parsed)
            for record in chain.residue_records
        ),
        key=lambda record: (
            record.residue_id.chain_id,
            record.residue_id.auth_seq_id,
            record.residue_id.insertion_code or "",
        ),
    )
    atoms: list[PocketAtom] = []
    diagnostics: list[Diagnostic] = []
    used_ids: set[str] = set()
    excluded_hydrogens = 0
    excluded_unknown_radii = 0
    excluded_missing_residue = 0
    for record in selected_records:
        residue_id = record.residue_id
        residue_label = _residue_label(residue_id)
        for index, atom in enumerate(
            sorted(record.atoms, key=lambda item: (item.name, item.source_atom_id or "", item.coordinate))
        ):
            element = atom.element.strip().upper()
            if element in _HYDROGEN_ELEMENTS:
                excluded_hydrogens += 1
                continue
            atom_id = _atom_id(residue_label, atom, index, used_ids)
            radius = vdw_radius_angstrom(element)
            if radius is None or not math.isfinite(radius) or radius <= 0.0:
                excluded_unknown_radii += 1
                diagnostics.append(
                    _diagnostic(
                        "pocket.detect.unknown_radius",
                        f"No validated pocket radius is available for element {element!r}; atom was excluded.",
                        source_id=residue_label,
                        atom_id=atom_id,
                        residue_id=_residue_label(residue_id),
                        remediation="Review the element annotation or exclude the atom from pocket geometry.",
                    )
                )
                continue
            atoms.append(
                PocketAtom(
                    atom_id=atom_id,
                    residue_id=residue_id,
                    coordinate=cast(tuple[float, float, float], tuple(atom.coordinate)),
                    radius_angstrom=radius,
                )
            )
    return (
        tuple(sorted(atoms, key=lambda item: item.atom_id)),
        tuple(_sorted_diagnostics(diagnostics)),
        {
            "excluded_hydrogens": excluded_hydrogens,
            "excluded_unknown_radii": excluded_unknown_radii,
            "excluded_missing_residue": excluded_missing_residue,
        },
    )


def _chain_is_selected(chain: ProteinChain, parsed: ParsedStructure) -> bool:
    if str(chain.model_id) != parsed.selection.model_id:
        return False
    locators = parsed.selection.chain_locators
    if not locators:
        return True
    author_chain_id = chain.author_chain_id if chain.author_chain_id is not None else chain.chain_id
    return any(
        (locator.author_chain_id is None or locator.author_chain_id == author_chain_id)
        and (locator.label_chain_id is None or locator.label_chain_id == chain.label_chain_id)
        and (locator.entity_id is None or locator.entity_id == chain.entity_id)
        for locator in locators
    )


def _atom_id(component_id: str, atom: AtomRecord, index: int, used_ids: set[str]) -> str:
    base = atom.source_atom_id or f"{component_id}:{atom.name}:{index}"
    candidate = str(base).strip()
    if not candidate:
        candidate = f"{component_id}:{atom.name}:{index}"
    if candidate in used_ids:
        candidate = f"{component_id}:{candidate}:{index}"
    suffix = 2
    original = candidate
    while candidate in used_ids:
        candidate = f"{original}:{suffix}"
        suffix += 1
    used_ids.add(candidate)
    return candidate


def _build_provenance(
    parsed: ParsedStructure,
    settings: PocketDetectionSettings,
    *,
    input_atoms: int,
) -> MethodProvenance:
    settings_json = settings.to_json()
    parameters: dict[str, object] = {
        **settings_json,
        "input_atoms": input_atoms,
        "selection_id": parsed.selection.selection_id,
        "selection": {
            "model_id": parsed.selection.model_id,
            "author_chain_ids": parsed.selection.author_chain_ids,
            "label_chain_ids": parsed.selection.label_chain_ids,
            "chain_locators": tuple(
                {
                    "author_chain_id": locator.author_chain_id,
                    "label_chain_id": locator.label_chain_id,
                    "entity_id": locator.entity_id,
                }
                for locator in parsed.selection.chain_locators
            ),
            "altloc_policy": parsed.selection.altloc_policy.value,
            "assembly_scope": parsed.selection.assembly_scope.value,
        },
        "atom_scope": "selected_polymer_heavy_atoms",
        "component_scope": "selected_model_and_chain_polymer_residue_records",
        "solvent_exposure_policy": "grid_probe_connectivity_with_boundary_rejection",
    }
    return MethodProvenance(
        method_id=_METHOD_ID,
        method_version=_METHOD_VERSION,
        parameters=cast(Mapping[str, FrozenJSON], parameters),
        units=_settings_units(settings_json),
        backend_versions={
            "scipy": package_version("scipy"),
            "pocket_radii": POCKET_RADII_VERSION,
        },
        input_hashes={
            "raw_source": parsed.raw_source_hash,
            "logical_content": parsed.selection.content_id,
        },
        analyzed_representation=parsed.selection.assembly_scope.value,
    )


def _settings_units(payload: Mapping[str, object], prefix: str = "") -> dict[str, str]:
    units: dict[str, str] = {}
    for key, value in payload.items():
        path = f"{prefix}.{key}" if prefix else key
        if isinstance(value, Mapping):
            units.update(_settings_units(value, path))
            continue
        leaf = key.casefold()
        if "angstrom" in leaf or leaf.endswith("_radius") or "spacing" in leaf or "padding" in leaf or "margin" in leaf:
            units[path] = "angstrom"
        elif (
            "count" in leaf
            or "candidate" in leaf
            or "atom" in leaf
            or "simplex" in leaf
            or "cluster_size" in leaf
            or "cell" in leaf
            or "check" in leaf
        ):
            units[path] = "count"
    return units


def _diagnostic(
    code: str,
    message: str,
    *,
    source_id: str | None = None,
    atom_id: str | None = None,
    residue_id: str | None = None,
    remediation: str | None = None,
) -> Diagnostic:
    severity = DiagnosticSeverity.WARNING if code.endswith("unknown_radius") else DiagnosticSeverity.ERROR
    return Diagnostic(
        code=code,
        severity=severity,
        message=message,
        source_id=source_id,
        atom_id=atom_id,
        residue_id=residue_id,
        remediation=remediation,
    )


def _residue_label(residue_id: ResidueId) -> str:
    insertion = residue_id.insertion_code or ""
    return ":".join(
        (
            residue_id.structure_id,
            residue_id.model_id,
            residue_id.chain_id,
            f"{residue_id.auth_seq_id}{insertion}",
            residue_id.residue_name,
        )
    )


def _sorted_diagnostics(diagnostics: Sequence[Diagnostic]) -> tuple[Diagnostic, ...]:
    return tuple(
        sorted(
            diagnostics,
            key=lambda item: (
                item.code,
                item.source_id or "",
                item.residue_id or "",
                item.atom_id or "",
                item.message,
            ),
        )
    )


def _has_numerical_failure(diagnostics: Sequence[Diagnostic]) -> bool:
    return any(item.code in _NUMERICAL_FAILURE_CODES for item in diagnostics)


def _progress(callback: ProgressCallback | None, stage: str) -> None:
    if callback is not None:
        callback(stage)


def _check_cancel(cancel_event: Event | None) -> None:
    if cancel_event is not None and cancel_event.is_set():
        raise AnalysisCancelledError("Pocket detection cancelled by the user")


__all__ = [
    "PocketDetectionReport",
    "ProgressCallback",
    "StructurePocketService",
    "detect_blind_pockets",
]
