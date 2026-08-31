"""Deterministic focused-pocket selection over detected candidates."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from structlens.core.evidence import Availability, Diagnostic, DiagnosticSeverity
from structlens.core.models import ResidueId
from structlens.core.provenance import MethodProvenance

from .ligands import PocketLigandSupport
from .models import PocketCandidate


@dataclass(frozen=True, slots=True)
class FocusedPocketSelection:
    """One focused pocket-selection outcome with explicit availability."""

    availability: Availability
    candidate: PocketCandidate | None = None
    diagnostics: tuple[Diagnostic, ...] = ()
    selection_mode: str = "key_residues"
    ligand_support: PocketLigandSupport | None = None
    provenance: MethodProvenance | None = None

    def __post_init__(self) -> None:
        availability = self.availability
        if not isinstance(availability, Availability):
            availability = Availability(availability)
            object.__setattr__(self, "availability", availability)
        if self.candidate is not None and not isinstance(self.candidate, PocketCandidate):
            raise TypeError("candidate must be a PocketCandidate or None")
        diagnostics = tuple(self.diagnostics)
        if any(not isinstance(item, Diagnostic) for item in diagnostics):
            raise TypeError("diagnostics must contain Diagnostic values")
        object.__setattr__(self, "diagnostics", diagnostics)
        mode = str(self.selection_mode).strip()
        if not mode:
            raise ValueError("selection_mode must not be empty")
        object.__setattr__(self, "selection_mode", mode)
        if self.ligand_support is not None and not isinstance(self.ligand_support, PocketLigandSupport):
            raise TypeError("ligand_support must be PocketLigandSupport or None")
        if self.provenance is not None and not isinstance(self.provenance, MethodProvenance):
            raise TypeError("provenance must be MethodProvenance or None")
        if availability is Availability.AVAILABLE and self.candidate is None:
            raise ValueError("candidate is required when availability is available")
        if availability is not Availability.AVAILABLE and self.candidate is not None:
            raise ValueError("unavailable focused selections must not carry a candidate")
        if self.ligand_support is not None and mode != "ligand_radius":
            raise ValueError("ligand_support requires ligand_radius selection_mode")

    def to_json(self) -> dict[str, object]:
        return {
            "availability": self.availability.value,
            "candidate": self.candidate.to_json() if self.candidate is not None else None,
            "diagnostics": [item.to_json() for item in self.diagnostics],
            "selection_mode": self.selection_mode,
            "ligand_support": self.ligand_support.to_json() if self.ligand_support is not None else None,
            "provenance": self.provenance.to_json() if self.provenance is not None else None,
        }


def select_candidate_by_seed_residues(
    candidates: Sequence[PocketCandidate],
    seed_residues: Sequence[ResidueId],
    *,
    selection_mode: str = "key_residues",
) -> FocusedPocketSelection:
    """Select the candidate with the best deterministic seed-residue overlap."""

    ordered = _ordered_candidates(candidates)
    residues = tuple(seed_residues)
    if any(not isinstance(residue, ResidueId) for residue in residues):
        raise TypeError("seed_residues must contain ResidueId values")
    if not residues or not ordered:
        return _no_seed_overlap(selection_mode)
    scored = [(candidate, _seed_score(candidate, residues)) for candidate in ordered]
    best_candidate, best_score = max(scored, key=lambda item: item[1])
    if best_score[0] == 0:
        return _no_seed_overlap(selection_mode)
    return FocusedPocketSelection(
        availability=Availability.AVAILABLE,
        candidate=best_candidate,
        selection_mode=selection_mode,
    )


def select_candidate_by_ligand_support(
    candidate_support: Sequence[tuple[PocketCandidate, PocketLigandSupport]],
    *,
    ligand_id: str,
) -> FocusedPocketSelection:
    """Select the best ligand-supported candidate using deterministic tie-breaks."""

    ligand_token = str(ligand_id).strip()
    if not ligand_token:
        raise ValueError("ligand_id must not be empty")
    pairs = tuple(candidate_support)
    for item in pairs:
        if (
            not isinstance(item, tuple)
            or len(item) != 2
            or not isinstance(item[0], PocketCandidate)
            or not isinstance(item[1], PocketLigandSupport)
        ):
            raise TypeError("candidate_support must contain (PocketCandidate, PocketLigandSupport) pairs")
    if any(support.component_id != ligand_token for _, support in pairs):
        raise ValueError("ligand_id must match every PocketLigandSupport component_id")
    pairs = tuple(sorted(pairs, key=lambda item: item[0].candidate_id))
    if not pairs:
        return FocusedPocketSelection(
            availability=Availability.NOT_APPLICABLE,
            diagnostics=(
                _diagnostic(
                    "pocket.focus.no_eligible_ligand",
                    "The selected ligand is not eligible for focused pocket analysis.",
                    source_id=ligand_token,
                ),
            ),
            selection_mode="ligand_radius",
        )
    best_candidate, best_support = max(
        pairs,
        key=lambda item: (
            item[1].atom_coverage_fraction,
            item[1].lining_residue_overlap_fraction,
            -item[1].ligand_center_distance_angstrom,
        ),
    )
    if best_support.atom_coverage_fraction == 0.0 and best_support.lining_residue_overlap_fraction == 0.0:
        return FocusedPocketSelection(
            availability=Availability.NOT_DETECTED,
            diagnostics=(
                _diagnostic(
                    "pocket.focus.no_ligand_supported_candidate",
                    "No focused pocket candidate covered the selected ligand or its contact residues.",
                    source_id=ligand_token,
                ),
            ),
            selection_mode="ligand_radius",
            ligand_support=best_support,
        )
    return FocusedPocketSelection(
        availability=Availability.AVAILABLE,
        candidate=best_candidate,
        selection_mode="ligand_radius",
        ligand_support=best_support,
    )


def _ordered_candidates(candidates: Sequence[PocketCandidate]) -> tuple[PocketCandidate, ...]:
    ordered = tuple(candidates)
    if any(not isinstance(item, PocketCandidate) for item in ordered):
        raise TypeError("candidates must contain PocketCandidate values")
    return tuple(sorted(ordered, key=lambda item: item.candidate_id))


def _seed_score(
    candidate: PocketCandidate,
    seed_residues: Sequence[ResidueId],
) -> tuple[int, float, float]:
    seed_set = set(seed_residues)
    lining = set(candidate.lining_residues)
    overlap = len(seed_set & lining)
    union = len(seed_set | lining)
    jaccard = overlap / union if union else 0.0
    coverage = overlap / len(seed_set) if seed_set else 0.0
    return (overlap, jaccard, coverage)


def _diagnostic(
    code: str,
    message: str,
    *,
    source_id: str | None = None,
) -> Diagnostic:
    return Diagnostic(
        code=code,
        severity=DiagnosticSeverity.WARNING,
        message=message,
        source_id=source_id,
    )


def _no_seed_overlap(selection_mode: str) -> FocusedPocketSelection:
    return FocusedPocketSelection(
        availability=Availability.NOT_DETECTED,
        diagnostics=(
            _diagnostic(
                "pocket.focus.no_seed_overlap",
                "No focused pocket candidate overlapped the selected seed residues.",
            ),
        ),
        selection_mode=selection_mode,
    )


__all__ = [
    "FocusedPocketSelection",
    "select_candidate_by_ligand_support",
    "select_candidate_by_seed_residues",
]
