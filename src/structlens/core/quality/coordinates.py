"""Deterministic atomic and residue-level coordinate quality checks.

The functions in this module are intentionally pure.  They screen evidence;
they do not repair atoms, drop warnings, or make a structural interpretation.
"""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from types import MappingProxyType

from structlens.core.evidence.status import Availability, Diagnostic, DiagnosticSeverity
from structlens.core.models.residue import ResidueId
from structlens.core.models.structure import AtomRecord, ProteinStructure, ResidueRecord
from structlens.core.provenance import MethodProvenance

from .models import CoordinateAtom, CoordinateQCSettings, CoordinateResidue, StructureQualityReport
from .radii import vdw_radius_angstrom

_BACKBONE = frozenset(("N", "CA", "C", "O"))
_DIAGNOSTIC_ORDER = MappingProxyType(
    {
        "element.missing": 5,
        "coordinate.invalid": 10,
        "coordinate.non_finite": 10,
        "occupancy.invalid": 20,
        "b_factor.invalid": 30,
        "element.unknown_radius": 40,
        "altloc.occupancy_sum": 50,
        "atom.duplicate_identity": 60,
        "backbone.missing_atoms": 70,
    }
)

@dataclass(frozen=True, slots=True)
class _AtomEntry:
    atom: CoordinateAtom
    residue_id: ResidueId | None
    ordinal: int


def _as_atom(value: object) -> CoordinateAtom:
    if isinstance(value, CoordinateAtom):
        return value
    if isinstance(value, AtomRecord):
        return CoordinateAtom.from_atom_record(value)
    raise TypeError("coordinate inputs must contain AtomRecord or CoordinateAtom values")


def _residue_text(residue_id: ResidueId | None) -> str | None:
    if residue_id is None:
        return None
    insertion = residue_id.insertion_code or ""
    return ":".join(
        (
            str(residue_id.structure_id),
            str(residue_id.model_id),
            str(residue_id.chain_id),
            f"{residue_id.auth_seq_id}{insertion}",
            str(residue_id.residue_name),
        )
    )


def _atom_text(entry: _AtomEntry) -> str:
    atom = entry.atom
    residue = _residue_text(entry.residue_id)
    name = atom.name or "<unnamed>"
    altloc = f"[{atom.altloc}]" if atom.altloc else ""
    if residue is not None:
        return f"{residue}:{name}{altloc}"
    if atom.source_atom_id:
        return atom.source_atom_id
    return f"{name}{altloc}"


def _identity_key(entry: _AtomEntry) -> tuple[str | None, str, str | None]:
    return (_residue_text(entry.residue_id), entry.atom.name, entry.atom.altloc)


def _entry_sort_key(entry: _AtomEntry) -> tuple[str, str, int]:
    return (_residue_text(entry.residue_id) or "", _atom_text(entry), entry.ordinal)


def _entries_for_source(
    source: object,
) -> tuple[tuple[_AtomEntry, ...], tuple[ResidueId, ...], int, frozenset[ResidueId]]:
    """Flatten supported structure/residue/atom inputs without scientific loss."""

    entries: list[_AtomEntry] = []
    residues: list[ResidueId] = []
    polymer_residues: set[ResidueId] = set()
    chain_count = 0

    def add_residue(
        residue_id: ResidueId, atoms: Iterable[object], *, is_polymer: bool = True
    ) -> None:
        residues.append(residue_id)
        if is_polymer:
            polymer_residues.add(residue_id)
        for atom in atoms:
            entries.append(_AtomEntry(_as_atom(atom), residue_id, len(entries)))

    if isinstance(source, ProteinStructure):
        chain_count = len(source.chains)
        for chain in source.chains:
            if chain.residue_records:
                for record in chain.residue_records:
                    add_residue(record.residue_id, record.atoms)
            else:
                for residue_id in chain.residues:
                    residues.append(residue_id)
                    polymer_residues.add(residue_id)
    elif isinstance(source, ResidueRecord):
        add_residue(source.residue_id, source.atoms)
        chain_count = 1
    elif isinstance(source, CoordinateResidue):
        add_residue(source.residue_id, source.atoms, is_polymer=source.is_polymer)
        chain_count = 1
    elif isinstance(source, (AtomRecord, CoordinateAtom)):
        entries.append(_AtomEntry(_as_atom(source), None, 0))
    else:
        if isinstance(source, (str, bytes)):
            raise TypeError("coordinate input must not be text")
        if not isinstance(source, Iterable):
            raise TypeError("coordinate input must be a structure, residue, atom, or iterable")
        values = tuple(source)
        if values and all(isinstance(value, (AtomRecord, CoordinateAtom)) for value in values):
            entries.extend(_AtomEntry(_as_atom(value), None, index) for index, value in enumerate(values))
        else:
            for value in values:
                if isinstance(value, ResidueRecord):
                    add_residue(value.residue_id, value.atoms)
                elif isinstance(value, CoordinateResidue):
                    add_residue(value.residue_id, value.atoms, is_polymer=value.is_polymer)
                else:
                    raise TypeError("iterable coordinate input must contain atoms or residues")
            chain_count = 1 if residues else 0
    return tuple(entries), tuple(residues), chain_count, frozenset(polymer_residues)


def _diagnostic(
    code: str,
    severity: DiagnosticSeverity,
    message: str,
    *,
    source_id: str | None,
    entry: _AtomEntry | None = None,
    residue_id: ResidueId | None = None,
    remediation: str,
) -> Diagnostic:
    context = entry
    if context is not None:
        residue_id = context.residue_id
    return Diagnostic(
        code=code,
        severity=severity,
        message=message,
        source_id=source_id,
        atom_id=_atom_text(context) if context is not None else None,
        residue_id=_residue_text(residue_id),
        remediation=remediation,
    )


def _as_float(value: object) -> float | None:
    if not isinstance(value, (str, int, float)):
        return None
    try:
        return float(value)
    except (TypeError, ValueError, OverflowError):
        return None


def _atomic_diagnostics(
    entries: Sequence[_AtomEntry],
    *,
    source_id: str | None,
    settings: CoordinateQCSettings,
) -> list[Diagnostic]:
    ordered = tuple(sorted(entries, key=_entry_sort_key))
    diagnostics: list[Diagnostic] = []

    for entry in ordered:
        if not entry.atom.element:
            diagnostics.append(
                _diagnostic(
                    "element.missing",
                    DiagnosticSeverity.ERROR,
                    f"Atom {_atom_text(entry)!r} has no explicit element symbol.",
                    source_id=source_id,
                    entry=entry,
                    remediation="Provide the source element symbol; do not infer it from the atom name.",
                )
            )

    for entry in ordered:
        coordinate = entry.atom.coordinate
        if len(coordinate) != 3:
            diagnostics.append(
                _diagnostic(
                    "coordinate.invalid",
                    DiagnosticSeverity.ERROR,
                    f"Atom {_atom_text(entry)!r} does not contain exactly three coordinate values.",
                    source_id=source_id,
                    entry=entry,
                    remediation="Inspect the source atom record and provide an x, y, and z coordinate.",
                )
            )
            continue
        values = tuple(_as_float(value) for value in coordinate)
        if any(value is None for value in values):
            diagnostics.append(
                _diagnostic(
                    "coordinate.invalid",
                    DiagnosticSeverity.ERROR,
                    f"Atom {_atom_text(entry)!r} contains a non-numeric coordinate value.",
                    source_id=source_id,
                    entry=entry,
                    remediation="Inspect the source atom record and provide numeric coordinates.",
                )
            )
        elif any(not math.isfinite(value) for value in values if value is not None):
            diagnostics.append(
                _diagnostic(
                    "coordinate.non_finite",
                    DiagnosticSeverity.ERROR,
                    f"Atom {_atom_text(entry)!r} contains a non-finite coordinate.",
                    source_id=source_id,
                    entry=entry,
                    remediation="Remove NaN or infinite coordinate values before structural analysis.",
                )
            )

    for entry in ordered:
        occupancy = _as_float(entry.atom.occupancy) if entry.atom.occupancy is not None else None
        if entry.atom.occupancy is not None and (
            occupancy is None
            or not math.isfinite(occupancy)
            or occupancy < settings.occupancy_minimum
            or occupancy > settings.occupancy_maximum
        ):
            diagnostics.append(
                _diagnostic(
                    "occupancy.invalid",
                    DiagnosticSeverity.ERROR,
                    f"Atom {_atom_text(entry)!r} has occupancy outside the inclusive "
                    f"[{settings.occupancy_minimum:g}, {settings.occupancy_maximum:g}] range.",
                    source_id=source_id,
                    entry=entry,
                    remediation="Check the source occupancy and keep it finite and within [0, 1].",
                )
            )

    for entry in ordered:
        b_factor = _as_float(entry.atom.b_factor) if entry.atom.b_factor is not None else None
        if entry.atom.b_factor is not None and (
            b_factor is None or not math.isfinite(b_factor) or b_factor < settings.b_factor_minimum
        ):
            diagnostics.append(
                _diagnostic(
                    "b_factor.invalid",
                    DiagnosticSeverity.ERROR,
                    f"Atom {_atom_text(entry)!r} has a non-finite B-factor or a value below "
                    f"the declared minimum {settings.b_factor_minimum:g}.",
                    source_id=source_id,
                    entry=entry,
                    remediation="Check the source B-factor and keep it finite and non-negative.",
                )
            )

    for entry in ordered:
        if entry.atom.element and vdw_radius_angstrom(entry.atom.element) is None:
            diagnostics.append(
                _diagnostic(
                    "element.unknown_radius",
                    DiagnosticSeverity.WARNING,
                    f"Atom {_atom_text(entry)!r} uses element {entry.atom.element!r} without a known radius.",
                    source_id=source_id,
                    entry=entry,
                    remediation="Provide a supported element or treat radius-dependent checks as unavailable.",
                )
            )

    seen_identity: dict[tuple[str | None, str, str | None], _AtomEntry] = {}
    seen_source_ids: dict[tuple[str | None, str], _AtomEntry] = {}
    for entry in ordered:
        identity = _identity_key(entry)
        duplicate = identity in seen_identity
        if not duplicate:
            seen_identity[identity] = entry
        source_atom_id = entry.atom.source_atom_id
        if source_atom_id:
            model_id = entry.residue_id.model_id if entry.residue_id is not None else None
            source_key = (model_id, source_atom_id)
            if source_key in seen_source_ids:
                duplicate = True
            else:
                seen_source_ids[source_key] = entry
        if duplicate:
            diagnostics.append(
                _diagnostic(
                    "atom.duplicate_identity",
                    DiagnosticSeverity.ERROR,
                    f"Atom identity {_atom_text(entry)!r} occurs more than once.",
                    source_id=source_id,
                    entry=entry,
                    remediation="Deduplicate source atom identities before analysis.",
                )
            )
    return diagnostics


def _residue_diagnostics(
    residues: Sequence[tuple[ResidueId, tuple[_AtomEntry, ...]]],
    *,
    source_id: str | None,
    settings: CoordinateQCSettings,
    polymer_residues: frozenset[ResidueId] = frozenset(),
) -> list[Diagnostic]:
    diagnostics: list[Diagnostic] = []
    for residue_id, entries in sorted(residues, key=lambda item: _residue_text(item[0]) or ""):
        altloc_groups: dict[str, list[float]] = defaultdict(list)
        altloc_labels: dict[str, set[str]] = defaultdict(set)
        for entry in entries:
            if entry.atom.altloc:
                altloc_labels[entry.atom.name].add(entry.atom.altloc)
                occupancy = _as_float(entry.atom.occupancy) if entry.atom.occupancy is not None else None
                if occupancy is not None and math.isfinite(occupancy) and 0.0 <= occupancy <= 1.0:
                    altloc_groups[entry.atom.name].append(occupancy)
        for name in sorted(altloc_labels):
            if len(altloc_labels[name]) < 2:
                continue
            total = sum(altloc_groups[name])
            if total > 1.0 + settings.altloc_occupancy_tolerance:
                representative = next(
                    entry for entry in entries if entry.atom.name == name and entry.atom.altloc is not None
                )
                diagnostics.append(
                    _diagnostic(
                        "altloc.occupancy_sum",
                        DiagnosticSeverity.WARNING,
                        f"Alternate locations for atom {name!r} sum to {total:.6g}, above 1.0.",
                        source_id=source_id,
                        entry=representative,
                        remediation="Review alternate-location occupancies and source conformer assignment.",
                    )
                )
        if residue_id not in polymer_residues:
            continue
        atom_names = {entry.atom.name for entry in entries}
        missing = tuple(sorted(_BACKBONE - atom_names))
        if missing:
            diagnostics.append(
                _diagnostic(
                    "backbone.missing_atoms",
                    DiagnosticSeverity.WARNING,
                    f"Residue {_residue_text(residue_id)!r} is missing backbone atom(s): {', '.join(missing)}.",
                    source_id=source_id,
                    residue_id=residue_id,
                    remediation="Confirm that the residue is complete before backbone geometry analysis.",
                )
            )
    return diagnostics


def _sort_diagnostics(diagnostics: Iterable[Diagnostic]) -> tuple[Diagnostic, ...]:
    return tuple(
        sorted(
            diagnostics,
            key=lambda item: (
                _DIAGNOSTIC_ORDER.get(item.code, 999),
                item.residue_id or "",
                item.atom_id or "",
                item.code,
                item.message,
            ),
        )
    )


def _entries_by_residue(
    entries: Sequence[_AtomEntry], residue_ids: Iterable[ResidueId] = ()
) -> tuple[tuple[ResidueId, tuple[_AtomEntry, ...]], ...]:
    grouped: dict[ResidueId, list[_AtomEntry]] = defaultdict(list)
    for residue_id in residue_ids:
        if residue_id not in grouped:
            grouped[residue_id] = []
    for entry in entries:
        if entry.residue_id is not None:
            grouped[entry.residue_id].append(entry)
    return tuple((residue_id, tuple(values)) for residue_id, values in grouped.items())


def diagnostics_for_atoms(
    atoms: Iterable[AtomRecord | CoordinateAtom],
    *,
    source_id: str | None = None,
    settings: CoordinateQCSettings | None = None,
) -> tuple[Diagnostic, ...]:
    """Run coordinate, occupancy, B-factor, element, and identity checks."""

    policy = settings or CoordinateQCSettings()
    entries = tuple(_AtomEntry(_as_atom(atom), None, index) for index, atom in enumerate(atoms))
    return _sort_diagnostics(_atomic_diagnostics(entries, source_id=source_id, settings=policy))


def check_residue_quality(
    residue: ResidueRecord | CoordinateResidue,
    *,
    source_id: str | None = None,
    settings: CoordinateQCSettings | None = None,
) -> tuple[Diagnostic, ...]:
    """Run atomic and residue-inventory checks for one residue."""

    policy = settings or CoordinateQCSettings()
    entries, residue_ids, _, polymer_residues = _entries_for_source(residue)
    diagnostics = _atomic_diagnostics(entries, source_id=source_id, settings=policy)
    diagnostics.extend(
        _residue_diagnostics(
            _entries_by_residue(entries, residue_ids),
            source_id=source_id,
            settings=policy,
            polymer_residues=polymer_residues,
        )
    )
    return _sort_diagnostics(diagnostics)


def check_coordinate_quality(
    source: object,
    *,
    source_id: str | None = None,
    settings: CoordinateQCSettings | None = None,
    provenance: MethodProvenance | None = None,
) -> StructureQualityReport:
    """Screen a structure or boundary snapshot and return a typed report."""

    policy = settings or CoordinateQCSettings()
    entries, residue_ids, chain_count, polymer_residues = _entries_for_source(source)
    diagnostics = _atomic_diagnostics(entries, source_id=source_id, settings=policy)
    diagnostics.extend(
        _residue_diagnostics(
            _entries_by_residue(entries, residue_ids),
            source_id=source_id,
            settings=policy,
            polymer_residues=polymer_residues,
        )
    )
    ordered = _sort_diagnostics(diagnostics)
    errors = sum(item.severity is DiagnosticSeverity.ERROR for item in ordered)
    warnings = sum(item.severity is DiagnosticSeverity.WARNING for item in ordered)
    if errors:
        availability = Availability.INVALID_INPUT
    elif not entries and not residue_ids:
        availability = Availability.NOT_APPLICABLE
    else:
        availability = Availability.AVAILABLE
    counts = {
        "atoms": len(entries),
        "residues": len(set(residue_ids)),
        "chains": chain_count,
        "errors": errors,
        "warnings": warnings,
    }
    return StructureQualityReport(
        availability=availability,
        diagnostics=ordered,
        counts=counts,
        settings=policy,
        provenance=provenance,
    )


# Public aliases make the intent clear at call sites without maintaining
# separate implementations.
check_coordinates = check_coordinate_quality
validate_coordinates = diagnostics_for_atoms


__all__ = [
    "check_coordinate_quality",
    "check_coordinates",
    "check_residue_quality",
    "diagnostics_for_atoms",
    "validate_coordinates",
    "vdw_radius_angstrom",
]
