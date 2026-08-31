"""Spatial heavy-atom overlap screening.

This module deliberately implements a *screening* measurement rather than a
MolProbity-compatible validation score.  It uses an explicit, versioned table
of van der Waals radii and a :class:`scipy.spatial.cKDTree` to find only
spatially possible pairs.  Hydrogens are not added or inferred, so the result
must be reported as ``heavy-atom overlap screening`` wherever it is shown.

The local value objects keep this module independent of the quality report
models.  An application adapter can convert parser/model atoms into
``ClashAtom`` without changing the screening contract.
"""

from __future__ import annotations

import math
import re
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import cKDTree  # type: ignore[import-untyped]

from structlens.core.models import AtomRecord

from .radii import (
    VAN_DER_WAALS_RADII_ANGSTROM,
    VDW_RADII,
    VDW_RADII_VERSION,
)

CLASH_ALGORITHM_VERSION = "structlens-heavy-atom-overlap-1"
_HYDROGEN_ELEMENTS = frozenset({"H", "D", "T"})
_RESIDUE_NUMBER_RE = re.compile(r"^\s*(-?\d+)")
_MIN_PLAUSIBLE_PEPTIDE_CN_ANGSTROM = 1.0
_MAX_PLAUSIBLE_PEPTIDE_CN_ANGSTROM = 2.0
_MIN_DISULFIDE_SS_ANGSTROM = 1.8
_MAX_DISULFIDE_SS_ANGSTROM = 2.3


@dataclass(frozen=True, slots=True)
class ClashAtom:
    """Immutable atom context used by the overlap screen.

    ``residue_id`` and ``component_id`` are opaque source identities.  When
    either identity is shared by two atoms, that pair is excluded because
    covalent geometry within one residue/component is not a clash signal.
    ``residue_number`` (or ``residue_index``) plus ``chain_id`` enables the
    explicit adjacent peptide C--N exclusion.  Without that context, no
    adjacency is guessed.
    """

    atom_id: str
    element: str
    coordinate: Sequence[float]
    atom_name: str | None = None
    residue_id: str | None = None
    component_id: str | None = None
    chain_id: str | None = None
    residue_number: int | str | None = None
    residue_index: int | None = None
    residue_name: str | None = None

    def __post_init__(self) -> None:
        atom_id = str(self.atom_id).strip()
        element = str(self.element).strip().upper()
        if not atom_id:
            raise ValueError("atom_id must not be empty")
        if not element:
            raise ValueError("atom element must not be empty")
        try:
            coordinate = tuple(float(value) for value in self.coordinate)
        except (TypeError, ValueError) as exc:
            raise ValueError("atom coordinate must contain numeric values") from exc
        if len(coordinate) != 3 or any(not math.isfinite(value) for value in coordinate):
            raise ValueError("atom coordinate must contain three finite values")
        if self.residue_index is not None:
            if isinstance(self.residue_index, bool) or not isinstance(self.residue_index, int):
                raise ValueError("residue_index must be an integer")
        if self.residue_number is not None:
            if isinstance(self.residue_number, bool) or not isinstance(self.residue_number, (int, str)):
                raise ValueError("residue_number must be an integer or string")
            if isinstance(self.residue_number, str) and not self.residue_number.strip():
                raise ValueError("residue_number must not be empty")
        if self.residue_index is not None and self.residue_number is not None:
            parsed = _parse_residue_number(self.residue_number)
            if parsed is not None and parsed != self.residue_index:
                raise ValueError("residue_number and residue_index disagree")
        object.__setattr__(self, "atom_id", atom_id)
        object.__setattr__(self, "element", element)
        object.__setattr__(self, "coordinate", coordinate)
        if self.atom_name is not None:
            object.__setattr__(self, "atom_name", str(self.atom_name).strip().upper() or None)
        for field_name in ("residue_id", "component_id", "chain_id", "residue_name"):
            value = getattr(self, field_name)
            if value is not None:
                normalized = str(value).strip()
                if field_name == "residue_name":
                    normalized = normalized.upper()
                object.__setattr__(self, field_name, normalized or None)
        if isinstance(self.residue_number, str):
            object.__setattr__(self, "residue_number", self.residue_number.strip())


@dataclass(frozen=True, slots=True)
class ClashObservation:
    """One ordered heavy-atom van der Waals overlap observation."""

    atom_a_id: str
    atom_b_id: str
    element_a: str
    element_b: str
    distance_angstrom: float
    overlap_angstrom: float
    vdw_radius_a_angstrom: float
    vdw_radius_b_angstrom: float
    residue_a_id: str | None = None
    residue_b_id: str | None = None
    component_a_id: str | None = None
    component_b_id: str | None = None

    def __post_init__(self) -> None:
        for field_name in ("atom_a_id", "atom_b_id", "element_a", "element_b"):
            value = str(getattr(self, field_name)).strip().upper() if field_name.startswith("element") else str(getattr(self, field_name)).strip()
            if not value:
                raise ValueError(f"{field_name} must not be empty")
            object.__setattr__(self, field_name, value)
        for field_name in (
            "distance_angstrom",
            "overlap_angstrom",
            "vdw_radius_a_angstrom",
            "vdw_radius_b_angstrom",
        ):
            numeric_value = float(getattr(self, field_name))
            if not math.isfinite(numeric_value) or numeric_value < 0.0:
                raise ValueError(f"{field_name} must be finite and non-negative")
            object.__setattr__(self, field_name, numeric_value)
        for field_name in ("residue_a_id", "residue_b_id", "component_a_id", "component_b_id"):
            value = getattr(self, field_name)
            if value is not None:
                object.__setattr__(self, field_name, str(value).strip() or None)


@dataclass(frozen=True, slots=True)
class ClashDiagnostic:
    """A deterministic warning explaining why an atom was not screened."""

    code: str
    severity: str
    message: str
    atom_id: str | None = None
    element: str | None = None
    context: tuple[tuple[str, str], ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        code = str(self.code).strip().upper()
        severity = str(self.severity).strip().lower()
        message = str(self.message).strip()
        if not code or not severity or not message:
            raise ValueError("diagnostic code, severity, and message are required")
        context = tuple((str(key), str(value)) for key, value in self.context)
        if any(not key.strip() for key, _ in context):
            raise ValueError("diagnostic context keys must be non-empty")
        object.__setattr__(self, "code", code)
        object.__setattr__(self, "severity", severity)
        object.__setattr__(self, "message", message)
        object.__setattr__(self, "context", context)
        if self.atom_id is not None:
            object.__setattr__(self, "atom_id", str(self.atom_id).strip() or None)
        if self.element is not None:
            object.__setattr__(self, "element", str(self.element).strip().upper() or None)


@dataclass(frozen=True, slots=True)
class ClashScreeningResult:
    """Immutable result for heavy-atom overlap screening.

    The measurement is pairwise overlap depth in Å, not a normalized score.
    It excludes hydrogen atoms and known intramolecular pairs, and reports
    unknown-radius atoms as diagnostics instead of assigning an arbitrary
    radius.  It has no claim of crystallographic validation or MolProbity
    equivalence.
    """

    observations: tuple[ClashObservation, ...] = field(default_factory=tuple)
    diagnostics: tuple[ClashDiagnostic, ...] = field(default_factory=tuple)
    input_atom_count: int = 0
    screened_atom_count: int = 0
    hydrogen_count: int = 0
    unknown_atom_count: int = 0
    overlap_tolerance_angstrom: float = 0.4
    method_name: str = "heavy-atom overlap screening"
    algorithm_version: str = CLASH_ALGORITHM_VERSION
    vdw_radii_version: str = VDW_RADII_VERSION

    def __post_init__(self) -> None:
        observations = tuple(self.observations)
        diagnostics = tuple(self.diagnostics)
        if any(not isinstance(item, ClashObservation) for item in observations):
            raise TypeError("observations must contain ClashObservation values")
        if any(not isinstance(item, ClashDiagnostic) for item in diagnostics):
            raise TypeError("diagnostics must contain ClashDiagnostic values")
        for field_name in (
            "input_atom_count",
            "screened_atom_count",
            "hydrogen_count",
            "unknown_atom_count",
        ):
            value = getattr(self, field_name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ValueError(f"{field_name} must be a non-negative integer")
        tolerance = float(self.overlap_tolerance_angstrom)
        if not math.isfinite(tolerance) or tolerance < 0.0:
            raise ValueError("overlap_tolerance_angstrom must be finite and non-negative")
        object.__setattr__(self, "observations", observations)
        object.__setattr__(self, "diagnostics", diagnostics)
        object.__setattr__(self, "overlap_tolerance_angstrom", tolerance)

    @property
    def has_warnings(self) -> bool:
        return bool(self.diagnostics)


def _parse_residue_number(value: int | str | None) -> int | None:
    if value is None:
        return None
    if isinstance(value, int) and not isinstance(value, bool):
        return value
    match = _RESIDUE_NUMBER_RE.match(str(value))
    return int(match.group(1)) if match else None


def _effective_residue_number(atom: ClashAtom) -> int | None:
    if atom.residue_index is not None:
        return atom.residue_index
    if atom.residue_number is not None:
        return _parse_residue_number(atom.residue_number)
    return None


def _is_same_residue_or_component(first: ClashAtom, second: ClashAtom) -> bool:
    return (
        first.residue_id is not None
        and first.residue_id == second.residue_id
    ) or (
        first.component_id is not None
        and first.component_id == second.component_id
    )


def _are_adjacent_residues(first: ClashAtom, second: ClashAtom) -> bool:
    if first.chain_id is None or first.chain_id != second.chain_id:
        return False
    first_number = _effective_residue_number(first)
    second_number = _effective_residue_number(second)
    return first_number is not None and second_number is not None and abs(first_number - second_number) == 1


def _is_adjacent_polymer_pair(first: ClashAtom, second: ClashAtom, distance: float) -> bool:
    if not _are_adjacent_residues(first, second):
        return False
    # Sequential residues are covalently connected and excluded as a whole.
    # The one exception is an implausibly short C--N contact: hiding that pair
    # would suppress the most direct coordinate-corruption signal.
    if {first.atom_name, second.atom_name} == {"C", "N"}:
        return _MIN_PLAUSIBLE_PEPTIDE_CN_ANGSTROM <= distance <= _MAX_PLAUSIBLE_PEPTIDE_CN_ANGSTROM
    return True


def _is_probable_disulfide(first: ClashAtom, second: ClashAtom, distance: float) -> bool:
    return (
        first.element == second.element == "S"
        and first.atom_name == second.atom_name == "SG"
        and first.residue_name == second.residue_name == "CYS"
        and _MIN_DISULFIDE_SS_ANGSTROM <= distance <= _MAX_DISULFIDE_SS_ANGSTROM
    )


def _excluded_pairs(values: Iterable[Sequence[str]]) -> frozenset[tuple[str, str]]:
    normalized: set[tuple[str, str]] = set()
    for value in values:
        pair = tuple(str(item).strip() for item in value)
        if len(pair) != 2 or not pair[0] or not pair[1] or pair[0] == pair[1]:
            raise ValueError("excluded atom pairs must contain two distinct non-empty atom IDs")
        first, second = sorted((pair[0], pair[1]))
        normalized.add((first, second))
    return frozenset(normalized)


def _coerce_atom(value: ClashAtom | AtomRecord, index: int) -> ClashAtom:
    if isinstance(value, ClashAtom):
        return value
    if isinstance(value, AtomRecord):
        atom_id = value.source_atom_id or f"atom-{index:08d}"
        return ClashAtom(atom_id, value.element, value.coordinate, atom_name=value.name)
    raise TypeError("screening inputs must be ClashAtom or AtomRecord values")


def screen_heavy_atom_overlaps(
    atoms: Iterable[ClashAtom | AtomRecord],
    *,
    overlap_tolerance_angstrom: float = 0.4,
    excluded_atom_pairs: Iterable[Sequence[str]] = (),
) -> ClashScreeningResult:
    """Screen known heavy atoms for van der Waals overlap.

    A pair is reported when its full van der Waals overlap is greater than or
    equal to ``tolerance``. The stored overlap remains the full
    ``radius_a + radius_b - distance``; the threshold is never subtracted from
    the reported measurement. ``tolerance`` is a declared geometric screening
    threshold, not a probabilistic score parameter.
    ``cKDTree.query_pairs`` bounds candidate generation by the largest
    possible radius sum, avoiding an all-pairs distance loop for sparse input.
    """

    tolerance = float(overlap_tolerance_angstrom)
    if not math.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError("overlap_tolerance_angstrom must be finite and non-negative")
    explicit_exclusions = _excluded_pairs(excluded_atom_pairs)

    normalized = tuple(_coerce_atom(value, index) for index, value in enumerate(atoms))
    diagnostics: list[ClashDiagnostic] = []
    known_atoms: list[ClashAtom] = []
    radii: list[float] = []
    seen_ids: set[str] = set()
    hydrogen_count = 0
    unknown_count = 0
    for item in normalized:
        if item.atom_id in seen_ids:
            diagnostics.append(
                ClashDiagnostic(
                    "DUPLICATE_ATOM_ID",
                    "warning",
                    "Duplicate atom identity was skipped after its first occurrence.",
                    atom_id=item.atom_id,
                    element=item.element,
                )
            )
            continue
        seen_ids.add(item.atom_id)
        if item.element in _HYDROGEN_ELEMENTS:
            hydrogen_count += 1
            continue
        radius = VDW_RADII.get(item.element)
        if radius is None:
            unknown_count += 1
            diagnostics.append(
                ClashDiagnostic(
                    "UNKNOWN_VDW_RADIUS",
                    "warning",
                    "No validated van der Waals radius is available; atom was skipped.",
                    atom_id=item.atom_id,
                    element=item.element,
                )
            )
            continue
        known_atoms.append(item)
        radii.append(radius)

    observations: list[ClashObservation] = []
    if known_atoms and tolerance < 2.0 * max(radii):
        coordinates = np.asarray([item.coordinate for item in known_atoms], dtype=float)
        tree = cKDTree(coordinates)
        max_query_distance = 2.0 * max(radii) - tolerance
        pairs = tree.query_pairs(max_query_distance, output_type="ndarray")
        for first_index, second_index in sorted((int(a), int(b)) for a, b in pairs):
            first = known_atoms[first_index]
            second = known_atoms[second_index]
            distance = float(np.linalg.norm(coordinates[first_index] - coordinates[second_index]))
            atom_pair = tuple(sorted((first.atom_id, second.atom_id)))
            if (
                _is_same_residue_or_component(first, second)
                or atom_pair in explicit_exclusions
                or _is_adjacent_polymer_pair(first, second, distance)
                or _is_probable_disulfide(first, second, distance)
            ):
                continue
            overlap = radii[first_index] + radii[second_index] - distance
            if overlap + 1.0e-12 < tolerance:
                continue
            if first.atom_id <= second.atom_id:
                atom_a, atom_b = first, second
                radius_a, radius_b = radii[first_index], radii[second_index]
            else:
                atom_a, atom_b = second, first
                radius_a, radius_b = radii[second_index], radii[first_index]
            observations.append(
                ClashObservation(
                    atom_a_id=atom_a.atom_id,
                    atom_b_id=atom_b.atom_id,
                    element_a=atom_a.element,
                    element_b=atom_b.element,
                    distance_angstrom=distance,
                    overlap_angstrom=overlap,
                    vdw_radius_a_angstrom=radius_a,
                    vdw_radius_b_angstrom=radius_b,
                    residue_a_id=atom_a.residue_id,
                    residue_b_id=atom_b.residue_id,
                    component_a_id=atom_a.component_id,
                    component_b_id=atom_b.component_id,
                )
            )

    observations.sort(key=lambda item: (item.atom_a_id, item.atom_b_id))
    diagnostics.sort(key=lambda item: (item.code, item.atom_id or "", item.element or ""))
    return ClashScreeningResult(
        observations=tuple(observations),
        diagnostics=tuple(diagnostics),
        input_atom_count=len(normalized),
        screened_atom_count=len(known_atoms),
        hydrogen_count=hydrogen_count,
        unknown_atom_count=unknown_count,
        overlap_tolerance_angstrom=tolerance,
    )


# Short aliases keep the domain vocabulary convenient while retaining the
# fully explicit function name in public documentation.
screen_clashes = screen_heavy_atom_overlaps
detect_heavy_atom_overlaps = screen_heavy_atom_overlaps
detect_clashes = screen_heavy_atom_overlaps


__all__ = [
    "CLASH_ALGORITHM_VERSION",
    "VDW_RADII",
    "VDW_RADII_VERSION",
    "VAN_DER_WAALS_RADII_ANGSTROM",
    "ClashAtom",
    "ClashDiagnostic",
    "ClashObservation",
    "ClashScreeningResult",
    "detect_clashes",
    "detect_heavy_atom_overlaps",
    "screen_clashes",
    "screen_heavy_atom_overlaps",
]
