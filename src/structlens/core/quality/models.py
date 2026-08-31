"""Immutable contracts used by the coordinate-quality boundary.

``AtomRecord`` intentionally rejects malformed scientific values.  A parser
must therefore report those values before it constructs an ``AtomRecord``.
``CoordinateAtom`` is the small, immutable boundary snapshot for that purpose:
it preserves finite and non-finite scalar values alike, while the QC functions
decide which values are scientifically usable.  Valid ``AtomRecord`` objects
are accepted directly too.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from numbers import Number
from types import MappingProxyType
from typing import Any, cast

from structlens.core.evidence.status import Availability, Diagnostic, DiagnosticSeverity
from structlens.core.models.residue import ResidueId
from structlens.core.models.structure import AtomRecord, ResidueRecord
from structlens.core.provenance import MethodProvenance


def _scalar(value: object) -> object:
    """Copy a scalar without coercing invalid numerical evidence away."""

    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    if isinstance(value, Number):
        try:
            return float(cast(Any, value))
        except (TypeError, ValueError, OverflowError) as exc:
            raise TypeError(f"coordinate evidence contains non-real value {value!r}") from exc
    raise TypeError(f"coordinate evidence contains unsupported value {type(value).__name__}")


def _sequence(value: object) -> tuple[object, ...]:
    if isinstance(value, (str, bytes)):
        return tuple(_scalar(item) for item in value)
    if isinstance(value, Iterable):
        values: tuple[object, ...] = tuple(_scalar(item) for item in value)
        return values
    return (_scalar(value),)


@dataclass(frozen=True, slots=True)
class CoordinateAtom:
    """Lenient immutable atom snapshot used before strict scientific parsing.

    This value is deliberately not a replacement for :class:`AtomRecord`.
    It exists so an importer can preserve an invalid coordinate, occupancy, or
    B-factor and let the QC gate emit a stable diagnostic instead of losing the
    source context to a constructor exception.
    """

    name: str
    element: str
    coordinate: tuple[object, ...] | Sequence[object]
    altloc: str | None = None
    occupancy: object | None = None
    b_factor: object | None = None
    source_atom_id: str | None = None
    source_serial: int | str | None = None
    formal_charge: object | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "name", "" if self.name is None else str(self.name).strip())
        object.__setattr__(self, "element", "" if self.element is None else str(self.element).strip().upper())
        object.__setattr__(self, "coordinate", _sequence(self.coordinate))
        if self.altloc is not None:
            altloc = str(self.altloc).strip()
            object.__setattr__(self, "altloc", altloc or None)
        if self.occupancy is not None:
            object.__setattr__(self, "occupancy", _scalar(self.occupancy))
        if self.b_factor is not None:
            object.__setattr__(self, "b_factor", _scalar(self.b_factor))
        if self.source_atom_id is not None:
            source_atom_id = str(self.source_atom_id).strip()
            object.__setattr__(self, "source_atom_id", source_atom_id or None)
        if self.source_serial is not None and isinstance(self.source_serial, str):
            source_serial = self.source_serial.strip()
            object.__setattr__(self, "source_serial", source_serial or None)
        if self.formal_charge is not None:
            object.__setattr__(self, "formal_charge", _scalar(self.formal_charge))

    @classmethod
    def from_atom_record(cls, atom: AtomRecord) -> CoordinateAtom:
        """Take a lenient snapshot of a validated ``AtomRecord``."""

        if not isinstance(atom, AtomRecord):
            raise TypeError("atom must be an AtomRecord")
        return cls(
            name=atom.name,
            element=atom.element,
            coordinate=atom.coordinate,
            altloc=atom.altloc,
            occupancy=atom.occupancy,
            b_factor=atom.b_factor,
            source_atom_id=atom.source_atom_id,
            source_serial=atom.source_serial,
            formal_charge=atom.formal_charge,
        )


@dataclass(frozen=True, slots=True)
class CoordinateResidue:
    """Lenient immutable residue snapshot for parser-boundary QC tests."""

    residue_id: ResidueId
    atoms: tuple[CoordinateAtom, ...] | Sequence[CoordinateAtom] = field(default_factory=tuple)
    is_polymer: bool = True

    def __post_init__(self) -> None:
        if not isinstance(self.residue_id, ResidueId):
            raise TypeError("residue_id must be a ResidueId")
        atoms = tuple(self.atoms)
        if any(not isinstance(atom, CoordinateAtom) for atom in atoms):
            raise TypeError("atoms must contain CoordinateAtom values")
        if not isinstance(self.is_polymer, bool):
            raise TypeError("is_polymer must be a boolean")
        object.__setattr__(self, "atoms", atoms)

    @classmethod
    def from_residue_record(cls, residue: ResidueRecord) -> CoordinateResidue:
        if not isinstance(residue, ResidueRecord):
            raise TypeError("residue must be a ResidueRecord")
        return cls(residue.residue_id, tuple(CoordinateAtom.from_atom_record(atom) for atom in residue.atoms))


@dataclass(frozen=True, slots=True)
class CoordinateQCSettings:
    """Numerical policy for deterministic coordinate-quality screens."""

    altloc_occupancy_tolerance: float = 0.01
    b_factor_minimum: float = 0.0
    occupancy_minimum: float = 0.0
    occupancy_maximum: float = 1.0
    min_cn_distance_angstrom: float = 0.0
    max_cn_distance_angstrom: float = 2.0
    min_ca_distance_angstrom: float = 2.5
    max_ca_distance_angstrom: float = 4.5
    heavy_atom_overlap_tolerance_angstrom: float = 0.4

    def __post_init__(self) -> None:
        for name in (
            "altloc_occupancy_tolerance",
            "b_factor_minimum",
            "occupancy_minimum",
            "occupancy_maximum",
            "min_cn_distance_angstrom",
            "max_cn_distance_angstrom",
            "min_ca_distance_angstrom",
            "max_ca_distance_angstrom",
            "heavy_atom_overlap_tolerance_angstrom",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and non-negative")
            object.__setattr__(self, name, value)
        if self.occupancy_minimum > self.occupancy_maximum:
            raise ValueError("occupancy_minimum must not exceed occupancy_maximum")
        if self.min_cn_distance_angstrom > self.max_cn_distance_angstrom:
            raise ValueError("minimum C-N distance must not exceed maximum C-N distance")
        if self.min_ca_distance_angstrom > self.max_ca_distance_angstrom:
            raise ValueError("minimum C-alpha distance must not exceed maximum C-alpha distance")

    def to_json(self) -> dict[str, float]:
        """Return a fresh JSON-ready copy of the scientific QC policy."""

        return {
            "altloc_occupancy_tolerance": self.altloc_occupancy_tolerance,
            "b_factor_minimum": self.b_factor_minimum,
            "occupancy_minimum": self.occupancy_minimum,
            "occupancy_maximum": self.occupancy_maximum,
            "min_cn_distance_angstrom": self.min_cn_distance_angstrom,
            "max_cn_distance_angstrom": self.max_cn_distance_angstrom,
            "min_ca_distance_angstrom": self.min_ca_distance_angstrom,
            "max_ca_distance_angstrom": self.max_ca_distance_angstrom,
            "heavy_atom_overlap_tolerance_angstrom": self.heavy_atom_overlap_tolerance_angstrom,
        }


def _freeze_counts(counts: Mapping[str, int]) -> Mapping[str, int]:
    if not isinstance(counts, Mapping):
        raise TypeError("counts must be a mapping")
    normalized: dict[str, int] = {}
    for key, value in counts.items():
        if not isinstance(key, str) or not key.strip():
            raise ValueError("count names must be non-empty strings")
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ValueError("counts must be non-negative integers")
        normalized[key] = value
    return MappingProxyType(normalized)


@dataclass(frozen=True, slots=True)
class StructureQualityReport:
    """Immutable result of the coordinate-quality gate."""

    availability: Availability
    diagnostics: tuple[Diagnostic, ...] = field(default_factory=tuple)
    counts: Mapping[str, int] = field(default_factory=dict)
    settings: CoordinateQCSettings = field(default_factory=CoordinateQCSettings)
    provenance: MethodProvenance | None = None

    def __post_init__(self) -> None:
        availability = self.availability
        if not isinstance(availability, Availability):
            try:
                availability = Availability(availability)
            except (TypeError, ValueError) as exc:
                raise ValueError(f"unknown availability: {self.availability!r}") from exc
            object.__setattr__(self, "availability", availability)
        diagnostics = tuple(self.diagnostics)
        if any(not isinstance(item, Diagnostic) for item in diagnostics):
            raise TypeError("diagnostics must contain Diagnostic values")
        if not isinstance(self.settings, CoordinateQCSettings):
            raise TypeError("settings must be CoordinateQCSettings")
        if self.provenance is not None and not isinstance(self.provenance, MethodProvenance):
            raise TypeError("provenance must be MethodProvenance or None")
        object.__setattr__(self, "diagnostics", diagnostics)
        object.__setattr__(self, "counts", _freeze_counts(self.counts))

    @property
    def status(self) -> Availability:
        """Compatibility alias for consumers that call result state ``status``."""

        return self.availability

    @property
    def error_count(self) -> int:
        return sum(item.severity is DiagnosticSeverity.ERROR for item in self.diagnostics)

    @property
    def warning_count(self) -> int:
        return sum(item.severity is DiagnosticSeverity.WARNING for item in self.diagnostics)

    def to_json(self) -> dict[str, object]:
        """Return a fresh, deterministic, JSON-ready report representation."""

        return {
            "availability": self.availability.value,
            "diagnostics": [item.to_json() for item in self.diagnostics],
            "counts": dict(self.counts),
            "settings": self.settings.to_json(),
            "provenance": self.provenance.to_json() if self.provenance is not None else None,
        }


# Names used by early integrations remain descriptive aliases, not separate
# contracts.  This makes the boundary easy to adopt without duplicating state.
CoordinateQualityReport = StructureQualityReport
QualitySettings = CoordinateQCSettings
RawCoordinateAtom = CoordinateAtom
RawCoordinateResidue = CoordinateResidue


__all__ = [
    "CoordinateAtom",
    "CoordinateQualityReport",
    "CoordinateQCSettings",
    "CoordinateResidue",
    "QualitySettings",
    "RawCoordinateAtom",
    "RawCoordinateResidue",
    "StructureQualityReport",
]
