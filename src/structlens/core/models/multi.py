"""Immutable contracts for multi-structure and integrated analyses.

These value objects deliberately contain results rather than GUI concerns.  A
reference-vs-many result keeps every target's residue locators intact; the
reference position is an anchor for display, never a replacement for source
numbering.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from enum import Enum
from types import MappingProxyType

from structlens.core.metrics.sequence_metrics import SequenceAlignmentMetrics
from structlens.core.metrics.structural_metrics import StructuralMetrics

from .correspondence import ResidueCorrespondence
from .mutation import MutationEvent
from .residue import ResidueId


class ComparisonMode(str, Enum):
    """Supported comparison topologies."""

    PAIRWISE = "pairwise"
    REFERENCE_VS_MANY = "reference_vs_many"
    ALL_VS_ALL = "all_vs_all"
    MULTIPLE_STRUCTURE_ALIGNMENT = "multiple_structure_alignment"


@dataclass(frozen=True, slots=True)
class StructuralTransform:
    """A reproducible rigid-body transform in row-major coordinates."""

    rotation: tuple[tuple[float, float, float], ...] = (
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
    )
    translation: tuple[float, float, float] = (0.0, 0.0, 0.0)

    def __post_init__(self) -> None:
        rotation = tuple(tuple(float(value) for value in row) for row in self.rotation)
        translation = tuple(float(value) for value in self.translation)
        if len(rotation) != 3 or any(len(row) != 3 for row in rotation):
            raise ValueError("rotation must be a 3x3 matrix")
        if len(translation) != 3:
            raise ValueError("translation must contain three values")
        object.__setattr__(self, "rotation", rotation)
        object.__setattr__(self, "translation", translation)


@dataclass(frozen=True, slots=True)
class TargetAnalysis:
    """One target's independent analysis anchored to a common reference."""

    target_id: str
    correspondence: tuple[ResidueCorrespondence, ...]
    mutations: tuple[MutationEvent, ...]
    sequence_metrics: SequenceAlignmentMetrics
    structural_metrics: StructuralMetrics | None = None
    transform: StructuralTransform = field(default_factory=StructuralTransform)
    provenance: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "correspondence", tuple(self.correspondence))
        object.__setattr__(self, "mutations", tuple(self.mutations))
        object.__setattr__(self, "provenance", MappingProxyType(dict(self.provenance)))

    @property
    def correspondences(self) -> tuple[ResidueCorrespondence, ...]:
        """Plural spelling used by application code and serialized payloads."""

        return self.correspondence


@dataclass(frozen=True, slots=True)
class ReferenceVsManyAnalysis:
    """Independent target analyses sharing one reference numbering system."""

    reference_id: str
    targets: Mapping[str, TargetAnalysis]
    comparison_mode: ComparisonMode = ComparisonMode.REFERENCE_VS_MANY

    def __post_init__(self) -> None:
        if not self.reference_id:
            raise ValueError("reference_id must not be empty")
        normalized = dict(self.targets)
        if set(normalized) != {analysis.target_id for analysis in normalized.values()}:
            raise ValueError("target mapping keys must match TargetAnalysis.target_id")
        object.__setattr__(self, "targets", MappingProxyType(normalized))
        if not isinstance(self.comparison_mode, ComparisonMode):
            object.__setattr__(self, "comparison_mode", ComparisonMode(self.comparison_mode))

    @property
    def target_ids(self) -> tuple[str, ...]:
        return tuple(self.targets)


@dataclass(frozen=True, slots=True)
class PairwiseMatrix:
    """A matrix that stores each unordered pair once and mirrors on read."""

    metric_name: str
    structure_ids: tuple[str, ...]
    values: Mapping[tuple[str, str], float | None]
    unit: str | None = None

    def __post_init__(self) -> None:
        structure_ids = tuple(dict.fromkeys(self.structure_ids))
        if len(structure_ids) != len(self.structure_ids):
            raise ValueError("structure_ids must be unique")
        allowed = set(structure_ids)
        canonical: dict[tuple[str, str], float | None] = {}
        for (left, right), value in self.values.items():
            if left not in allowed or right not in allowed:
                raise ValueError("matrix value references an unknown structure")
            if left == right:
                continue
            key = (left, right) if left < right else (right, left)
            canonical.setdefault(key, None if value is None else float(value))
        object.__setattr__(self, "structure_ids", structure_ids)
        object.__setattr__(self, "values", MappingProxyType(canonical))

    def value(self, left: str, right: str) -> float | None:
        if left == right:
            return None
        key = (left, right) if left < right else (right, left)
        return self.values.get(key)

    @classmethod
    def from_pairs(
        cls,
        metric_name: str,
        structure_ids: Sequence[str],
        values: Mapping[tuple[str, str], float | None],
        *,
        unit: str | None = None,
    ) -> PairwiseMatrix:
        return cls(metric_name, tuple(structure_ids), values, unit)


@dataclass(frozen=True, slots=True)
class MultiStructurePosition:
    """One reference-aligned position in a multiple-structure alignment."""

    alignment_index: int
    reference_residue: ResidueId | None
    residues: Mapping[str, ResidueId | None]
    coverage: float
    ca_positional_variability_angstrom: float | None
    per_structure_deviation_angstrom: Mapping[str, float | None] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.alignment_index < 0:
            raise ValueError("alignment_index must be non-negative")
        if not 0.0 <= self.coverage <= 1.0:
            raise ValueError("coverage must be between 0 and 1")
        object.__setattr__(self, "residues", MappingProxyType(dict(self.residues)))
        object.__setattr__(
            self,
            "per_structure_deviation_angstrom",
            MappingProxyType(dict(self.per_structure_deviation_angstrom)),
        )


@dataclass(frozen=True, slots=True)
class MultipleStructureAnalysis:
    """Reference-frame summary from a true multiple-structure alignment."""

    structure_ids: tuple[str, ...]
    aligned_positions: tuple[MultiStructurePosition, ...]
    transforms: Mapping[str, StructuralTransform]
    provenance: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        structure_ids = tuple(dict.fromkeys(self.structure_ids))
        transforms = dict(self.transforms)
        if not set(transforms).issubset(structure_ids):
            raise ValueError("transform references an unknown structure")
        object.__setattr__(self, "structure_ids", structure_ids)
        object.__setattr__(self, "aligned_positions", tuple(self.aligned_positions))
        object.__setattr__(self, "transforms", MappingProxyType(transforms))
        object.__setattr__(self, "provenance", MappingProxyType(dict(self.provenance)))


@dataclass(frozen=True, slots=True)
class AnalysisSelection:
    """Shared selection state consumed by residues, charts, and PyMOL export."""

    reference_residue: ResidueId | None = None
    target_id: str | None = None
    target_residue: ResidueId | None = None
    pair: tuple[str, str] | None = None
    key_group: str | None = None


@dataclass(frozen=True, slots=True)
class AllVsAllAnalysis:
    """All-vs-all pair results plus symmetric scientific matrices."""

    structure_ids: tuple[str, ...]
    pair_results: Mapping[tuple[str, str], object]
    matrices: Mapping[str, PairwiseMatrix]

    def __post_init__(self) -> None:
        object.__setattr__(self, "structure_ids", tuple(self.structure_ids))
        object.__setattr__(self, "pair_results", MappingProxyType(dict(self.pair_results)))
        object.__setattr__(self, "matrices", MappingProxyType(dict(self.matrices)))


__all__ = [
    "AllVsAllAnalysis",
    "AnalysisSelection",
    "ComparisonMode",
    "MultiStructurePosition",
    "MultipleStructureAnalysis",
    "PairwiseMatrix",
    "ReferenceVsManyAnalysis",
    "StructuralTransform",
    "TargetAnalysis",
]
