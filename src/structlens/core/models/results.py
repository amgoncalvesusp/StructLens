"""Immutable aggregate returned by one reference-versus-target analysis."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType

from structlens.core.models.correspondence import ResidueCorrespondence
from structlens.core.models.mutation import MutationEvent


@dataclass(frozen=True, slots=True)
class AnalysisResult:
    reference_id: str
    target_id: str
    correspondences: tuple[ResidueCorrespondence, ...]
    mutations: tuple[MutationEvent, ...]
    sequence_identity: float
    sequence_coverage: float
    alignment_decision: str
    strict_rmsd_angstrom: float | None = None
    refined_rmsd_angstrom: float | None = None
    mapped_residue_count: int = 0
    refined_residue_count: int | None = None
    excluded_alignment_indices: tuple[int, ...] = ()
    tm_score: float | None = None
    provenance: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "provenance", MappingProxyType(dict(self.provenance)))

    @property
    def mutation_count(self) -> int:
        return sum(event.kind.value != "conserved" for event in self.mutations)

    @property
    def insertion_count(self) -> int:
        return sum(event.kind.value == "insertion" for event in self.mutations)

    @property
    def deletion_count(self) -> int:
        return sum(event.kind.value == "deletion" for event in self.mutations)


__all__ = ["AnalysisResult"]
