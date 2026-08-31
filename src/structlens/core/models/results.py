"""Immutable aggregate returned by one reference-versus-target analysis."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import TYPE_CHECKING

from structlens.core.models.correspondence import ResidueCorrespondence
from structlens.core.models.mutation import MutationEvent
from structlens.core.provenance import MethodProvenance

if TYPE_CHECKING:
    from .multi import StructuralTransform


@dataclass(frozen=True, slots=True)
class AnalysisResult:
    reference_id: str
    target_id: str
    correspondences: tuple[ResidueCorrespondence, ...]
    mutations: tuple[MutationEvent, ...]
    sequence_identity: float
    sequence_coverage: float
    alignment_decision: str
    sequence_similarity: float | None = None
    strict_rmsd_angstrom: float | None = None
    refined_rmsd_angstrom: float | None = None
    mapped_residue_count: int = 0
    refined_residue_count: int | None = None
    excluded_alignment_indices: tuple[int, ...] = ()
    tm_score: float | None = None
    # Keep the v0.3 flat mapping separate from the typed scientific contract.
    provenance: Mapping[str, str] = field(default_factory=dict)
    transform: StructuralTransform | None = None
    method_provenance: MethodProvenance | None = None

    def __post_init__(self) -> None:
        if isinstance(self.provenance, MethodProvenance):
            raise TypeError("pass typed provenance via method_provenance, not legacy provenance")
        object.__setattr__(self, "provenance", MappingProxyType(dict(self.provenance)))

    @property
    def typed_provenance(self) -> MethodProvenance | None:
        """Compatibility alias for consumers migrating to typed provenance."""

        return self.method_provenance

    @property
    def mutation_count(self) -> int:
        """Count sequence changes only.

        A non-standard residue such as MSE has no canonical one-letter code, so
        it is neither conserved nor a substitution. Counting it as a mutation
        made a structure compared against itself report a change it does not
        have. Those positions stay in ``mutations`` as non-standard evidence.
        """

        return sum(event.kind.value not in {"conserved", "nonstandard"} for event in self.mutations)

    @property
    def insertion_count(self) -> int:
        return sum(event.kind.value == "insertion" for event in self.mutations)

    @property
    def deletion_count(self) -> int:
        return sum(event.kind.value == "deletion" for event in self.mutations)


__all__ = ["AnalysisResult"]
