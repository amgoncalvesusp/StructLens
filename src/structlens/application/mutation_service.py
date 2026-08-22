"""Application facade for mutation descriptors."""

from collections.abc import Sequence

from structlens.core.models import MutationEvent, ResidueCorrespondence
from structlens.core.mutations.detector import detect_mutations


class MutationService:
    def detect(self, correspondences: Sequence[ResidueCorrespondence]) -> tuple[MutationEvent, ...]:
        return tuple(detect_mutations(correspondences))


__all__ = ["MutationService"]
