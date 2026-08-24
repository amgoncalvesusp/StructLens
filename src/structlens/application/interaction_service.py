"""Application facade for interaction detection and reference comparison."""

from __future__ import annotations

from collections.abc import Callable, Iterable

from structlens.core.interactions import (
    InteractionDifference,
    InteractionRecord,
    InteractionThresholds,
    compare_interactions,
    detect_interactions,
)
from structlens.core.models import ResidueId, ResidueRecord


class InteractionAnalysisService:
    def detect(
        self,
        residues: Iterable[ResidueRecord],
        thresholds: InteractionThresholds | None = None,
        *,
        structure_id: str | None = None,
    ) -> tuple[InteractionRecord, ...]:
        return detect_interactions(residues, thresholds, structure_id=structure_id)

    def compare(
        self,
        reference: Iterable[InteractionRecord],
        target: Iterable[InteractionRecord],
        position_for: Callable[[ResidueId], str | None],
    ) -> tuple[InteractionDifference, ...]:
        return compare_interactions(reference, target, position_for)


__all__ = ["InteractionAnalysisService"]
