"""Report-driven facade that keeps visualization state separate from science."""

from typing import cast, overload

from structlens.core.models import AnalysisResult, ResidueCorrespondence
from structlens.core.reports import AnalysisReport, CorrespondenceSnapshot
from structlens.plugin.visualization.renderer import (
    VisualizationRenderer,
    VisualizationState,
)


class VisualizationService:
    def __init__(self, renderer: VisualizationRenderer | None = None) -> None:
        self.renderer = renderer or VisualizationRenderer()

    @overload
    def select(
        self, result: AnalysisReport, state: VisualizationState
    ) -> tuple[CorrespondenceSnapshot, ...]: ...

    @overload
    def select(
        self, result: AnalysisResult, state: VisualizationState
    ) -> tuple[ResidueCorrespondence, ...]: ...

    def select(
        self, result: AnalysisReport | AnalysisResult, state: VisualizationState
    ) -> tuple[CorrespondenceSnapshot, ...] | tuple[ResidueCorrespondence, ...]:
        if isinstance(result, AnalysisReport):
            if result.analysis is None:
                return ()
            return cast(
                tuple[CorrespondenceSnapshot, ...],
                self.renderer.filtered_correspondences(result.analysis.correspondences, state),
            )
        return self.renderer.filtered_correspondences(result.correspondences, state)


__all__ = ["VisualizationService"]
