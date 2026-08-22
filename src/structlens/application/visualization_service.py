"""Application facade that keeps visualization state separate from analysis."""

from structlens.core.models import AnalysisResult, ResidueCorrespondence
from structlens.plugin.visualization.renderer import (
    VisualizationRenderer,
    VisualizationState,
)


class VisualizationService:
    def __init__(self, renderer: VisualizationRenderer | None = None) -> None:
        self.renderer = renderer or VisualizationRenderer()

    def select(
        self, result: AnalysisResult, state: VisualizationState
    ) -> tuple[ResidueCorrespondence, ...]:
        return self.renderer.filtered_correspondences(result.correspondences, state)


__all__ = ["VisualizationService"]
