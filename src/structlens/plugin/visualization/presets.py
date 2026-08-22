"""Named visualization presets shared by headless and Qt surfaces."""

from .renderer import VisualizationRenderer, VisualizationState


def preset_state(name: str) -> VisualizationState:
    return VisualizationRenderer().apply_preset(name)


__all__ = ["preset_state"]
