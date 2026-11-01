"""Optional Qt chart backend loading kept outside the panel composition."""

from __future__ import annotations

import importlib
from typing import Any

from structlens.application.chart_plot import matrix_image_kwargs


def load_chart_classes() -> tuple[Any, Any] | None:
    """Return the Qt canvas and Matplotlib figure classes when installed."""

    try:
        canvas_module = importlib.import_module("matplotlib.backends.backend_qtagg")
        figure_module = importlib.import_module("matplotlib.figure")
    except ImportError:
        return None
    return getattr(canvas_module, "FigureCanvasQTAgg"), getattr(figure_module, "Figure")  # noqa: B009


__all__ = ["load_chart_classes", "matrix_image_kwargs"]
