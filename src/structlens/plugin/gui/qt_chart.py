"""Optional Qt chart backend loading kept outside the panel composition."""

from __future__ import annotations

import importlib
from collections.abc import Iterable
from typing import Any

from structlens.application.chart_data import MatrixDataset


def matrix_image_kwargs(dataset: MatrixDataset, values: Iterable[float | None]) -> dict[str, Any]:
    """Return display-only limits without changing authoritative matrix values."""

    kwargs: dict[str, Any] = {"aspect": "auto", "interpolation": "nearest"}
    numeric_values = tuple(float(value) for value in values if value is not None)
    if not numeric_values:
        return kwargs

    minimum = min(numeric_values)
    maximum = max(numeric_values)
    descriptor = " ".join(
        (
            dataset.chart_id,
            dataset.title,
            dataset.row_label,
            dataset.column_label,
            dataset.interpretation,
        )
    ).casefold()
    is_delta = (
        "delta" in descriptor
        or "Δ" in descriptor
        or "distance difference" in descriptor
        or "distance_difference" in descriptor
    )
    if is_delta:
        bound = max(abs(minimum), abs(maximum))
        if bound > 0.0:
            kwargs.update(vmin=-bound, vmax=bound)
    elif minimum < maximum:
        kwargs.update(vmin=minimum, vmax=maximum)
    return kwargs


def load_chart_classes() -> tuple[Any, Any] | None:
    """Return the Qt canvas and Matplotlib figure classes when installed."""

    try:
        canvas_module = importlib.import_module("matplotlib.backends.backend_qtagg")
        figure_module = importlib.import_module("matplotlib.figure")
    except ImportError:
        return None
    return getattr(canvas_module, "FigureCanvasQTAgg"), getattr(figure_module, "Figure")  # noqa: B009


__all__ = ["load_chart_classes", "matrix_image_kwargs"]
