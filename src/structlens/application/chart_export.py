"""Publication-oriented exports for chart datasets."""

from __future__ import annotations

import math
from io import BytesIO
from pathlib import Path

from openpyxl import Workbook

from structlens.core.reports.safe_io import atomic_write_bytes

from .chart_data import ChartDataset, DistanceChartCells, MatrixDataset
from .chart_plot import draw_chart

_MAX_CHART_BYTES = 100 * 1024 * 1024
_IMAGE_FORMATS = {".jpeg": "jpeg", ".jpg": "jpeg", ".png": "png", ".tif": "tiff", ".tiff": "tiff"}


def export_chart_xlsx(dataset: ChartDataset | MatrixDataset, path: str | Path) -> None:
    """Write the plotted values and labels as real numeric XLSX cells."""

    workbook = Workbook()
    sheet = workbook.active
    if sheet is None:
        raise RuntimeError("Could not create an XLSX worksheet")
    sheet.title = dataset.chart_id[:31]
    if isinstance(dataset, ChartDataset):
        sheet.append(["series", dataset.x_label, _unit_label(dataset.y_label, dataset.unit), "label"])
        for series in dataset.series:
            for index, (x, y) in enumerate(series.points):
                label = series.labels[index] if index < len(series.labels) else ""
                sheet.append([series.name, x, y, label])
    elif isinstance(dataset.cells, DistanceChartCells):
        matrix = dataset.cells.matrix
        if len(matrix.reference_positions) >= 16_384:
            raise ValueError("Distance map exceeds the XLSX column limit")
        sheet.append(["Target − reference (Å)", *matrix.reference_positions])
        for label, values, mask in zip(matrix.reference_positions, matrix.delta_angstrom, matrix.valid_mask, strict=True):
            sheet.append([label, *(v if valid else None for v, valid in zip(values, mask, strict=True))])
    else:
        sheet.append(["row", "column", "Target − reference (Å)" if dataset.chart_id == "distance_difference_heatmap" else "value", "text", "status"])
        for cell in dataset.cells:
            sheet.append([cell.row, cell.column, cell.value, cell.text, cell.status])
    # Labels originate in structures; workbook cells must remain literal text.
    for row in sheet.iter_rows():
        for output_cell in row:
            if output_cell.data_type == "f":
                output_cell.data_type = "s"
    sheet.freeze_panes = "A2"
    output = BytesIO()
    workbook.save(output)
    atomic_write_bytes(
        path,
        output.getvalue(),
        max_bytes=_MAX_CHART_BYTES,
        label="chart XLSX",
        replace_existing=True,
    )


def export_chart_image(
    dataset: ChartDataset | MatrixDataset,
    path: str | Path,
    *,
    dpi: int = 300,
    width_inches: float = 7.0,
    height_inches: float = 4.5,
    white_background: bool = True,
) -> None:
    """Render a labelled chart at a real publication DPI, never by upscaling."""

    if dpi not in {300, 600}:
        raise ValueError("Publication chart export supports 300 or 600 dpi")
    if (
        not math.isfinite(width_inches)
        or not math.isfinite(height_inches)
        or not 0.0 < width_inches <= 20.0
        or not 0.0 < height_inches <= 20.0
    ):
        raise ValueError("Chart dimensions must be finite, positive, and at most 20 inches")
    suffix = Path(path).suffix.casefold()
    image_format = _IMAGE_FORMATS.get(suffix)
    if image_format is None:
        raise ValueError("Chart image path must use JPEG, PNG, or TIFF")

    try:
        from matplotlib.figure import Figure
    except ImportError as error:  # pragma: no cover - environment dependent
        raise RuntimeError("Chart image export requires the optional 'charts' dependency (matplotlib).") from error
    figure = Figure(figsize=(width_inches, height_inches), dpi=dpi)
    try:
        figure.patch.set_facecolor("white" if white_background else "#111827")
        axes = figure.add_subplot(111)
        axes.set_facecolor("white" if white_background else "#111827")
        draw_chart(axes, dataset)
        figure.tight_layout()
        output = BytesIO()
        figure.savefig(output, format=image_format, dpi=dpi, facecolor=figure.get_facecolor())
        atomic_write_bytes(
            path,
            output.getvalue(),
            max_bytes=_MAX_CHART_BYTES,
            label="chart image",
            replace_existing=True,
        )
    finally:
        figure.clear()


__all__ = ["export_chart_image", "export_chart_xlsx"]


def _unit_label(label: str, unit: str | None) -> str:
    if unit is None or unit in label:
        return label
    return f"{label} ({unit})"
