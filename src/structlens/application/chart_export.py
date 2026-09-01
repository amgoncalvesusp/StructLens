"""Publication-oriented exports for chart datasets."""

from __future__ import annotations

import math
from io import BytesIO
from pathlib import Path

from openpyxl import Workbook

from structlens.core.reports.safe_io import atomic_write_bytes

from .chart_data import ChartDataset, MatrixDataset

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
    else:
        sheet.append(["row", "column", "value", "text", "status"])
        for cell in dataset.cells:
            sheet.append([cell.row, cell.column, cell.value, cell.text, cell.status])
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
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as error:  # pragma: no cover - environment dependent
        raise RuntimeError("Chart image export requires the optional 'charts' dependency (matplotlib).") from error
    figure, axes = plt.subplots(figsize=(width_inches, height_inches), dpi=dpi)
    try:
        figure.patch.set_facecolor("white" if white_background else "#111827")
        axes.set_facecolor("white" if white_background else "#111827")
        if isinstance(dataset, ChartDataset):
            y_label = _unit_label(dataset.y_label, dataset.unit)
            for series in dataset.series:
                points = [(x, y) for x, y in series.points if y is not None]
                if points:
                    axes.plot(
                        [point[0] for point in points],
                        [point[1] for point in points],
                        marker="o",
                        label=series.name,
                    )
            axes.set_xlabel(dataset.x_label)
            axes.set_ylabel(y_label)
            axes.set_title(dataset.title)
            if len(dataset.series) > 1:
                axes.legend()
        else:
            rows = list(dict.fromkeys(cell.row for cell in dataset.cells))
            columns = list(dict.fromkeys(cell.column for cell in dataset.cells))
            values = {(cell.row, cell.column): cell.value for cell in dataset.cells}
            image = [[_matrix_value(values.get((row, column))) for column in columns] for row in rows]
            axes.imshow(image, aspect="auto", interpolation="nearest")
            axes.set_xticks(range(len(columns)), columns, rotation=45, ha="right")
            axes.set_yticks(range(len(rows)), rows)
            axes.set_xlabel(dataset.column_label)
            axes.set_ylabel(dataset.row_label)
            axes.set_title(dataset.title)
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
        plt.close(figure)


__all__ = ["export_chart_image", "export_chart_xlsx"]


def _unit_label(label: str, unit: str | None) -> str:
    if unit is None or unit in label:
        return label
    return f"{label} ({unit})"


def _matrix_value(value: float | None) -> float:
    return float("nan") if value is None else float(value)
