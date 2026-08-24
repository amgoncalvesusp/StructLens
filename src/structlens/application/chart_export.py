"""Publication-oriented exports for chart datasets."""

from __future__ import annotations

from pathlib import Path

from openpyxl import Workbook

from .chart_data import ChartDataset, MatrixDataset


def export_chart_xlsx(dataset: ChartDataset | MatrixDataset, path: str | Path) -> None:
    """Write the plotted values and labels as real numeric XLSX cells."""

    workbook = Workbook()
    sheet = workbook.active
    if sheet is None:
        raise RuntimeError("Could not create an XLSX worksheet")
    sheet.title = dataset.chart_id[:31]
    if isinstance(dataset, ChartDataset):
        sheet.append(["series", dataset.x_label, dataset.y_label, "label"])
        for series in dataset.series:
            for index, (x, y) in enumerate(series.points):
                label = series.labels[index] if index < len(series.labels) else ""
                sheet.append([series.name, x, y, label])
    else:
        sheet.append(["row", "column", "value", "text", "status"])
        for cell in dataset.cells:
            sheet.append([cell.row, cell.column, cell.value, cell.text, cell.status])
    sheet.freeze_panes = "A2"
    workbook.save(path)


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

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as error:  # pragma: no cover - environment dependent
        raise RuntimeError(
            "Chart image export requires the optional 'charts' dependency (matplotlib)."
        ) from error
    figure, axes = plt.subplots(figsize=(width_inches, height_inches), dpi=dpi)
    figure.patch.set_facecolor("white" if white_background else "#111827")
    axes.set_facecolor("white" if white_background else "#111827")
    if isinstance(dataset, ChartDataset):
        for series in dataset.series:
            points = [(x, y) for x, y in series.points if y is not None]
            if points:
                axes.plot([point[0] for point in points], [point[1] for point in points], marker="o", label=series.name)
        axes.set_xlabel(dataset.x_label)
        axes.set_ylabel(dataset.y_label)
        axes.set_title(dataset.title)
        if len(dataset.series) > 1:
            axes.legend()
    else:
        rows = list(dict.fromkeys(cell.row for cell in dataset.cells))
        columns = list(dict.fromkeys(cell.column for cell in dataset.cells))
        values = {
            (cell.row, cell.column): cell.value
            for cell in dataset.cells
        }
        image: list[list[float]] = []
        for row in rows:
            image_row: list[float] = []
            for column in columns:
                value = values.get((row, column))
                image_row.append(float("nan") if value is None else value)
            image.append(image_row)
        axes.imshow(image, aspect="auto", interpolation="nearest")
        axes.set_xticks(range(len(columns)), columns, rotation=45, ha="right")
        axes.set_yticks(range(len(rows)), rows)
        axes.set_xlabel(dataset.column_label)
        axes.set_ylabel(dataset.row_label)
        axes.set_title(dataset.title)
    figure.tight_layout()
    figure.savefig(path, dpi=dpi, facecolor=figure.get_facecolor())
    plt.close(figure)


__all__ = ["export_chart_image", "export_chart_xlsx"]
