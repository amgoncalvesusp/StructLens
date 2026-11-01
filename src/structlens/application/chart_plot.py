"""Render the same authoritative data in a Qt preview and publication images."""

from collections.abc import Iterable
from typing import Any

from .chart_data import ChartDataset, DistanceChartCells, MatrixDataset


def draw_chart(axes: Any, dataset: ChartDataset | MatrixDataset) -> None:
    if isinstance(dataset, ChartDataset):
        for series in dataset.series:
            if not series.points:
                continue
            axes.plot(
                [x for x, _ in series.points],
                [float("nan") if y is None else y for _, y in series.points],
                marker="o",
                linewidth=1.4,
                label=series.name,
                picker=5,
            )
        axes.set_xlabel(dataset.x_label)
        axes.set_ylabel(
            dataset.y_label
            if not dataset.unit or dataset.unit in dataset.y_label
            else f"{dataset.y_label} ({dataset.unit})"
        )
        if len(dataset.series) > 1 and axes.lines:
            axes.legend()
    elif isinstance(dataset.cells, DistanceChartCells):
        import numpy as np
        from matplotlib import colormaps

        matrix = dataset.cells.matrix
        positions = matrix.reference_positions
        masked: np.ma.MaskedArray[Any, Any] = np.ma.array(matrix.delta_angstrom, mask=np.logical_not(matrix.valid_mask))
        bound = float(np.ma.max(np.ma.abs(masked))) if masked.count() else 0.0
        bound = bound or 1.0
        mesh = axes.imshow(
            masked,
            aspect="equal",
            interpolation="nearest",
            vmin=-bound,
            vmax=bound,
            cmap=colormaps["RdBu_r"].with_extremes(bad="#aeb6bf"),
        )
        axes.figure.colorbar(mesh, ax=axes, label="Target − reference (Å)")
        step = max(1, (len(positions) + 13) // 14)
        axes.set_xticks(range(0, len(positions), step), positions[::step], rotation=45, ha="right")
        axes.set_yticks(range(0, len(positions), step), positions[::step])
        axes.set_xlabel(dataset.column_label)
        axes.set_ylabel(dataset.row_label)
    else:
        rows = list(dict.fromkeys(cell.row for cell in dataset.cells))
        columns = list(dict.fromkeys(cell.column for cell in dataset.cells))
        values = {(cell.row, cell.column): cell.value for cell in dataset.cells}
        if rows and columns:
            image = [[float("nan") if (value := values.get((r, c))) is None else value for c in columns] for r in rows]
            kwargs = matrix_image_kwargs(dataset, values.values())
            if dataset.chart_id == "mutation_conservation_matrix":
                from matplotlib.colors import ListedColormap

                palette = ListedColormap(["#d07922", "#3977a4"]).with_extremes(bad="#aeb6bf")
                mesh = axes.imshow(image, aspect="auto", interpolation="nearest", cmap=palette, vmin=0.0, vmax=1.0)
                legend = axes.figure.colorbar(mesh, ax=axes, ticks=[0, 1])
                legend.set_ticks([0, 1], labels=["Substitution", "Conserved"])
            else:
                axes.imshow(image, **kwargs)
            for labels, set_ticks, args in (
                (rows, axes.set_yticks, {}),
                (columns, axes.set_xticks, {"rotation": 45, "ha": "right"}),
            ):
                step = max(1, (len(labels) + 13) // 14)
                set_ticks(range(0, len(labels), step), labels[::step], **args)
        axes.set_xlabel(dataset.column_label)
        axes.set_ylabel(dataset.row_label)
    axes.set_title(dataset.title)


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
