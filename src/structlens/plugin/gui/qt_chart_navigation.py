"""Availability and residue navigation for measured report charts."""

from typing import Any

from structlens.application.chart_data import (
    ChartDataset,
    DistanceChartCells,
    MatrixDataset,
    structural_deviation_profile,
)
from structlens.application.report_geometry import reference_position
from structlens.core.models import AnalysisResult

from .qt_context import QtMixinContext

_UNAVAILABLE = {
    "Structural deviation profile": "Compare structures with mapped Cα coordinates.",
    "MSA conservation profile": "Compare sequences with at least one measurable conservation column.",
    "Distance-difference heatmap": "Compare structures with measured internal Cα distances.",
    "Mutation / conservation matrix": "Requires at least one conserved or substituted canonical amino acid at a reference position.",
    "Sequence logo": "Inspect conservation on the MSA conservation profile. A sequence-logo dataset is not supplied by this report.",
    "Pairwise similarity heatmap": "Requires a matrix of comparisons across multiple targets.",
    "Sequence–structure relationship": "Requires multiple targets with both sequence and structure measurements.",
    "Structural conservation profile": "Requires several targets aligned to a shared reference.",
    "Key-residue comparison": "Requires a chart dataset for selected reference residues across targets.",
    "Interaction difference matrix": "Inspect measured interaction differences in the report details. This report supplies no interaction chart.",
    "Site comparison": "Inspect site measurements on Sites & pockets. This report supplies no site-comparison chart.",
    "Evidence Card figure": "Select a row on Residues to read its evidence card. This report supplies no card figure.",
}


class ChartNavigationMixin(QtMixinContext):
    def _available_chart(self, label: str, result: AnalysisResult | None = None) -> ChartDataset | MatrixDataset | None:
        dataset = self._chart_datasets.get(label)
        if dataset is None and label == "Structural deviation profile":
            legacy = result or self.model.analysis
            if legacy is not None:
                dataset = structural_deviation_profile(legacy)
        if dataset is None and label == "MSA conservation profile":
            dataset = self._msa_chart_dataset
        if isinstance(dataset, ChartDataset):
            return dataset if any(y is not None for series in dataset.series for _, y in series.points) else None
        if isinstance(dataset, MatrixDataset):
            if isinstance(dataset.cells, DistanceChartCells):
                return dataset if dataset.cells.has_measurements else None
            return dataset if any(cell.value is not None for cell in dataset.cells) else None
        return None

    def _chart_unavailable_reason(self, label: str) -> str:
        return f"The {label} dataset is unavailable. " + _UNAVAILABLE.get(label, "No measured data was loaded for this chart.")

    def _clear_chart_selection(self) -> None:
        self.chart_selection_label.clear()
        for button in self.chart_residue_buttons:
            button.hide()
            button.setProperty("alignment_index", None)

    def _open_chart_endpoint(self, button: Any) -> None:
        index = button.property("alignment_index")
        if isinstance(index, int):
            self._open_residue_filter("all", index)

    def _wire_chart_selection(self, canvas: Any, axes: Any, dataset: ChartDataset | MatrixDataset) -> None:
        if isinstance(dataset, ChartDataset):
            series = [s for s in dataset.series if s.points]

            def pick(event: Any) -> None:
                if canvas is not self._chart_canvas or event.artist not in axes.lines or event.mouseevent.button != 1:
                    return
                line = event.artist
                line_index = list(axes.lines).index(line)
                if line_index >= len(series) or len(event.ind) == 0:
                    return
                # Insertion codes can share the same numeric author position.
                point = min(
                    event.ind,
                    key=lambda i: sum(
                        (a - b) ** 2
                        for a, b in zip(
                            line.get_transform().transform(line.get_xydata()[i]),
                            (event.mouseevent.x, event.mouseevent.y),
                            strict=True,
                        )
                    ),
                )
                values = series[line_index]
                if point < len(values.metadata) and isinstance(
                    index := values.metadata[point].get("alignment_index"), int
                ):
                    self._open_residue_filter("all", index)

            canvas.mpl_connect("pick_event", pick)
        elif dataset.chart_id == "distance_difference_heatmap":

            def pick_cell(event: Any) -> None:
                if canvas is self._chart_canvas and event.inaxes is axes and event.button == 1:
                    self._select_distance_cell(event.xdata, event.ydata)

            canvas.mpl_connect("button_press_event", pick_cell)

    def _select_distance_cell(self, x: float, y: float) -> None:
        report = self.model.report
        if report is None or report.distance_map is None or report.analysis is None:
            return
        distance = report.distance_map
        labels = distance.reference_positions
        if not (-0.5 <= x < len(labels) - 0.5 and -0.5 <= y < len(labels) - 0.5):
            return
        row, column = int(y + 0.5), int(x + 0.5)
        pair = labels[row], labels[column]
        lookup = {
            reference_position(c.reference): c.alignment_index for c in report.analysis.correspondences if c.reference
        }
        if distance.valid_mask[row][column]:
            value = (
                f"reference {distance.reference_distances_angstrom[row][column]:.3f} Å · "
                f"target {distance.target_distances_angstrom[row][column]:.3f} Å · "
                f"Δ {distance.delta_angstrom[row][column]:+.3f} Å"
            )
        else:
            value = "unavailable (missing coordinates)"
        self.chart_selection_label.setText(f"{pair[0]} ↔ {pair[1]} · {value}")
        for button, label in zip(self.chart_residue_buttons, pair, strict=True):
            button.setText(f"Inspect {label}")
            button.setProperty("alignment_index", lookup.get(label))
            button.setEnabled(label in lookup)
            button.show()

    def _highlight_chart_residue(self, alignment_index: int | None) -> None:
        canvas = self._chart_canvas
        if canvas is None or not canvas.figure.axes:
            return
        for artist in getattr(canvas, "_structlens_highlights", ()):
            artist.remove()
        marks = []
        axes = canvas.figure.axes[0]
        dataset = self._available_chart(self.chart_combo.currentText())
        if isinstance(dataset, ChartDataset) and alignment_index is not None:
            for series in dataset.series:
                for point, meta in zip(series.points, series.metadata, strict=False):
                    if meta.get("alignment_index") == alignment_index and point[1] is not None:
                        marks.append(axes.scatter(*point, facecolor="none", edgecolor="#d87800", s=90, zorder=5))
        elif isinstance(dataset, MatrixDataset) and dataset.chart_id == "distance_difference_heatmap":
            report = self.model.report
            if report is not None and report.analysis is not None and report.distance_map is not None:
                item = next(
                    (item for item in report.analysis.correspondences if item.alignment_index == alignment_index), None
                )
                labels = report.distance_map.reference_positions
                if item is not None and item.reference is not None and reference_position(item.reference) in labels:
                    position = labels.index(reference_position(item.reference))
                    marks.extend(
                        (
                            axes.axvline(position, color="#333333", linestyle=":"),
                            axes.axhline(position, color="#333333", linestyle=":"),
                        )
                    )
        canvas._structlens_highlights = marks
        canvas.draw_idle()
