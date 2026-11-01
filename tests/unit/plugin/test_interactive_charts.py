"""Scientific chart navigation uses only measured report positions."""

import os
from dataclasses import replace

import numpy as np
import pytest
from openpyxl import load_workbook

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from PySide6 import QtWidgets  # noqa: E402

from structlens.application.chart_export import export_chart_xlsx  # noqa: E402
from structlens.plugin.gui.main_panel import build_qt_panel  # noqa: E402
from tests.unit.plugin._task13_round4_fixtures import rich_report  # noqa: E402
from tests.unit.plugin.test_gui_workflow_usability import finish_report, load_pair  # noqa: E402

MouseEvent = pytest.importorskip("matplotlib.backend_bases").MouseEvent


def test_identity_matrix_uses_canonical_events_and_masks_nonstandard_positions(control, tmp_path):
    load_pair(control)
    report = finish_report(control)
    assert report.analysis.mutation_count == 0
    control.chart_combo.setCurrentText("Mutation / conservation matrix")
    matrix = control._chart_datasets["Mutation / conservation matrix"]
    cells = {cell.column: cell for cell in matrix.cells}
    assert cells["A:101"].value is None and cells["A:101"].status == "nonstandard"
    assert cells["A:100A"].value == 1.0 and cells["A:100A"].text == "G100AG"
    mesh = control._chart_canvas.figure.axes[0].images[0]
    assert mesh.get_clim() == (0.0, 1.0)
    assert mesh.get_array().mask[0, 2]
    out = tmp_path / "mutations.xlsx"
    export_chart_xlsx(matrix, out)
    rows = {row[1]: row for row in load_workbook(out).active.values}
    assert rows["A:101"][2] is None and rows["A:101"][4] == "nonstandard"
    assert rows["A:100A"][3] == "G100AG"


@pytest.fixture
def control():
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    panel = build_qt_panel()
    yield panel._structlens_controller
    panel._structlens_controller.close()
    panel.deleteLater()
    app.processEvents()


def test_distance_chart_masks_missing_cells_and_selects_both_residues(control, tmp_path):
    request, report = rich_report(tmp_path)
    distance = report.distance_map
    labels = distance.reference_positions[:2]
    assert len(labels) == 2
    report = replace(
        report,
        distance_map=replace(
            distance,
            reference_positions=labels,
            reference_distances_angstrom=((0.0, 4.0), (4.0, 0.0)),
            target_distances_angstrom=((0.0, 3.0), (3.0, 0.0)),
            delta_angstrom=((0.0, -1.0), (-1.0, 0.0)),
            valid_mask=((False, True), (True, True)),
        ),
    )
    control._populate_report(report, request=request)
    control.chart_combo.setCurrentText("Distance-difference heatmap")
    dataset = control._chart_datasets["Distance-difference heatmap"]
    canvas = control._chart_canvas
    axes = canvas.figure.axes[0]
    img = axes.images[0]
    assert img.get_clim() == (-1.0, 1.0)
    assert img.get_array().mask[0, 0]
    output = tmp_path / "map.xlsx"
    export_chart_xlsx(dataset, output)
    rows = list(load_workbook(output).active.values)
    assert rows[1][1] is None and rows[1][2] == -1.0
    canvas.draw()
    x, y = axes.transData.transform((1, 0))
    canvas.callbacks.process("button_press_event", MouseEvent("button_press_event", canvas, x, y, button=1))
    assert all(label in control.chart_selection_label.text() for label in labels)
    assert "-1.000" in control.chart_selection_label.text()
    for button, label in zip(control.chart_residue_buttons, labels, strict=True):
        button.click()
        assert control.pages.currentWidget().objectName() == "pageResidues"
        assert label in control.evidence_card_label.text()


def test_missing_points_break_profile_and_click_opens_reported_alignment_index(control, tmp_path):
    request, report = rich_report(tmp_path)
    corr = report.analysis.correspondences
    assert len(corr) >= 3
    report = replace(
        report,
        analysis=replace(
            report.analysis,
            correspondences=(
                replace(corr[0], ca_displacement_angstrom=0.2),
                replace(corr[1], ca_displacement_angstrom=None),
                replace(corr[2], ca_displacement_angstrom=0.4),
                *corr[3:],
            ),
        ),
    )
    control._populate_report(report, request=request)
    canvas = control._chart_canvas
    line = canvas.figure.axes[0].lines[0]
    assert np.isnan(line.get_xydata()[1, 1])
    canvas.draw()
    x, y = line.get_transform().transform(line.get_xydata()[2])
    canvas.callbacks.process("button_press_event", MouseEvent("button_press_event", canvas, x, y, button=1))
    assert control.pages.currentWidget().objectName() == "pageResidues"
    assert control.residue_table.item(control.residue_table.currentRow(), 0).text() == str(corr[2].alignment_index)


def test_chart_unavailability_offers_reason_and_disables_export(control):
    assert not any(button.isEnabled() for button in control.chart_export_buttons)
    control.chart_combo.setCurrentText("Pairwise similarity heatmap")
    assert "target" in control.chart_preview_status.text().lower()
    assert not any(button.isEnabled() for button in control.chart_export_buttons)


def test_availability_never_scans_lazy_distance_cells_on_selection(control, tmp_path, monkeypatch):
    request, report = rich_report(tmp_path)
    control._populate_report(report, request=request)
    name = "Distance-difference heatmap"
    data = control._chart_datasets[name]

    def forbidden(*_args):
        pytest.fail("Checking availability materialized a distance cell")

    monkeypatch.setattr(type(data.cells), "__getitem__", forbidden)
    for _ in range(25):
        assert control._available_chart(name) is data


def test_msa_navigation_resolves_reference_identity_even_when_indices_differ(control, tmp_path):
    request, report = rich_report(tmp_path)
    analysis = report.analysis
    report = replace(
        report,
        analysis=replace(
            analysis,
            correspondences=tuple(replace(c, alignment_index=c.alignment_index + 30) for c in analysis.correspondences),
        ),
    )
    control._populate_report(report, request=request)
    control.chart_combo.setCurrentText("MSA conservation profile")
    series = control._chart_datasets["MSA conservation profile"].series[0]
    visible = next(
        i
        for i, (_, y) in enumerate(series.points)
        if y is not None and series.metadata[i]["alignment_index"] is not None
    )
    index = series.metadata[visible]["alignment_index"]
    canvas = control._chart_canvas
    canvas.draw()
    line = canvas.figure.axes[0].lines[0]
    x, y = line.get_transform().transform(line.get_xydata()[visible])
    canvas.callbacks.process("button_press_event", MouseEvent("button_press_event", canvas, x, y, button=1))
    assert control.pages.currentWidget().objectName() == "pageResidues"
    assert int(control.residue_table.item(control.residue_table.currentRow(), 0).text()) == index
    assert index >= 30
