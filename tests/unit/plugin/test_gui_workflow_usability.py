"""User-visible readiness, pending settings and result navigation."""

import os
from concurrent.futures import Future
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from PySide6 import QtWidgets  # noqa: E402

from structlens.application.report_service import ReportService  # noqa: E402
from structlens.plugin.gui.main_panel import SCIENTIFIC_SECTIONS, build_qt_panel  # noqa: E402


@pytest.fixture
def controller():
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    panel = build_qt_panel(command=None)
    control = panel._structlens_controller
    yield control
    control.close()
    panel.deleteLater()
    app.processEvents()


def load_pair(control):
    for role in ("reference", "target"):
        assert control._load_source(role, Path("tests/fixtures/parsing/enriched.pdb"))
    control.mode_combo.setCurrentText("Sequence")


def finish_report(control):
    manual_pairs = ()
    if control.mode_combo.currentData() == "manual":
        manual_pairs = tuple(control._manual_pairs(
            control._selected_chain(control.reference_structure, control.reference_chain_combo, control.reference_model_combo),
            control._selected_chain(control.target_structure, control.target_chain_combo, control.target_model_combo),
        ))
    request = control._report_controller.build_request_from_snapshots(
        control.reference_loaded_source.snapshot,
        control.target_loaded_source.snapshot,
        reference_selection=control.reference_loaded_source.selection(model_id="1", chain_id="A"),
        target_selection=control.target_loaded_source.selection(model_id="1", chain_id="A"),
        analysis_settings=control._settings(),
        site_definitions=control._site_definitions,
        manual_pairs=manual_pairs,
    )
    report = ReportService().analyze(request)
    control._populate_report(report, request=request)
    return report


def test_navigation_follows_preparation_result_details_and_output(controller):
    assert SCIENTIFIC_SECTIONS == (
        "Project", "Structures", "Results", "Sequences", "Residues", "Sites", "Charts", "PyMOL", "Export"
    )
    for index, name in enumerate(SCIENTIFIC_SECTIONS):
        controller.nav.setCurrentRow(index)
        assert controller.pages.currentWidget().objectName() == f"page{name}"


def test_compare_requires_loaded_sources_and_recovers_after_path_edit(controller):
    assert not controller.compare_button.isEnabled()
    load_pair(controller)
    assert controller.compare_button.isEnabled()
    assert "enriched.pdb" in controller.workflow_context.text()
    assert "Sequence" in controller.workflow_context.text()
    previous = controller.reference_edit.text()
    controller.reference_edit.setText(previous + ".edited")
    assert not controller.compare_button.isEnabled()
    assert controller.load_sources_button.text() == "Load edited sources"
    controller.reference_edit.setText(previous)
    assert controller.compare_button.isEnabled()


def test_completed_report_opens_summary_and_settings_are_visibly_pending(controller):
    load_pair(controller)
    report = finish_report(controller)
    assert controller.pages.currentWidget().objectName() == "pageResults"
    assert controller.result_labels["sequence_identity"].text() == "100.0%"
    assert controller.result_labels["strict_rmsd_angstrom"].text() == "0.000"
    assert controller.header_status.text() == "Complete"
    assert report.report_id[:12] in controller.export_state_label.text()
    previous = controller.cutoff_spin.value()
    controller.cutoff_spin.setValue(previous + 1)
    assert controller.header_status.text() == "Changes pending"
    assert "last completed" in controller.export_state_label.text().lower()
    assert controller.model.report is report
    assert controller._canonical_binding().report is report
    controller.cutoff_spin.setValue(previous)
    assert controller.header_status.text() == "Complete"


def test_export_error_cannot_unlock_an_active_comparison(controller):
    load_pair(controller)
    future = Future()
    future.set_running_or_notify_cancel()
    controller._future = future
    controller._set_busy(True)
    try:
        controller._show_error("Export failed: output is not writable")
        assert controller.model.busy
        assert not controller.compare_button.isEnabled()
        assert not controller.reference_edit.isEnabled()
        controller._start_analysis()
        assert controller._future is future
    finally:
        future.set_result(None)
        controller._future = None


def test_manual_project_restores_pairs_before_repeating_comparison(controller, monkeypatch, tmp_path):
    load_pair(controller)
    controller.mode_combo.setCurrentIndex(3)
    controller.manual_edit.setPlainText("A:100 -> A:100\nA:100A -> A:100A")
    report = finish_report(controller)
    path = tmp_path / "manual.json"
    monkeypatch.setattr(QtWidgets.QFileDialog, "getSaveFileName", lambda *_args: (str(path), ""))
    controller._save_project()
    assert path.exists()
    controller.manual_edit.clear()
    monkeypatch.setattr(QtWidgets.QFileDialog, "getOpenFileName", lambda *_args: (str(path), ""))
    controller._open_project()
    assert controller.manual_edit.toPlainText() == "A:100 -> A:100\nA:100A -> A:100A"
    assert controller.header_status.text() == "Complete"
    assert finish_report(controller).report_id == report.report_id


def test_site_pending_and_busy_state_preserve_completed_report(controller):
    load_pair(controller)
    report = finish_report(controller)
    residue = controller.reference_structure.chains[0].residue_records[0].residue_id
    controller.site_residues_edit.setText(f"{residue.chain_id}:{residue.auth_seq_id}")
    controller._define_site_from_controls()
    assert "Compare structures" in controller.site_status_label.text()
    assert controller.header_status.text() == "Changes pending"
    controller._set_busy(True)
    assert not controller.compare_button.isEnabled()
    assert not controller.reference_edit.isEnabled()
    controller._analysis_cancelled()
    assert controller.compare_button.isEnabled()
    assert controller.reference_edit.isEnabled()
    assert controller.model.report is report
    assert controller.header_status.text() == "Changes pending"


def test_pymol_report_limitation_is_visible_before_clicking(controller):
    load_pair(controller)
    finish_report(controller)
    assert not controller.pymol_open_button.isEnabled()
    assert not controller.pymol_export_button.isEnabled()
    assert "unavailable" in controller.pymol_report_status.text().lower()
