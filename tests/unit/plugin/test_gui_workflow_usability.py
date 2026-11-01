"""User-visible readiness, pending settings and result navigation."""

import os
from concurrent.futures import Future
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from PySide6 import QtCore, QtTest, QtWidgets  # noqa: E402

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


def test_fitted_report_enables_native_pymol_scene(controller):
    load_pair(controller)
    finish_report(controller)
    assert controller.pymol_open_button.isEnabled()
    assert controller.pymol_export_button.isEnabled()
    assert ".pml" in controller.pymol_report_status.text().lower()


def compare_changed_pair(controller, tmp_path, *, reverse=False, two_substitutions=False):
    # Different author numbers and an insertion code make row-index mixups visible.
    paths = []
    for name, start, residues in (
        ("reference", 101, ("ALA", "VAL", "SER", "LYS", "GLY", "PHE")),
        ("target", 201, ("ALA", "VAL", "THR", "LYS", "ALA" if two_substitutions else "GLY")),
    ):
        lines = [
            f"ATOM  {i + 1:5d}  CA  {residue:3} A{start + i:4d}{'B' if i == 2 else ' '}"
            f"   {3.8 * i:8.3f}{float(i % 2):8.3f}{0.0:8.3f}  1.00 15.00           C\n"
            for i, residue in enumerate(residues)
        ]
        path = tmp_path / f"{name}.pdb"
        path.write_text("".join(lines) + "END\n", encoding="ascii")
        paths.append(path)
    reference, target = reversed(paths) if reverse else paths
    for role, path in (("reference", reference), ("target", target)):
        assert controller._load_source(role, path)
    controller.mode_combo.setCurrentText("Sequence")
    request = controller._report_controller.build_request(
        reference, target, reference_chain_id="A", target_chain_id="A", analysis_settings=controller._settings()
    )
    report = ReportService().analyze(request)
    controller._populate_report(report, request=request)


def visible_residue_rows(control):
    table = control.residue_table
    return [i for i in range(table.rowCount()) if not table.isRowHidden(i)]


def test_author_search_selects_exact_model_position_and_keeps_gap(controller, tmp_path):
    compare_changed_pair(controller, tmp_path)
    controller.nav.setCurrentRow(SCIENTIFIC_SECTIONS.index("Residues"))
    controller.residue_search.setText("A:203B")
    assert visible_residue_rows(controller) == [2]
    assert "A:103B SER" in controller.evidence_card_label.text()
    controller.residue_search.setText("A:203")
    assert not visible_residue_rows(controller)
    controller.residue_search.setText("A:106")
    assert visible_residue_rows(controller) == [5]
    assert "Target gap" in controller.evidence_card_label.text()
    controller._open_residue_filter("substitution", 2)
    assert not controller.residue_search.text()
    assert visible_residue_rows(controller) == [2]


def test_alignment_click_uses_report_correspondence_and_tracks_focus(controller, tmp_path):
    compare_changed_pair(controller, tmp_path)
    strip = controller.alignment_strip
    assert strip.item(0, 2).text() == "S" and strip.item(1, 2).text() == "T"
    assert strip.item(1, 5).text() == "–"
    controller.nav.setCurrentRow(SCIENTIFIC_SECTIONS.index("Sequences"))
    strip.cellClicked.emit(1, 5)
    assert controller.pages.currentWidget().objectName() == "pageResidues"
    assert controller.residue_table.currentRow() == 5
    assert "Target gap" in controller.evidence_card_label.text()
    controller._open_residue_filter("substitution", 2)
    assert strip.currentColumn() == 2


def test_alignment_tooltips_quote_deposited_identifiers(controller, tmp_path):
    from dataclasses import replace

    compare_changed_pair(controller, tmp_path)
    row = controller.model.report.analysis.correspondences[0]
    unusual_chain = replace(row.reference, chain_id="<img src='file:///example.png'>")
    controller._set_alignment_rows((replace(row, reference=unusual_chain),))
    tooltip = controller.alignment_strip.item(0, 0).toolTip()
    assert "<img" not in tooltip
    assert "&lt;img" in tooltip


def test_alignment_window_only_marks_the_selected_report_column(controller, tmp_path):
    from dataclasses import replace

    compare_changed_pair(controller, tmp_path)
    row = controller.model.report.analysis.correspondences[0]
    controller._set_alignment_rows(tuple(replace(row, alignment_index=i) for i in range(30)))
    controller._sync_alignment_selection(5)
    assert controller.alignment_strip.currentColumn() == 5
    controller.alignment_window.setValue(2)
    assert not controller.alignment_strip.selectedItems()
    controller.alignment_window.setValue(1)
    assert controller.alignment_strip.currentColumn() == 5


def test_identical_chains_show_no_mutation_and_no_drilldown(controller):
    load_pair(controller)
    finish_report(controller)
    count = controller.result_labels["mutation_count"]
    assert count.text() == "0"
    assert not count.isEnabled()
    assert all(not button.isEnabled() for button in controller.mutation_kind_buttons.values())
    assert controller.results_table.item(0, 9).text() == "0"


def test_summary_click_filters_selects_and_updates_evidence_with_keyboard(controller, tmp_path):
    compare_changed_pair(controller, tmp_path)
    controller.widget.resize(1200, 780)
    controller.widget.show()
    count = controller.result_labels["mutation_count"]
    assert count.text() == "2"
    assert controller.results_table.item(0, 9).text() == "2"
    count.click()
    QtWidgets.QApplication.processEvents()
    assert controller.pages.currentWidget().objectName() == "pageResidues"
    assert controller.residue_filter_combo.currentData() == "changes"
    table = controller.residue_table
    change_rows = visible_residue_rows(controller)
    assert [table.item(row, 3).text() for row in change_rows] == ["substitution", "deletion"]
    assert table.currentRow() == change_rows[0]
    assert "A:103B" in controller.evidence_card_label.text()
    assert "A:203B" in controller.evidence_card_label.text()
    evidence_left = controller.evidence_card_label.mapTo(controller.widget, QtCore.QPoint(0, 0)).x()
    table_right = table.mapTo(controller.widget, QtCore.QPoint(table.width(), 0)).x()
    assert evidence_left >= table_right
    controller.widget.resize(800, 600)
    QtWidgets.QApplication.processEvents()
    page_scroll = controller.pages.currentWidget().findChild(QtWidgets.QScrollArea)
    assert page_scroll.horizontalScrollBar().maximum() == 0
    card = controller.evidence_card_label
    assert card.height() >= card.heightForWidth(card.width())

    table.setFocus()
    QtTest.QTest.keyClick(table, QtCore.Qt.Key.Key_Down)
    assert table.currentRow() == change_rows[1]
    assert "A:106" in controller.evidence_card_label.text()
    assert "A:203B" not in controller.evidence_card_label.text()

    combo = controller.residue_filter_combo
    combo.setCurrentIndex(combo.findData("all"))
    assert len(visible_residue_rows(controller)) == 6
    table.scrollToTop()
    first_cell = table.visualItemRect(table.item(0, 0)).center()
    QtTest.QTest.mouseClick(table.viewport(), QtCore.Qt.MouseButton.LeftButton, pos=first_cell)
    assert table.currentRow() == 0
    assert "A:101" in controller.evidence_card_label.text()


@pytest.mark.parametrize(
    ("reverse", "kind", "position"),
    [(False, "substitution", "A:103B"), (False, "deletion", "A:106"), (True, "insertion", "A:106")],
)
def test_mutation_kind_button_selects_exact_affected_position(controller, tmp_path, reverse, kind, position):
    compare_changed_pair(controller, tmp_path, reverse=reverse)
    button = controller.mutation_kind_buttons[kind]
    assert button.isEnabled()
    button.click()
    table = controller.residue_table
    assert controller.pages.currentWidget().objectName() == "pageResidues"
    assert controller.residue_filter_combo.currentData() == kind
    rows = visible_residue_rows(controller)
    assert len(rows) == 1
    assert rows == [table.currentRow()]
    assert table.item(rows[0], 3).text() == kind
    assert position in controller.evidence_card_label.text()
    if kind == "insertion":
        selection = controller._selection_controller
        assert selection.lookup_correspondence().target.auth_seq_id == "106"
        with pytest.raises(KeyError):
            selection.lookup_evidence_card()
        assert "Reference gap" in controller.evidence_card_label.text()
        assert "Evidence card unavailable" in controller.evidence_card_label.text()


def test_empty_residue_filter_clears_the_selected_report_position(controller, tmp_path):
    compare_changed_pair(controller, tmp_path)
    controller.mutation_kind_buttons["deletion"].click()
    assert controller.model.selection is not None
    combo = controller.residue_filter_combo
    combo.setCurrentIndex(combo.findData("insertion"))
    assert visible_residue_rows(controller) == []
    assert controller.residue_table.currentRow() == -1
    assert not controller.residue_table.selectedItems()
    assert "A:106" not in controller.evidence_card_label.text()
    assert controller.model.selection is None


def test_panel_transaction_restores_selected_change_and_filter(controller, tmp_path):
    compare_changed_pair(controller, tmp_path)
    controller.mutation_kind_buttons["deletion"].click()
    report = controller.model.report
    original_card = controller.evidence_card_label.text()
    assert "A:106" in original_card
    snapshot = controller._capture_panel_transaction_state()

    controller._invalidate_analysis_views()
    load_pair(controller)
    finish_report(controller)
    assert controller.model.report.report_id != report.report_id
    assert controller.result_labels["mutation_count"].text() == "0"

    controller._restore_panel_transaction_state(snapshot)
    assert controller.model.report == report
    assert controller.pages.currentWidget().objectName() == "pageResidues"
    assert controller.result_labels["mutation_count"].text() == "2"
    assert controller.residue_filter_combo.currentData() == "deletion"
    assert visible_residue_rows(controller) == [controller.residue_table.currentRow()]
    assert controller.evidence_card_label.text() == original_card


def test_clicking_a_specific_sequence_change_opens_that_position(controller, tmp_path):
    compare_changed_pair(controller, tmp_path, two_substitutions=True)
    controller.widget.resize(1200, 780)
    controller.widget.show()
    controller.nav.setCurrentRow(SCIENTIFIC_SECTIONS.index("Sequences"))
    QtWidgets.QApplication.processEvents()
    mutations = controller.mutation_table
    changed = [row for row in range(mutations.rowCount()) if mutations.item(row, 1).text() == "substitution"]
    assert len(changed) == 2
    item = mutations.item(changed[1], 0)
    mutations.scrollToItem(item)
    QtTest.QTest.mouseClick(
        mutations.viewport(), QtCore.Qt.MouseButton.LeftButton, pos=mutations.visualItemRect(item).center()
    )
    assert controller.pages.currentWidget().objectName() == "pageResidues"
    assert controller.residue_filter_combo.currentData() == "substitution"
    rows = visible_residue_rows(controller)
    assert len(rows) == 2
    assert controller.residue_table.currentRow() == rows[1]
    assert "A:105" in controller.evidence_card_label.text()
    assert "A:205" in controller.evidence_card_label.text()

    controller.nav.setCurrentRow(SCIENTIFIC_SECTIONS.index("Sequences"))
    conserved = mutations.item(0, 0)
    mutations.scrollToItem(conserved)
    QtWidgets.QApplication.processEvents()
    pos = mutations.visualItemRect(conserved).center()
    QtTest.QTest.mouseClick(mutations.viewport(), QtCore.Qt.MouseButton.LeftButton, pos=pos)
    QtTest.QTest.mouseDClick(mutations.viewport(), QtCore.Qt.MouseButton.LeftButton, pos=pos)
    assert controller.pages.currentWidget().objectName() == "pageResidues"
    assert controller.residue_filter_combo.currentData() == "all"
    assert controller.residue_table.currentRow() == 0
    assert "A:101" in controller.evidence_card_label.text()
    assert "A:201" in controller.evidence_card_label.text()
