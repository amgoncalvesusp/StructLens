"""The ligand contact diagram and residue inspector select the same position."""

import json
from pathlib import Path

import pytest

from structlens.plugin.gui.qt_compat import load_qt


def test_contact_panel_selection_and_no_coordinate_contacts(monkeypatch):
    monkeypatch.setenv("QT_QPA_PLATFORM", "offscreen")
    qt = load_qt()
    if qt is None:
        pytest.skip("Qt unavailable")
    from structlens.plugin.gui.main_panel import build_qt_panel
    from tests.unit.plugin.test_gui_workflow_usability import finish_report, load_pair

    app = qt.widgets.QApplication.instance() or qt.widgets.QApplication([])
    panel = build_qt_panel(command=None)
    c = panel._structlens_controller
    try:
        load_pair(c)
        finish_report(c)
        combo = c.site_ligand_combo
        atp = next(i for i in range(combo.count()) if "ATP" in str(combo.itemData(i)))
        combo.setCurrentIndex(atp)
        c._refresh_ligand_contacts()
        assert c.ligand_contact_table.rowCount() > 0
        assert "4.5" in c.ligand_contact_status.text()
        c.ligand_contact_table.setCurrentCell(0, 0)
        selected = c.ligand_contact_table.item(0, 0).text()
        assert "GLY" in selected or "MSE" in selected
        c.contact_evidence_button.click()
        assert c.pages.currentWidget().objectName() == "pageResidues"
        assert c.model.selection.reference_residue.residue_name in selected
        assert c.ligand_contact_scene.selectedItems()
        c._sync_ligand_contact_selection(None)
        assert not c.ligand_contact_scene.selectedItems()
        c._invalidate_analysis_views()
        assert c.ligand_contact_table.rowCount() == 0
    finally:
        c.close()
        panel.deleteLater()
        app.processEvents()


def test_external_pose_uses_reported_profile_and_links_click_to_residue(tmp_path, monkeypatch):
    from PySide6 import QtCore, QtTest, QtWidgets

    from structlens.plugin.gui.main_panel import SCIENTIFIC_SECTIONS, build_qt_panel
    from tests.unit.application.test_ligand_contacts import document
    from tests.unit.plugin.test_gui_workflow_usability import finish_report

    monkeypatch.setenv("QT_QPA_PLATFORM", "offscreen")
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    panel = build_qt_panel(command=None)
    c = panel._structlens_controller
    receptor = tmp_path / "one-model.pdb"
    fixture = Path("tests/fixtures/parsing/enriched.pdb").read_text()
    receptor.write_text(fixture.split("ENDMDL")[0] + "ENDMDL\nEND\n")
    try:
        for role in ("reference", "target"):
            assert c._load_source(role, receptor)
        c.mode_combo.setCurrentText("Sequence")
        finish_report(c)
        assert c.contact_import_button.isEnabled()
        raw = document(c.model.binding.request.reference_snapshot.raw_sha256)
        raw["interactions"][0]["receptor_residue"] = "ALA100A"
        raw["interactions"][0]["chemistry_basis"] = '<img src="file://example/hidden">'
        raw["interactions"].append({**raw["interactions"][0], "receptor_residue": "ASN444Z"})
        raw["interactions"].extend({**raw["interactions"][0], "ligand_atom": f"L{i}"} for i in range(2, 12))
        raw["poses"][0]["total_interactions"] = 12
        raw["totals"]["interactions"] = 12
        result = tmp_path / "result.json"
        result.write_text(json.dumps(raw))
        poses = tmp_path / "poses.pdbqt"
        poses.write_bytes(b"pose")
        paths = iter(((str(result), ""), (str(poses), "")))
        monkeypatch.setattr(QtWidgets.QFileDialog, "getOpenFileName", lambda *args: next(paths))
        c.contact_import_button.click()
        assert c.contact_source.currentData() == "p1"
        assert "LUNA-like × DSV-like (conservative)" in c.ligand_contact_status.text()
        assert c.ligand_contact_table.item(0, 1).text() == "hbond"
        assert c.ligand_contact_table.rowCount() == 8
        c.contact_next.click()
        assert c.ligand_contact_table.rowCount() == 4
        assert c.ligand_contact_table.item(0, 3).text() == "L8"
        c.contact_prev.click()
        c.nav.setCurrentRow(SCIENTIFIC_SECTIONS.index("Sites"))
        panel.resize(1100, 760)
        panel.show()
        app.processEvents()
        node = c._contact_nodes[0]
        point = c.ligand_contact_view.mapFromScene(node.sceneBoundingRect().center())
        QtTest.QTest.mouseClick(c.ligand_contact_view.viewport(), QtCore.Qt.MouseButton.LeftButton, pos=point)
        assert c.ligand_contact_table.currentRow() == 0
        assert c.residue_table.item(c.residue_table.currentRow(), 1).text() == "A:100 ALA"
        assert "A:100 ALA" in c.evidence_card_label.text()
        assert "<img" not in c.ligand_contact_table.item(0, 5).toolTip()
        snapshot = c._capture_panel_transaction_state()
        c._invalidate_analysis_views()
        assert c.contact_source.count() == 1
        c._restore_panel_transaction_state(snapshot)
        assert c.contact_source.currentData() == "p1"
        assert c.ligand_contact_table.currentRow() == 0
        assert c._docklens_contacts is not None
        c.ligand_contact_table.setCurrentCell(1, 0)
        assert c.model.selection is None
        assert c.residue_table.currentRow() == -1
        assert not c.contact_evidence_button.isEnabled()
        assert "does not map uniquely" in c.evidence_card_label.text()
        c.contact_source.setCurrentIndex(c.contact_source.findData("p2"))
        assert c.ligand_contact_table.rowCount() == 0
        assert "0 reported contacts" in c.ligand_contact_status.text()
    finally:
        c.close()
        panel.close()
        panel.deleteLater()
        app.processEvents()
