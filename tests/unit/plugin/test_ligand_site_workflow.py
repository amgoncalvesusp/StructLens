"""Protein-site membership and desktop ligand selection agree on exact deposited components."""

import os
from dataclasses import replace
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from PySide6 import QtWidgets  # noqa: E402

from structlens.application.report_service import ReportService  # noqa: E402
from structlens.core.models import AlignmentMode, AnalysisSettings  # noqa: E402
from structlens.core.parsing import load_structure_evidence  # noqa: E402
from structlens.core.sites import SiteDefinition, SiteDefinitionMode  # noqa: E402
from structlens.plugin.gui.main_panel import build_qt_panel  # noqa: E402
from structlens.plugin.gui.report_controller import AnalysisReportController  # noqa: E402

FIXTURE = Path(__file__).parents[2] / "fixtures" / "parsing" / "enriched.pdb"


@pytest.fixture
def panel():
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    widget = build_qt_panel(command=None)
    try:
        yield widget._structlens_controller
    finally:
        widget._structlens_controller.close()
        widget.deleteLater()
        app.processEvents()


def test_ligand_site_uses_same_members_in_metrics_and_evidence_cards(tmp_path):
    # A second, distant copy prevents attribution of one ligand's site to its namesake.
    source = FIXTURE.read_text()
    atp = next(line for line in source.splitlines() if "ATP A 202" in line)
    distant = atp[:6] + "   55" + atp[11:22] + " 305" + atp[26:30] + "  80.000" + atp[38:]
    path = tmp_path / "ligands.pdb"
    path.write_text(source.replace(atp, atp + "\n" + distant))
    request = AnalysisReportController().build_request(
        path, path, reference_model_id="1", reference_chain_id="A", target_model_id="1", target_chain_id="A",
        analysis_settings=AnalysisSettings(alignment_mode=AlignmentMode.SEQUENCE),
    )
    parsed = load_structure_evidence(request.reference_snapshot, selection=request.reference_selection)
    atps = [component for component in parsed.components if component.residue_name == "ATP"]
    assert len(atps) == 2
    sites = tuple(
        SiteDefinition(f"atp-{i}", mode=SiteDefinitionMode.LIGAND_RADIUS, ligand_id=c.component_id, radius_angstrom=4.0)
        for i, c in enumerate(atps)
    )
    report = ReportService().analyze(replace(request, site_definitions=sites))
    assert len(report.sites) == 1
    metric = report.sites[0]
    assert metric.mapped_residue_count == 2 and metric.coverage_fraction == 1.0
    members = {c.residue_id.auth_seq_id + (c.residue_id.insertion_code or "") for c in report.evidence_cards if c.site.metrics}
    assert members == {"100A", "101"}
    assert all(card.site.metrics in ((), (metric,)) for card in report.evidence_cards)


def test_site_picker_preview_names_and_reference_changes(panel):
    assert panel._load_source("reference", FIXTURE)
    panel.site_mode_combo.setCurrentText("Ligand radius")
    index = next(i for i in range(panel.site_ligand_combo.count()) if "ATP" in panel.site_ligand_combo.itemText(i))
    panel.site_ligand_combo.setCurrentIndex(index)
    ligand_id = panel.site_ligand_combo.currentData()
    panel.site_radius_spin.setValue(4.0)
    assert "A:100A" in panel.site_preview_label.text() and "A:101" in panel.site_preview_label.text()
    assert "A:100 ALA" not in panel.site_preview_label.text()
    panel.site_name_edit.setText("ATP pocket")
    panel.site_define_button.click()
    panel.site_mode_combo.setCurrentText("Key residues")
    panel.site_residues_edit.setText("A:100")
    panel.site_name_edit.setText("Catalytic")
    panel.site_define_button.click()
    assert [site.name for site in panel._site_definitions] == ["ATP pocket", "Catalytic"]
    assert panel._site_definitions[0].ligand_id == ligand_id
    panel.site_define_button.click()
    assert len(panel._site_definitions) == 2
    assert panel._load_source("target", FIXTURE)
    assert len(panel._site_definitions) == 2
    assert panel._load_source("reference", FIXTURE)
    assert len(panel._site_definitions) == 2
    panel.configured_site_combo.setCurrentIndex(1)
    panel.site_remove_button.click()
    assert [site.name for site in panel._site_definitions] == ["ATP pocket"]
    assert panel._load_source("reference", FIXTURE.with_name("numbering_altloc.pdb"))
    assert not panel._site_definitions
    assert panel.site_ligand_combo.findData(ligand_id) == -1


def test_open_project_restores_named_ligand_sites(panel, tmp_path, monkeypatch):
    for role in ("reference", "target"):
        assert panel._load_source(role, FIXTURE)
    panel.site_mode_combo.setCurrentText("Ligand radius")
    panel.site_name_edit.setText("ATP shell")
    panel.site_radius_spin.setValue(4)
    panel.site_define_button.click()
    definitions = panel._site_definitions
    request = panel._report_controller.build_request(
        FIXTURE, FIXTURE, reference_model_id="1", reference_chain_id="A", target_model_id="1", target_chain_id="A",
        site_definitions=definitions, analysis_settings=panel._settings(),
    )
    panel._populate_report(ReportService().analyze(request), request=request)
    path = tmp_path / "project.json"
    monkeypatch.setattr(QtWidgets.QFileDialog, "getSaveFileName", lambda *_: (str(path), ""))
    panel._save_project()
    assert path.exists()
    assert panel._load_source("reference", FIXTURE.with_name("numbering_altloc.pdb"))
    monkeypatch.setattr(QtWidgets.QFileDialog, "getOpenFileName", lambda *_: (str(path), ""))
    panel._open_project()
    assert panel.configured_site_combo.currentText() == "ATP shell"
    assert panel.site_ligand_combo.currentData() == definitions[0].ligand_id
    assert panel._site_definitions == definitions


def test_ref_model_switch_clears_inapplicable_ligand_and_preview(panel):
    assert panel._load_source("reference", FIXTURE)
    panel.site_mode_combo.setCurrentText("Ligand radius")
    assert panel.site_ligand_combo.currentData() is not None
    panel.site_define_button.click()
    assert panel._site_definitions
    panel.reference_model_combo.setCurrentIndex(panel.reference_model_combo.findData("2"))
    assert panel.site_ligand_combo.currentData() is None
    assert not panel._site_definitions
    assert "ligand identifier" in panel.site_preview_label.text()
