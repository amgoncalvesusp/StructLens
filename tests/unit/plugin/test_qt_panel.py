import os
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
qt_widgets = pytest.importorskip("PySide6.QtWidgets")

from structlens.application.chart_data import ChartDataset, ChartSeries  # noqa: E402
from structlens.application.msa_service import align_sequences  # noqa: E402
from structlens.core.models import AnalysisResult, ResidueId  # noqa: E402
from structlens.core.msa import AnalysisSequence, MSASettings, SequenceResidueRef  # noqa: E402
from structlens.core.sites import SiteMetrics  # noqa: E402
from structlens.plugin.gui.main_panel import SCIENTIFIC_SECTIONS, build_qt_panel  # noqa: E402


@pytest.fixture(scope="module")
def application():
    app = qt_widgets.QApplication.instance() or qt_widgets.QApplication([])
    yield app
    app.quit()


def test_qt_panel_builds_operate_workflow(application) -> None:
    panel = build_qt_panel()
    controller = panel._structlens_controller

    assert panel.objectName() == "structlensPanel"
    assert controller.nav.count() == 9
    assert controller.pages.count() == 9
    assert controller.compare_button.text() == "Compare"
    assert controller.mode_combo.count() == 4
    assert controller.comparison_combo.count() == 1
    assert controller.comparison_combo.currentData() == "pairwise"
    assert controller.mutation_table.columnCount() == 8
    assert controller.residue_table.columnCount() == 10
    assert controller.msa_table.columnCount() == 3
    assert controller.nav.item(4).text() == "Sites"
    assert controller.nav.item(6).text() == "PyMOL"
    assert controller.nav.item(8).text() == "Export"

    controller.close()
    panel.deleteLater()


def test_chart_exports_follow_the_selected_profile(application) -> None:
    panel = build_qt_panel()
    controller = panel._structlens_controller

    assert controller.chart_export_buttons
    assert all(button.isEnabled() for button in controller.chart_export_buttons)

    controller.chart_combo.setCurrentText("Mutation / conservation matrix")
    assert all(not button.isEnabled() for button in controller.chart_export_buttons)

    controller.chart_combo.setCurrentText("Structural deviation profile")
    assert all(button.isEnabled() for button in controller.chart_export_buttons)

    controller.close()
    panel.deleteLater()


def test_analysis_and_manual_recovery_land_on_the_relevant_pages(application) -> None:
    panel = build_qt_panel()
    controller = panel._structlens_controller
    controller._show_error = lambda _: None

    controller._analysis_finished(
        AnalysisResult(
            reference_id="reference",
            target_id="target",
            correspondences=(),
            mutations=(),
            sequence_identity=0.0,
            sequence_coverage=0.0,
            alignment_decision="test",
        )
    )
    assert controller.nav.currentRow() == SCIENTIFIC_SECTIONS.index("Structures")
    assert controller.structure_result_table.rowCount() == 1
    assert controller.results_table.rowCount() == 1

    fixture = Path("tests/fixtures/parsing/numbering_altloc.pdb")
    assert controller._load_source("reference", fixture) is True
    assert controller._load_source("target", fixture) is True
    controller.mode_combo.setCurrentIndex(3)
    controller.manual_edit.setPlainText("invalid pair")
    controller._start_analysis()
    assert controller.nav.currentRow() == SCIENTIFIC_SECTIONS.index("Structures")

    controller.close()
    panel.deleteLater()


def test_standalone_mode_keeps_file_workflow_without_pymol_actions(application) -> None:
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller

    assert controller.command is None
    assert any(
        "Standalone mode" in label.text()
        for label in panel.findChildren(qt_widgets.QLabel)
    )
    assert not any(
        button.text() == "Apply to PyMOL"
        for button in panel.findChildren(qt_widgets.QPushButton)
    )

    controller.close()
    panel.deleteLater()


def test_failed_source_reload_clears_previous_structure(application) -> None:
    panel = build_qt_panel()
    controller = panel._structlens_controller
    fixture = Path("tests/fixtures/parsing/numbering_altloc.pdb")
    assert controller._load_source("reference", fixture) is True
    assert controller.reference_structure is not None
    controller._show_error = lambda _: None

    assert controller._load_source("reference", Path("missing-input.pdb")) is False
    assert controller.reference_structure is None
    assert controller.reference_chain_combo.count() == 0

    controller.close()
    panel.deleteLater()


def test_results_page_compiles_all_analysis_results_and_preserves_history(application) -> None:
    panel = build_qt_panel()
    controller = panel._structlens_controller
    first = AnalysisResult(
        reference_id="reference",
        target_id="target-a",
        correspondences=(),
        mutations=(),
        sequence_identity=0.90,
        sequence_coverage=0.80,
        alignment_decision="accepted",
    )
    second = AnalysisResult(
        reference_id="reference",
        target_id="target-b",
        correspondences=(),
        mutations=(),
        sequence_identity=0.75,
        sequence_coverage=0.70,
        alignment_decision="accepted",
    )
    controller._analysis_finished(first)
    controller._analysis_finished(second)
    assert controller.results_table.rowCount() == 2
    assert controller.results_table.item(0, 1).text() == "target-a"
    assert controller.results_table.item(1, 1).text() == "target-b"
    assert len(controller._analysis_history) == 2
    controller.close()
    panel.deleteLater()


def test_sites_and_charts_show_authoritative_result_areas(application) -> None:
    panel = build_qt_panel()
    controller = panel._structlens_controller
    metrics = SiteMetrics(
        site_id="active-site",
        structure_id="target",
        mapped_residue_count=4,
        coverage_fraction=0.8,
        global_frame_backbone_rmsd_angstrom=1.2,
        site_fitted_backbone_rmsd_angstrom=0.6,
        sasa_angstrom2=123.4,
        atomic_envelope_volume_angstrom3=456.7,
    )
    controller.set_site_metrics((metrics,))
    assert controller.site_metrics_table.rowCount() == 1
    assert controller.site_metrics_table.item(0, 0).text() == "active-site"
    dataset = ChartDataset("msa", "MSA", "column", "conservation", "fraction", (), "descriptive")
    controller.set_chart_datasets({"MSA conservation profile": dataset})
    controller.chart_combo.setCurrentText("MSA conservation profile")
    assert controller.chart_preview_status.text() != "Chart unavailable. Run the corresponding scientific service first."
    controller.close()
    panel.deleteLater()


def test_msa_viewer_renders_authoritative_alignment_rows(application) -> None:
    panel = build_qt_panel()
    controller = panel._structlens_controller
    def sequence(identifier: str, value: str) -> AnalysisSequence:
        return AnalysisSequence(
            identifier,
            "A",
            value,
            tuple(SequenceResidueRef(index, character, ResidueId(identifier, "1", "A", str(index), None, "ALA")) for index, character in enumerate(value)),
            "derived",
        )
    alignment = align_sequences((sequence("ref", "ABC"), sequence("target", "ABCGH")), MSASettings())
    controller.set_msa_result(alignment)
    assert controller.msa_table.rowCount() == 2
    assert controller.msa_table.item(0, 1).text() == "ABC--"
    assert "insertion columns" in controller.msa_summary_label.text()
    controller.close()
    panel.deleteLater()


def test_v03_xlsx_export_routes_staged_records(application, monkeypatch, tmp_path: Path) -> None:
    panel = build_qt_panel()
    controller = panel._structlens_controller
    controller.model = controller.model.with_analysis(
        AnalysisResult(
            reference_id="reference",
            target_id="target",
            correspondences=(),
            mutations=(),
            sequence_identity=1.0,
            sequence_coverage=1.0,
            alignment_decision="test",
        )
    )
    output = tmp_path / "result.xlsx"
    controller.set_v03_export_records(provenance=("fixture",))
    called: dict[str, object] = {}
    monkeypatch.setattr(
        controller.w.QFileDialog,
        "getSaveFileName",
        lambda *_args: (str(output), "XLSX (*.xlsx)"),
    )
    monkeypatch.setattr(
        "structlens.plugin.gui.qt_panel.export_v03_xlsx",
        lambda path, **kwargs: called.update(path=str(path), **kwargs),
    )
    controller._export_xlsx()
    assert called == {"path": str(output), "provenance": ("fixture",)}
    controller.close()
    panel.deleteLater()


def test_chart_exports_route_staged_v03_dataset(application, monkeypatch, tmp_path: Path) -> None:
    panel = build_qt_panel()
    controller = panel._structlens_controller
    controller.model = controller.model.with_analysis(
        AnalysisResult(
            reference_id="reference",
            target_id="target",
            correspondences=(),
            mutations=(),
            sequence_identity=1.0,
            sequence_coverage=1.0,
            alignment_decision="test",
        )
    )
    dataset = ChartDataset(
        "msa_conservation",
        "MSA conservation",
        "Alignment column",
        "Conservation",
        "fraction",
        (ChartSeries("reference", ((1.0, 1.0),)),),
        "Descriptive alignment conservation.",
    )
    controller.set_chart_datasets({"MSA conservation profile": dataset})
    controller.chart_combo.setCurrentText("MSA conservation profile")
    assert all(button.isEnabled() for button in controller.chart_export_buttons)
    output = tmp_path / "chart.xlsx"
    monkeypatch.setattr(controller.w.QFileDialog, "getSaveFileName", lambda *_args: (str(output), "XLSX (*.xlsx)"))
    captured: dict[str, object] = {}
    monkeypatch.setattr("structlens.plugin.gui.qt_panel.export_chart_xlsx", lambda value, path: captured.update(value=value, path=path))
    controller._export_chart_xlsx()
    assert captured["value"] == dataset
    controller.close()
    panel.deleteLater()


def test_new_analysis_invalidates_staged_v03_views(application) -> None:
    panel = build_qt_panel()
    controller = panel._structlens_controller
    controller.set_v03_export_records(provenance=("old",))
    controller.set_chart_datasets({"MSA conservation profile": ChartDataset("msa", "MSA", "x", "y", None, (), "")})
    controller.set_msa_result(align_sequences((
        AnalysisSequence("ref", "A", "A", (SequenceResidueRef(0, "A", ResidueId("ref", "1", "A", "0", None, "ALA")),), "derived"),
    ), MSASettings()))
    assert controller.msa_table.rowCount() == 1
    controller._analysis_finished(
        AnalysisResult(
            reference_id="new-reference",
            target_id="new-target",
            correspondences=(),
            mutations=(),
            sequence_identity=1.0,
            sequence_coverage=1.0,
            alignment_decision="test",
        )
    )
    assert controller._v03_export_records == {}
    assert controller._chart_datasets == {}
    assert controller.msa_table.rowCount() == 0
    controller.close()
    panel.deleteLater()


def test_pymol_bundle_export_forwards_only_staged_v03_payloads(
    application, monkeypatch, tmp_path: Path
) -> None:
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    fixture = Path("tests/fixtures/parsing/numbering_altloc.pdb")
    assert controller._load_source("reference", fixture) is True
    assert controller._load_source("target", fixture) is True

    controller.model = controller.model.with_analysis(
        AnalysisResult(
            reference_id="reference",
            target_id="target",
            correspondences=(),
            mutations=(),
            sequence_identity=1.0,
            sequence_coverage=1.0,
            alignment_decision="test",
        )
    )
    assert all(value is None for value in controller._v03_bundle_kwargs().values())

    payloads = {
        "msa_summary": {"columns": [{"index": 0}]},
        "conservation": {"columns": [{"reference_label": "A:100"}]},
        "interactions": {"differences": []},
        "sites": {"metrics": []},
        "evidence": {"cards": []},
        "vectors": {"vectors": []},
    }
    controller.set_v03_bundle_payloads(**payloads)
    forwarded: dict[str, object] = {}

    def fake_write_bundle(output_path, **kwargs):
        forwarded.update(kwargs)
        return Path(output_path)

    monkeypatch.setattr(
        "structlens.plugin.gui.qt_panel.write_pymol_bundle", fake_write_bundle
    )
    output = tmp_path / "analysis.structlens-pymol"
    monkeypatch.setattr(
        controller.w.QFileDialog,
        "getSaveFileName",
        lambda *_args: (str(output), "StructLens-PyMOL bundle (*.structlens-pymol)"),
    )

    controller._export_pymol_bundle()

    for name, payload in payloads.items():
        assert forwarded[name] == payload

    controller.close()
    panel.deleteLater()
