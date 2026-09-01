"""High-value coverage for the final Qt export/transaction boundaries."""

from __future__ import annotations

from pathlib import Path

import pytest
from PySide6.QtWidgets import QApplication

from structlens.application.chart_data import ChartDataset, ChartSeries
from structlens.application.project_state import ProjectState
from structlens.core.models import AnalysisResult
from structlens.plugin.gui.main_panel import build_qt_panel
from structlens.plugin.gui.model import CanonicalReportBinding, StructLensPanelModel
from structlens.plugin.gui.project_transaction import save_canonical_project, stage_project
from structlens.plugin.visualization.renderer import VisualizationState

from ._task13_round4_fixtures import rich_report
from ._task13_round5_fixtures import typed_report


@pytest.fixture(scope="module")
def application() -> QApplication:
    app = QApplication.instance() or QApplication([])
    yield app
    app.quit()


@pytest.fixture()
def controller(application):
    panel = build_qt_panel(command=None)
    value = panel._structlens_controller
    yield value
    value.close()
    panel.deleteLater()
    application.processEvents()


def _install_report(controller: object, request: object, report: object) -> None:
    reference_path = Path(request.reference_selection.path)  # type: ignore[attr-defined]
    target_path = Path(request.target_selection.path)  # type: ignore[attr-defined]
    assert controller._load_source("reference", reference_path)  # type: ignore[attr-defined]
    assert controller._load_source("target", target_path)  # type: ignore[attr-defined]
    controller._report_request = request  # type: ignore[attr-defined]
    controller._analysis_finished(report)  # type: ignore[attr-defined]


def _canonical_candidate(tmp_path: Path):
    saved = tmp_path / "saved"
    saved.mkdir()
    request, report = rich_report(saved)
    binding = CanonicalReportBinding.create(report, request)
    project_path = tmp_path / "saved.structlens.json"
    save_canonical_project(binding, project_path, visualization_state={})
    return stage_project(project_path), project_path


def _legacy_result() -> AnalysisResult:
    return AnalysisResult(
        reference_id="legacy-reference",
        target_id="legacy-target",
        correspondences=(),
        mutations=(),
        sequence_identity=0.75,
        sequence_coverage=0.8,
        alignment_decision="legacy",
    )


def test_canonical_commit_exception_restores_populated_report_and_sources(
    controller, monkeypatch, tmp_path: Path
) -> None:
    candidate, project_path = _canonical_candidate(tmp_path)
    current = tmp_path / "current"
    current.mkdir()
    current_request, current_report = typed_report(current)
    _install_report(controller, current_request, current_report)
    before_model = controller.model
    before_source = controller.reference_loaded_source
    before_footer = controller.footer_status.text()
    before_header = controller.header_status.text()

    def fail_commit(*_args: object, **_kwargs: object) -> None:
        raise RuntimeError("commit phase failed")

    monkeypatch.setattr(controller, "_commit_report", fail_commit)
    with pytest.raises(RuntimeError, match="commit phase failed"):
        controller._commit_canonical_project(candidate, project_path)

    assert controller.model is before_model
    assert controller.model.report is current_report
    assert controller._report_request is current_request
    assert controller.reference_loaded_source is not None
    assert before_source is not None
    assert controller.reference_loaded_source.path == before_source.path
    assert controller.reference_loaded_source.snapshot == before_source.snapshot
    assert controller.footer_status.text() == before_footer
    assert controller.header_status.text() == before_header
    assert controller.results_table.rowCount() == len(controller._presentation_history)


def test_canonical_commit_exception_restores_legacy_analysis_state(controller, monkeypatch, tmp_path: Path) -> None:
    candidate, project_path = _canonical_candidate(tmp_path)
    fixture = Path(__file__).parents[2] / "fixtures" / "parsing" / "numbering_altloc.pdb"
    assert controller._load_source("reference", fixture)
    assert controller._load_source("target", fixture)
    legacy = _legacy_result()
    controller.model = controller.model.with_analysis(legacy)
    controller._report_request = None
    before_model = controller.model

    monkeypatch.setattr(
        controller,
        "_commit_report",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("legacy rollback")),
    )
    with pytest.raises(RuntimeError, match="legacy rollback"):
        controller._commit_canonical_project(candidate, project_path)

    assert controller.model is before_model
    assert controller.model.analysis is legacy
    assert controller.model.report is None
    assert controller.result_decision.text().startswith("<b>legacy</b>")
    assert controller.reference_loaded_source is not None
    assert controller.reference_loaded_source.path == fixture


def test_canonical_commit_exception_restores_empty_panel_state(controller, monkeypatch, tmp_path: Path) -> None:
    candidate, project_path = _canonical_candidate(tmp_path)
    before_model = controller.model
    assert before_model == StructLensPanelModel()

    monkeypatch.setattr(
        controller,
        "_commit_report",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("empty rollback")),
    )
    with pytest.raises(RuntimeError, match="empty rollback"):
        controller._commit_canonical_project(candidate, project_path)

    assert controller.model is before_model
    assert controller.model.report is None
    assert controller.model.analysis is None
    assert controller.reference_loaded_source is None
    assert controller.target_loaded_source is None
    assert controller.footer_status.text() == before_model.status


def test_apply_visualization_state_restores_all_controls_and_immutable_model(controller) -> None:
    state = VisualizationState(
        preset="Publication",
        local_radius_angstrom=8.5,
        show_labels=True,
        show_reference=False,
        show_target=True,
    )
    before = controller.model
    controller._apply_visualization_state(state)

    assert controller.model is not before
    assert controller.model.visualization_state == state
    assert controller.preset_combo.currentText() == "Publication"
    assert controller.radius_spin.value() == pytest.approx(8.5)
    assert controller.labels_check.isChecked() is True
    assert controller.reference_check.isChecked() is False
    assert controller.target_check.isChecked() is True


def test_save_project_legacy_uses_explicit_source_object_and_history_state(
    controller, monkeypatch, tmp_path: Path
) -> None:
    legacy = _legacy_result()
    controller.model = controller.model.with_analysis(legacy)
    controller.reference_object_name = "reference-object"
    controller.target_object_name = "target-object"
    controller._analysis_history = (legacy,)
    output = tmp_path / "legacy.json"
    monkeypatch.setattr(controller.w.QFileDialog, "getSaveFileName", lambda *_: (str(output), "JSON"))

    controller._save_project()

    saved = ProjectState.load(output)
    assert saved.analysis_results == (legacy,)
    assert saved.reference_source is None
    assert saved.target_sources == ()
    assert saved.source_objects == {"reference": "reference-object", "target": "target-object"}
    assert "Project saved" in controller.footer_status.text()


def test_save_project_cancel_does_not_construct_or_write_legacy_payload(
    controller, monkeypatch, tmp_path: Path
) -> None:
    controller.model = controller.model.with_analysis(_legacy_result())
    monkeypatch.setattr(controller.w.QFileDialog, "getSaveFileName", lambda *_: ("", ""))
    output = tmp_path / "cancelled.json"

    controller._save_project()

    assert not output.exists()


def test_legacy_export_wrappers_route_csv_and_empty_xlsx_records(controller, monkeypatch, tmp_path: Path) -> None:
    legacy = _legacy_result()
    controller.model = controller.model.with_analysis(legacy)
    calls: list[tuple[object, object]] = []
    monkeypatch.setattr(controller.w.QFileDialog, "getSaveFileName", lambda *_: (str(tmp_path / "result.out"), "out"))
    monkeypatch.setattr(
        "structlens.plugin.gui.qt_exports_mixin.export_analysis_csv",
        lambda value, path, **_kwargs: calls.append((value, path)),
    )
    monkeypatch.setattr(
        "structlens.plugin.gui.qt_exports_mixin.export_analysis_xlsx",
        lambda value, path, **_kwargs: calls.append((value, path)),
    )

    controller._export_csv()
    controller._v03_export_records = {}
    controller._export_xlsx()

    assert calls == [(legacy, str(tmp_path / "result.out")), (legacy, str(tmp_path / "result.out"))]


def test_canonical_open_and_chart_exports_fail_before_external_side_effects(
    controller, monkeypatch, tmp_path: Path
) -> None:
    request, report = rich_report(tmp_path)
    _install_report(controller, request, report)
    errors: list[str] = []
    controller._show_error = errors.append
    writers: list[object] = []
    monkeypatch.setattr(
        "structlens.plugin.gui.qt_exports_mixin.write_pymol_bundle",
        lambda *_args, **_kwargs: writers.append(True),
    )
    controller._open_in_pymol()
    assert writers == []
    assert "task 15" in errors[-1].lower()

    controller.model = StructLensPanelModel()
    controller._export_chart_xlsx()
    assert "comparison" in errors[-1].lower()


def test_canonical_binding_rejects_none_and_stale_presentation(controller, monkeypatch, tmp_path: Path) -> None:
    request, report = rich_report(tmp_path)
    _install_report(controller, request, report)
    controller.model = StructLensPanelModel()
    with pytest.raises(ValueError, match="no canonical report"):
        controller._canonical_binding()

    controller.model = controller.model.with_report(report)
    controller._report_request = request
    controller.model = controller.model.with_binding(CanonicalReportBinding.create(report, request))
    controller._report_presentation = None
    with pytest.raises(ValueError, match="presentation is stale"):
        controller._canonical_binding()


def test_chart_export_cancel_preserves_dataset_and_status(controller, monkeypatch, tmp_path: Path) -> None:
    legacy = _legacy_result()
    controller.model = controller.model.with_analysis(legacy)
    dataset = ChartDataset(
        "legacy-chart",
        "Legacy chart",
        "position",
        "score",
        "fraction",
        (ChartSeries("target", ((1.0, 0.5),)),),
        "descriptive",
    )
    controller._chart_datasets = {"MSA conservation profile": dataset}
    controller.chart_combo.setCurrentText("MSA conservation profile")
    monkeypatch.setattr(controller.w.QFileDialog, "getSaveFileName", lambda *_: ("", ""))
    controller._export_chart_xlsx()
    assert controller._chart_datasets["MSA conservation profile"] is dataset
    assert "Legacy chart" not in controller.footer_status.text()


def test_site_definition_ligand_and_residue_radius_require_their_inputs(controller, monkeypatch) -> None:
    fixture = Path(__file__).parents[2] / "fixtures" / "parsing" / "numbering_altloc.pdb"
    assert controller._load_source("reference", fixture)
    errors: list[str] = []
    monkeypatch.setattr(controller, "_show_error", errors.append)
    controller.site_mode_combo.setCurrentText("Ligand radius")
    controller.site_ligand_edit.clear()
    controller._define_site_from_controls()
    assert "ligand identifier" in errors[-1]
    controller.site_mode_combo.setCurrentText("Residue radius")
    controller.site_residues_edit.clear()
    controller._define_site_from_controls()
    assert "center position" in errors[-1]
