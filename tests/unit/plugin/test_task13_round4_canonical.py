"""Round-four RED tests for canonical bindings, projects, and integrations."""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import pytest
from PySide6.QtWidgets import QApplication

from structlens.application.project_state import ProjectState
from structlens.core.models import AnalysisResult
from structlens.plugin.gui import qt_panel
from structlens.plugin.gui.main_panel import build_qt_panel

from ._task13_round4_fixtures import rich_report


@pytest.fixture(scope="module")
def application() -> QApplication:
    app = QApplication.instance() or QApplication([])
    yield app
    app.quit()


def _load_report_on_panel(controller: object, request: object, report: object) -> None:
    assert hasattr(request, "reference_selection")
    assert hasattr(request, "target_selection")
    reference_path = Path(request.reference_selection.path)  # type: ignore[attr-defined]
    target_path = Path(request.target_selection.path)  # type: ignore[attr-defined]
    assert controller._load_source("reference", reference_path)  # type: ignore[attr-defined]
    assert controller._load_source("target", target_path)  # type: ignore[attr-defined]
    controller._report_request = request  # type: ignore[attr-defined]
    controller._analysis_finished(report)  # type: ignore[attr-defined]


def test_save_project_persists_verified_report_and_content_addressed_artifacts(
    application, tmp_path: Path, monkeypatch
) -> None:
    """A normal report save must contain its exact report and both snapshots."""

    request, report = rich_report(tmp_path)
    project_path = tmp_path / "accepted.structlens.json"
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    _load_report_on_panel(controller, request, report)
    monkeypatch.setattr(
        controller.w.QFileDialog,
        "getSaveFileName",
        lambda *_args: (str(project_path), "StructLens project (*.json)"),
    )

    controller._save_project()

    assert project_path.is_file()
    payload = json.loads(project_path.read_text(encoding="utf-8"))
    v04 = payload["v04"]
    assert v04["report_hash"] == report.report_id
    assert v04["report_verification"] == "verified"
    assert v04["reference_selection"]["selection_id"] == report.reference_selection.selection_id
    assert v04["target_selection"]["selection_id"] == report.target_selection.selection_id
    assert v04["report_path"]
    report_path = project_path.parent / v04["report_path"]
    assert report_path.is_file()
    artifact_root = Path(f"{project_path}.artifacts") / report.report_id
    assert (artifact_root / "report.json").is_file()
    assert tuple((artifact_root / "snapshots").glob("*.snapshot")), "snapshot evidence must be persisted"
    assert payload.get("analysis_results") == []
    controller.close()
    panel.deleteLater()


def test_open_project_roundtrips_the_same_report_and_snapshot_bound_selection(
    application, tmp_path: Path, monkeypatch
) -> None:
    """Open must rehydrate the canonical artifact instead of a legacy result."""

    request, report = rich_report(tmp_path)
    project_path = tmp_path / "roundtrip.structlens.json"
    first = build_qt_panel(command=None)
    first_controller = first._structlens_controller
    _load_report_on_panel(first_controller, request, report)
    monkeypatch.setattr(
        first_controller.w.QFileDialog,
        "getSaveFileName",
        lambda *_args: (str(project_path), "StructLens project (*.json)"),
    )
    first_controller._save_project()
    first_controller.close()
    first.deleteLater()

    second = build_qt_panel(command=None)
    second_controller = second._structlens_controller
    errors: list[str] = []
    second_controller._show_error = errors.append
    monkeypatch.setattr(
        second_controller.w.QFileDialog,
        "getOpenFileName",
        lambda *_args: (str(project_path), "StructLens project (*.json)"),
    )

    second_controller._open_project()

    assert not errors
    assert second_controller.model.report is not None
    assert second_controller.model.report.report_id == report.report_id
    assert second_controller.model.analysis is None
    assert second_controller._report_request is not None
    assert second_controller._report_request.reference_snapshot.content_id == request.reference_snapshot.content_id
    assert second_controller._report_request.target_snapshot.raw_sha256 == request.target_snapshot.raw_sha256
    second_controller.close()
    second.deleteLater()


def test_open_project_replaces_a_stale_controller_request_before_rebinding(
    application, tmp_path: Path, monkeypatch
) -> None:
    """Opening a canonical project must not reuse a previous report request."""

    primary = tmp_path / "primary"
    secondary = tmp_path / "secondary"
    primary.mkdir()
    secondary.mkdir()
    request, report = rich_report(primary)
    project_path = tmp_path / "stale-request.structlens.json"
    saved = build_qt_panel(command=None)
    saved_controller = saved._structlens_controller
    _load_report_on_panel(saved_controller, request, report)
    monkeypatch.setattr(
        saved_controller.w.QFileDialog,
        "getSaveFileName",
        lambda *_args: (str(project_path), "StructLens project (*.json)"),
    )
    saved_controller._save_project()
    saved_controller.close()
    saved.deleteLater()

    stale_request, stale_report = rich_report(secondary)
    reopened = build_qt_panel(command=None)
    reopened_controller = reopened._structlens_controller
    _load_report_on_panel(reopened_controller, stale_request, stale_report)
    errors: list[str] = []
    reopened_controller._show_error = errors.append
    monkeypatch.setattr(
        reopened_controller.w.QFileDialog,
        "getOpenFileName",
        lambda *_args: (str(project_path), "StructLens project (*.json)"),
    )

    reopened_controller._open_project()

    assert not errors
    assert reopened_controller.model.binding is not None
    assert reopened_controller.model.binding.report.report_id == report.report_id
    assert reopened_controller.model.binding.request is reopened_controller._report_request
    reopened_controller.close()
    reopened.deleteLater()


def test_normal_report_state_contains_one_frozen_verified_binding(application, tmp_path: Path) -> None:
    """Rendering a report must retain report, request, snapshots, and presentation together."""

    request, report = rich_report(tmp_path)
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    _load_report_on_panel(controller, request, report)
    binding = getattr(controller.model, "binding", None)

    assert binding is not None, "Task 13 requires one verified canonical report binding"
    assert binding.report is report
    assert binding.request is request
    assert binding.snapshots == (request.reference_snapshot, request.target_snapshot)
    assert binding.presentation.report_id == report.report_id
    with pytest.raises((AttributeError, TypeError)):
        binding.report = None
    controller.close()
    panel.deleteLater()


def test_save_project_without_verified_binding_fails_closed(application, tmp_path: Path, monkeypatch) -> None:
    """A render-only compatibility report may not be persisted as canonical evidence."""

    _request, report = rich_report(tmp_path)
    output = tmp_path / "unverified.structlens.json"
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    # Deliberately use the render-only adapter and omit the request/binding.
    controller.model = controller.model.with_report(report)
    errors: list[str] = []
    controller._show_error = errors.append
    monkeypatch.setattr(
        controller.w.QFileDialog,
        "getSaveFileName",
        lambda *_args: (str(output), "StructLens project (*.json)"),
    )

    controller._save_project()

    assert errors
    assert not output.exists()
    controller.close()
    panel.deleteLater()


def test_tampered_canonical_report_rolls_back_open_without_replacing_previous_state(
    application, tmp_path: Path, monkeypatch
) -> None:
    """A failed canonical open is transactional and never reports success."""

    request, report = rich_report(tmp_path)
    project_path = tmp_path / "tampered.structlens.json"
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    _load_report_on_panel(controller, request, report)
    monkeypatch.setattr(
        controller.w.QFileDialog,
        "getSaveFileName",
        lambda *_args: (str(project_path), "StructLens project (*.json)"),
    )
    controller._save_project()
    before_report = controller.model.report
    before_status = controller.footer_status.text()
    payload = json.loads(project_path.read_text(encoding="utf-8"))
    assert isinstance(payload["v04"].get("report_path"), str) and payload["v04"]["report_path"]
    report_path = project_path.parent / payload["v04"]["report_path"]
    report_path.write_bytes(report_path.read_bytes() + b"tampered")
    errors: list[str] = []
    controller._show_error = errors.append
    monkeypatch.setattr(
        controller.w.QFileDialog,
        "getOpenFileName",
        lambda *_args: (str(project_path), "StructLens project (*.json)"),
    )

    controller._open_project()

    assert errors
    assert controller.model.report is before_report
    assert controller.footer_status.text() == before_status
    controller.close()
    panel.deleteLater()


def test_traversal_in_canonical_project_reference_is_rejected_without_state_mutation(
    application, tmp_path: Path, monkeypatch
) -> None:
    """Project metadata cannot escape its content-addressed artifact root."""

    request, report = rich_report(tmp_path)
    project_path = tmp_path / "traversal.structlens.json"
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    _load_report_on_panel(controller, request, report)
    monkeypatch.setattr(
        controller.w.QFileDialog,
        "getSaveFileName",
        lambda *_args: (str(project_path), "StructLens project (*.json)"),
    )
    controller._save_project()
    before_report = controller.model.report
    payload = json.loads(project_path.read_text(encoding="utf-8"))
    payload["v04"]["report_path"] = "../outside/report.json"
    project_path.write_text(json.dumps(payload), encoding="utf-8")
    errors: list[str] = []
    controller._show_error = errors.append
    monkeypatch.setattr(
        controller.w.QFileDialog,
        "getOpenFileName",
        lambda *_args: (str(project_path), "StructLens project (*.json)"),
    )

    controller._open_project()

    assert errors
    assert controller.model.report is before_report
    controller.close()
    panel.deleteLater()


def test_missing_canonical_snapshot_rolls_back_open_without_replacing_previous_state(
    application, tmp_path: Path, monkeypatch
) -> None:
    """Every referenced snapshot is required before a project can commit open."""

    request, report = rich_report(tmp_path)
    project_path = tmp_path / "missing-snapshot.structlens.json"
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    _load_report_on_panel(controller, request, report)
    monkeypatch.setattr(
        controller.w.QFileDialog,
        "getSaveFileName",
        lambda *_args: (str(project_path), "StructLens project (*.json)"),
    )
    controller._save_project()
    payload = json.loads(project_path.read_text(encoding="utf-8"))
    assert isinstance(payload["v04"].get("report_path"), str) and payload["v04"]["report_path"]
    artifact_root = Path(f"{project_path}.artifacts") / report.report_id
    snapshots = tuple((artifact_root / "snapshots").glob("*.snapshot"))
    assert snapshots
    snapshots[0].unlink()
    before_report = controller.model.report
    errors: list[str] = []
    controller._show_error = errors.append
    monkeypatch.setattr(
        controller.w.QFileDialog,
        "getOpenFileName",
        lambda *_args: (str(project_path), "StructLens project (*.json)"),
    )

    controller._open_project()

    assert errors
    assert controller.model.report is before_report
    controller.close()
    panel.deleteLater()


def test_invalid_canonical_visualization_state_rolls_back_open_without_replacing_previous_state(
    application, tmp_path: Path, monkeypatch
) -> None:
    """Visualization metadata must be validated before canonical GUI mutation."""

    request, report = rich_report(tmp_path)
    project_path = tmp_path / "invalid-visualization.structlens.json"
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    _load_report_on_panel(controller, request, report)
    monkeypatch.setattr(
        controller.w.QFileDialog,
        "getSaveFileName",
        lambda *_args: (str(project_path), "StructLens project (*.json)"),
    )
    controller._save_project()
    payload = json.loads(project_path.read_text(encoding="utf-8"))
    payload["visualization_state"]["color_mode"] = "not-a-real-color-mode"
    project_path.write_text(json.dumps(payload), encoding="utf-8")
    before_report = controller.model.report
    before_status = controller.footer_status.text()
    errors: list[str] = []
    controller._show_error = errors.append
    monkeypatch.setattr(
        controller.w.QFileDialog,
        "getOpenFileName",
        lambda *_args: (str(project_path), "StructLens project (*.json)"),
    )

    controller._open_project()

    assert errors
    assert controller.model.report is before_report
    assert controller.footer_status.text() == before_status
    controller.close()
    panel.deleteLater()


def test_legacy_project_is_explicitly_unverified_and_does_not_become_canonical_report(
    application, tmp_path: Path, monkeypatch
) -> None:
    """Legacy AnalysisResult projects remain isolated compatibility state."""

    legacy = AnalysisResult(
        reference_id="legacy-reference",
        target_id="legacy-target",
        correspondences=(),
        mutations=(),
        sequence_identity=1.0,
        sequence_coverage=1.0,
        alignment_decision="legacy",
    )
    project_path = tmp_path / "legacy.structlens.json"
    ProjectState(analysis_results=(legacy,)).save(project_path)
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    errors: list[str] = []
    controller._show_error = errors.append
    monkeypatch.setattr(
        controller.w.QFileDialog,
        "getOpenFileName",
        lambda *_args: (str(project_path), "StructLens project (*.json)"),
    )

    controller._open_project()

    assert not errors
    assert controller.model.report is None
    assert controller.model.analysis == legacy
    assert (
        "legacy" in controller.footer_status.text().lower() or "unverified" in controller.footer_status.text().lower()
    )
    controller.close()
    panel.deleteLater()


def test_canonical_export_rejects_a_report_request_mismatch_before_invoking_exporter(
    application, tmp_path: Path, monkeypatch
) -> None:
    """JSON/CSV/XLSX routes must not export a report with stale provenance."""

    request, report = rich_report(tmp_path)
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    _load_report_on_panel(controller, request, report)
    stale = replace(request, analysis_settings=replace(request.analysis_settings, minimum_sequence_identity=0.99))
    controller._report_request = stale
    calls: list[tuple[object, object, object]] = []
    errors: list[str] = []
    controller._show_error = errors.append
    monkeypatch.setattr(
        controller.w.QFileDialog,
        "getSaveFileName",
        lambda *_args: (str(tmp_path / "should-not-exist.json"), "JSON (*.json)"),
    )

    def exporter(value: object, path: object, **kwargs: object) -> None:
        calls.append((value, path, kwargs))

    controller._export("json", exporter, "JSON")

    assert not calls
    assert errors and ("stale" in errors[-1].lower() or "verified" in errors[-1].lower())
    controller.close()
    panel.deleteLater()


def test_report_visualization_does_not_convert_to_legacy_analysis(application, tmp_path: Path, monkeypatch) -> None:
    """Normal PyMOL Apply consumes the canonical report boundary directly."""

    request, report = rich_report(tmp_path)
    panel = build_qt_panel(command=object())
    controller = panel._structlens_controller
    _load_report_on_panel(controller, request, report)
    controller.reference_object_name = "reference-object"
    controller.target_object_name = "target-object"
    calls: list[tuple[object, dict[str, object]]] = []

    class SpyPyMol:
        def apply(self, value: object, **kwargs: object) -> None:
            calls.append((value, kwargs))

        def reset(self) -> None:
            return None

    controller._pymol = SpyPyMol()

    def forbidden(*_args: object, **_kwargs: object) -> object:
        raise AssertionError("normal report visualization must not call legacy_analysis")

    monkeypatch.setattr("structlens.plugin.gui.qt_visualization_mixin._legacy_analysis", forbidden)
    controller._apply_visualization()

    assert calls
    assert calls[0][0] is report
    assert calls[0][1]["reference_object"] == "reference-object"
    assert calls[0][1]["target_object"] == "target-object"
    controller.close()
    panel.deleteLater()


def test_canonical_pymol_bundle_route_fails_closed_until_snapshot_native_writer_exists(
    application, tmp_path: Path, monkeypatch
) -> None:
    """Task 13 cannot silently hand a mutable-path legacy bundle to users."""

    request, report = rich_report(tmp_path)
    output = tmp_path / "blocked.structlens-pymol"
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    _load_report_on_panel(controller, request, report)
    calls: list[object] = []
    errors: list[str] = []
    controller._show_error = errors.append
    monkeypatch.setattr(
        controller.w.QFileDialog,
        "getSaveFileName",
        lambda *_args: (str(output), "StructLens-PyMOL bundle (*.structlens-pymol)"),
    )

    def forbidden_writer(*args: object, **kwargs: object) -> object:
        calls.append((args, kwargs))
        return output

    monkeypatch.setattr(qt_panel, "write_pymol_bundle", forbidden_writer)
    controller._export_pymol_bundle()

    assert not calls
    assert not output.exists()
    assert errors and "unavailable for this report format" in errors[-1].lower()
    controller.close()
    panel.deleteLater()


def test_tsv_export_is_a_real_user_reachable_button_using_report_and_snapshots(
    application, tmp_path: Path, monkeypatch
) -> None:
    """TSV must be wired to a visible action, not only a private generic helper."""

    request, report = rich_report(tmp_path)
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    _load_report_on_panel(controller, request, report)
    buttons = [button for button in controller.widget.findChildren(controller.w.QPushButton) if "TSV" in button.text()]
    assert buttons, "Export page must expose a TSV action"
    exporter = getattr(qt_panel, "export_analysis_tsv", None)
    assert callable(exporter), "TSV action must route to the canonical TSV exporter"
    captured: list[tuple[object, object, object]] = []
    monkeypatch.setattr(
        controller.w.QFileDialog,
        "getSaveFileName",
        lambda *_args: (str(tmp_path / "result.tsv"), "TSV (*.tsv)"),
    )
    monkeypatch.setattr(
        qt_panel,
        "export_analysis_tsv",
        lambda value, path, **kwargs: captured.append((value, path, kwargs)),
    )

    buttons[0].click()

    assert captured
    assert captured[0][0] is report
    assert captured[0][2]["snapshots"] == (request.reference_snapshot, request.target_snapshot)
    controller.close()
    panel.deleteLater()
