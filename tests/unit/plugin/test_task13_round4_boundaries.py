"""Round-four RED tests for immutable source and report boundaries."""

from __future__ import annotations

from pathlib import Path

import pytest
from PySide6.QtWidgets import QApplication

from structlens.application.dto import AnalysisReportRequest
from structlens.core.models import AlignmentMode, AnalysisSettings, ResidueId
from structlens.core.parsing import capture_snapshot
from structlens.core.sites import SiteDefinition, SiteDefinitionMode
from structlens.plugin.gui.main_panel import build_qt_panel
from structlens.plugin.gui.report_controller import AnalysisReportController

from ._task13_round4_fixtures import (
    non_first_request_paths,
    rich_report,
    selection,
)


@pytest.fixture(scope="module")
def application() -> QApplication:
    app = QApplication.instance() or QApplication([])
    yield app
    app.quit()


def test_snapshot_request_builder_preserves_exact_non_first_selection_and_evidence(tmp_path: Path) -> None:
    """A request builder must not recapture bytes or replace model/chain IDs."""

    reference_path, target_path = non_first_request_paths(tmp_path)
    reference_snapshot = capture_snapshot(reference_path)
    target_snapshot = capture_snapshot(target_path)
    reference_selection = selection(reference_snapshot, model="2", chain="B", path=reference_path)
    target_selection = selection(target_snapshot, model="2", chain="B", path=target_path)
    reference_residue = ResidueId("multi-reference", "2", "B", "1", None, "ALA")
    target_residue = ResidueId("multi-target", "2", "B", "1", None, "ALA")
    site = SiteDefinition(
        "locked-round4-site",
        mode=SiteDefinitionMode.KEY_RESIDUES,
        reference_residues=(reference_residue,),
    )
    settings = AnalysisSettings(alignment_mode=AlignmentMode.MANUAL)
    pairs = ((reference_residue, target_residue),)
    controller = AnalysisReportController(service=object())

    builder = getattr(controller, "build_request_from_snapshots", None)
    assert callable(builder), "Task 13 requires a snapshot-native request builder"
    built = builder(
        reference_snapshot,
        target_snapshot,
        reference_selection=reference_selection,
        target_selection=target_selection,
        analysis_settings=settings,
        site_definitions=(site,),
        manual_pairs=pairs,
    )

    assert isinstance(built, AnalysisReportRequest)
    assert built.reference_snapshot is reference_snapshot
    assert built.target_snapshot is target_snapshot
    assert built.reference_selection == reference_selection
    assert built.target_selection == target_selection
    assert built.reference_selection.model_id == "2"
    assert built.reference_selection.author_chain_ids == ("B",)
    assert built.target_selection.model_id == "2"
    assert built.target_selection.author_chain_ids == ("B",)
    assert built.site_definitions == (site,)
    assert built.manual_pairs == pairs

    # Replacing the mutable source after capture cannot change the already
    # constructed request or cause a second parse to silently select model 1.
    reference_path.write_bytes(b"MODEL        1\nENDMDL\nEND\n")
    target_path.unlink()
    assert built.reference_snapshot is reference_snapshot
    assert built.target_snapshot is target_snapshot
    assert built.reference_selection.model_id == "2"
    assert built.target_selection.author_chain_ids == ("B",)


def test_compare_never_resets_non_first_model_chain_or_reparses_source(
    application, tmp_path: Path, monkeypatch
) -> None:
    """Pressing Compare must use the loaded snapshot, not a fresh first-item parse."""

    reference_path, target_path = non_first_request_paths(tmp_path)
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    errors: list[str] = []
    controller._show_error = errors.append
    assert controller._load_source("reference", reference_path)
    assert controller._load_source("target", target_path)
    controller.reference_model_combo.setCurrentText("Model 2")
    controller.target_model_combo.setCurrentText("Model 2")
    controller._model_changed("reference")
    controller._model_changed("target")
    controller.reference_chain_combo.setCurrentText("Chain B · 2 residues")
    controller.target_chain_combo.setCurrentText("Chain B · 2 residues")
    assert controller.reference_model_combo.currentData() == "2"
    assert controller.target_model_combo.currentData() == "2"
    assert controller.reference_chain_combo.currentData() == "B"
    assert controller.target_chain_combo.currentData() == "B"

    report_root = tmp_path / "report-data"
    report_root.mkdir()
    _, report = rich_report(report_root)
    captured: list[AnalysisReportRequest] = []

    class Service:
        def analyze(self, built: AnalysisReportRequest) -> object:
            captured.append(built)
            return report

    controller._report_controller._service = Service()
    monkeypatch.setattr(
        "structlens.plugin.gui.qt_sources_mixin._ANALYSIS_EXECUTOR",
        _ImmediateExecutor(),
    )

    def forbidden_reload() -> bool:
        raise AssertionError("Compare must not recapture and repopulate source controls")

    monkeypatch.setattr(controller, "_load_sources_from_edits", forbidden_reload)
    controller._start_analysis()

    assert not errors
    assert captured, "Compare did not submit the canonical request"
    built = captured[-1]
    assert built.reference_selection.model_id == "2"
    assert built.reference_selection.author_chain_ids == ("B",)
    assert built.target_selection.model_id == "2"
    assert built.target_selection.author_chain_ids == ("B",)
    controller.close()
    panel.deleteLater()


def test_missing_source_after_compare_preserves_previous_report_and_selection(application, tmp_path: Path) -> None:
    """Input loss is a failed next run, not permission to clear valid evidence."""

    current_request, report = rich_report(tmp_path)
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    assert controller._load_source("reference", Path(current_request.reference_selection.path))
    assert controller._load_source("target", Path(current_request.target_selection.path))
    controller._report_request = current_request
    controller._analysis_finished(report)
    before = controller.model
    Path(current_request.reference_selection.path).unlink()  # type: ignore[arg-type]
    errors: list[str] = []
    controller._show_error = errors.append

    controller._start_analysis()

    assert errors
    assert controller.model.report is before.report
    assert controller.model.selection is before.selection
    assert controller._report_request is current_request
    controller.close()
    panel.deleteLater()


def test_snapshot_capture_failure_is_user_visible_and_does_not_escape_compare(
    application, tmp_path: Path, monkeypatch
) -> None:
    """Filesystem races at capture are contained at the GUI request boundary."""

    current_request, report = rich_report(tmp_path)
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    assert controller._load_source("reference", Path(current_request.reference_selection.path))
    assert controller._load_source("target", Path(current_request.target_selection.path))
    controller._report_request = current_request
    controller._analysis_finished(report)
    errors: list[str] = []
    controller._show_error = errors.append

    def race(*_args: object, **_kwargs: object) -> object:
        raise OSError("source disappeared during capture")

    monkeypatch.setattr("structlens.plugin.gui.report_controller.capture_snapshot", race)
    controller._start_analysis()

    assert errors
    assert "capture" in errors[-1].lower() or "source" in errors[-1].lower()
    assert controller.model.report is report
    controller.close()
    panel.deleteLater()


class _ImmediateFuture:
    def __init__(self, value: object) -> None:
        self._value = value

    def done(self) -> bool:
        return True

    def result(self) -> object:
        return self._value

    def cancel(self) -> bool:
        return True


class _ImmediateExecutor:
    def submit(self, function: object, built_request: object) -> _ImmediateFuture:
        assert callable(function)
        return _ImmediateFuture(function(built_request))


__all__ = ["test_snapshot_request_builder_preserves_exact_non_first_selection_and_evidence"]
