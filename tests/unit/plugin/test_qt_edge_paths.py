"""Behavioral regressions for Qt source, runner, export, and view boundaries."""

from __future__ import annotations

from pathlib import Path

import pytest
from PySide6.QtWidgets import QApplication

from structlens.core.errors import AnalysisCancelledError
from structlens.core.models import AnalysisResult
from structlens.plugin.gui.main_panel import build_qt_panel
from tests.unit.plugin._task13_round4_fixtures import rich_report
from tests.unit.plugin.test_gui_model import _minimal_report

FIXTURE = Path("tests/fixtures/parsing/numbering_altloc.pdb")


@pytest.fixture(scope="module")
def application():
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


def _legacy_result() -> AnalysisResult:
    return AnalysisResult(
        reference_id="reference",
        target_id="target",
        correspondences=(),
        mutations=(),
        sequence_identity=1.0,
        sequence_coverage=1.0,
        alignment_decision="accepted",
    )


class _Future:
    def __init__(self, value: object = None, *, done: bool = True, cancel: bool = False):
        self.value = value
        self.done_value = done
        self.cancel_value = cancel
        self.cancel_called = False

    def done(self) -> bool:
        return self.done_value

    def result(self) -> object:
        if isinstance(self.value, BaseException):
            raise self.value
        return self.value

    def cancel(self) -> bool:
        self.cancel_called = True
        return self.cancel_value


class _PyMol:
    def __init__(self, *, names=(), content: str | None = None, error: Exception | None = None):
        self.names = list(names)
        self.content = content or FIXTURE.read_text(encoding="utf-8")
        self.error = error
        self.saved: list[tuple[str, str]] = []
        self.loaded: list[tuple[str, str]] = []
        self.deleted: list[str] = []

    def get_names(self, kind: str) -> list[str]:
        assert kind == "objects"
        return list(self.names)

    def save(self, path: str, name: str) -> None:
        if self.error is not None:
            raise self.error
        self.saved.append((path, name))
        Path(path).write_text(self.content, encoding="utf-8")

    def load(self, path: str, name: str) -> None:
        if self.error is not None:
            raise self.error
        self.loaded.append((path, name))
        if name not in self.names:
            self.names.append(name)

    def delete(self, name: str) -> None:
        if self.error is not None:
            raise self.error
        self.deleted.append(name)
        self.names = [value for value in self.names if value != name]


def test_browse_cancel_does_not_replace_a_source(controller, monkeypatch) -> None:
    monkeypatch.setattr(controller.w.QFileDialog, "getOpenFileName", lambda *_: ("", ""))

    controller._browse_source("reference")

    assert controller.reference_structure is None
    assert controller.reference_edit.text() == ""


def test_browse_success_loads_structure_and_exposes_model_zero(controller, monkeypatch) -> None:
    monkeypatch.setattr(controller.w.QFileDialog, "getOpenFileName", lambda *_: (str(FIXTURE), "PDB"))

    controller._browse_source("reference")

    assert controller.reference_structure is not None
    assert controller.reference_edit.text() == str(FIXTURE)
    assert controller.reference_model_combo.itemData(0) == "0"
    assert controller.reference_chain_combo.count() == 1


def test_empty_coordinate_parse_clears_stale_source_and_reports_invalid_input(controller, tmp_path: Path) -> None:
    assert controller._load_source("reference", FIXTURE) is True
    invalid = tmp_path / "empty-records.pdb"
    invalid.write_text("remark without coordinate records\n", encoding="utf-8")
    errors: list[str] = []
    controller._show_error = errors.append

    loaded = controller._load_source("reference", invalid)

    assert loaded is False
    assert controller.reference_structure is None
    assert controller.reference_chain_combo.count() == 0
    assert errors and "Could not load reference" in errors[-1]


def test_model_change_repopulates_chains_for_the_selected_pdb_model(controller, tmp_path: Path) -> None:
    records = [line for line in FIXTURE.read_text(encoding="utf-8").splitlines() if line != "END"]
    multi_model = tmp_path / "two-models.pdb"
    multi_model.write_text(
        "MODEL        1\n" + "\n".join(records) + "\nENDMDL\nMODEL        2\n" + "\n".join(records) + "\nENDMDL\nEND\n",
        encoding="utf-8",
    )

    assert controller._load_source("reference", multi_model) is True
    controller.reference_model_combo.setCurrentIndex(1)
    controller._model_changed("reference")

    selected = controller._selected_chain(
        controller.reference_structure,
        controller.reference_chain_combo,
        controller.reference_model_combo,
    )
    assert selected is not None
    assert selected.model_id == "2"


@pytest.mark.parametrize(
    ("command", "message"),
    [
        (None, "available only inside a PyMOL session"),
        (object(), "cannot list and save objects"),
        (_PyMol(), "No PyMOL objects are available"),
    ],
)
def test_pymol_source_guards_are_user_visible(controller, monkeypatch, command, message) -> None:
    controller.command = command
    errors: list[str] = []
    monkeypatch.setattr(controller, "_show_error", errors.append)

    controller._use_pymol_object("reference")

    assert errors and message in errors[-1]


def test_pymol_object_cancel_does_not_save_or_change_source(controller, monkeypatch) -> None:
    command = _PyMol(names=("ligand_complex",))
    controller.command = command
    monkeypatch.setattr(controller.w.QInputDialog, "getItem", lambda *_: ("ligand_complex", False))

    controller._use_pymol_object("reference")

    assert command.saved == []
    assert controller.reference_structure is None


def test_pymol_object_save_failure_cleans_up_and_reports(controller, monkeypatch) -> None:
    command = _PyMol(names=("broken",), error=OSError("save failed"))
    controller.command = command
    monkeypatch.setattr(controller.w.QInputDialog, "getItem", lambda *_: ("broken", True))
    errors: list[str] = []
    monkeypatch.setattr(controller, "_show_error", errors.append)

    controller._use_pymol_object("reference")

    assert errors and "Could not read PyMOL object broken" in errors[-1]
    assert controller.reference_structure is None
    assert controller._temporary_paths == []


def test_named_pymol_object_success_uses_the_same_source_loader(controller) -> None:
    command = _PyMol()
    controller.command = command

    loaded = controller._load_named_pymol_object("target", "complex")

    assert loaded is True
    assert controller.target_structure is not None
    assert controller.target_object_name == "complex"
    assert command.saved and command.saved[-1][1] == "complex"


def test_pymol_object_success_is_parsed_and_materializes_a_snapshot_owned_object(controller, monkeypatch) -> None:
    command = _PyMol(names=("complex",))
    panel = controller.widget
    panel._structlens_controller.command = command
    monkeypatch.setattr(controller.w.QInputDialog, "getItem", lambda *_: ("complex", True))

    controller._use_pymol_object("reference")

    assert controller.reference_structure is not None
    assert command.saved and command.saved[-1][1] == "complex"
    assert controller.reference_object_name is not None
    assert controller.reference_object_name != "complex"
    assert controller.reference_object_name.startswith("structlens_panel_")
    assert command.loaded and command.loaded[-1][1] == controller.reference_object_name


@pytest.mark.parametrize("operation", ["missing-save", "load-error", "already-loaded"])
def test_pymol_file_boundary_handles_missing_or_existing_objects(controller, operation) -> None:
    if operation == "missing-save":
        controller.command = object()
        assert controller._load_named_pymol_object("reference", "obj") is False
        return

    if operation == "load-error":
        controller.command = _PyMol(error=OSError("load failed"))
        assert controller._load_file_into_pymol("reference", FIXTURE) is None
        return

    command = _PyMol(names=("structlens_panel_numbering_altloc_reference",))
    controller.command = command
    assert controller._load_file_into_pymol("reference", FIXTURE) == "structlens_panel_numbering_altloc_reference"
    assert command.deleted == ["structlens_panel_numbering_altloc_reference"]
    assert command.loaded and command.loaded[-1][1] == "structlens_panel_numbering_altloc_reference"


def test_manual_pair_parser_accepts_comments_and_exact_residue_identity(controller) -> None:
    assert controller._load_source("reference", FIXTURE) is True
    assert controller._load_source("target", FIXTURE) is True
    reference = controller._selected_chain(controller.reference_structure, controller.reference_chain_combo)
    target = controller._selected_chain(controller.target_structure, controller.target_chain_combo)
    assert reference is not None and target is not None
    residue = reference.residue_records[0].residue_id
    token = f"{residue.chain_id}:{residue.auth_seq_id}{residue.insertion_code or ''}"
    controller.manual_edit.setPlainText(f"# authoritative pair\n{token} -> {token}\n")

    pairs = controller._manual_pairs(reference, target)

    assert pairs == [(residue, target.residue_records[0].residue_id)]


@pytest.mark.parametrize(
    ("text", "message"),
    [("A:100", "must use ref"), ("A:999 -> A:999", "names a residue")],
)
def test_manual_pair_parser_rejects_malformed_or_unknown_residues(controller, text, message) -> None:
    assert controller._load_source("reference", FIXTURE) is True
    assert controller._load_source("target", FIXTURE) is True
    reference = controller._selected_chain(controller.reference_structure, controller.reference_chain_combo)
    target = controller._selected_chain(controller.target_structure, controller.target_chain_combo)
    assert reference is not None and target is not None
    controller.manual_edit.setPlainText(text)

    with pytest.raises(ValueError, match=message):
        controller._manual_pairs(reference, target)


def test_load_sources_from_edits_requires_both_coordinate_sources(controller) -> None:
    controller.reference_edit.setText(str(FIXTURE))
    controller.target_edit.clear()

    loaded = controller._load_sources_from_edits()

    assert loaded is False
    assert controller.reference_structure is not None
    assert controller.target_structure is None
    assert controller.target_chain_combo.count() == 0


def test_start_analysis_reports_invalid_request_without_submitting_a_future(controller, monkeypatch) -> None:
    assert controller._load_source("reference", FIXTURE) is True
    assert controller._load_source("target", FIXTURE) is True
    errors: list[str] = []
    monkeypatch.setattr(controller, "_show_error", errors.append)
    monkeypatch.setattr(
        controller._report_controller,
        "build_request",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(ValueError("invalid chain request")),
    )

    controller._start_analysis()

    assert errors == ["invalid chain request"]
    assert controller._future is None
    assert controller.pages.currentWidget().objectName() == "pageStructures"


def test_poll_pending_future_does_not_publish_a_partial_report(controller) -> None:
    pending = _Future(done=False)
    controller._future = pending
    controller.model = controller.model.with_busy("Building canonical analysis report…")

    controller._poll_analysis()

    assert controller._future is pending
    assert controller.model.report is None
    assert controller.model.busy is True


@pytest.mark.parametrize(
    ("value", "status"),
    [(_minimal_report(), "Analysis report complete"), (AnalysisCancelledError(), "Comparison cancelled")],
)
def test_poll_finished_future_has_truthful_report_or_cancel_state(controller, value, status) -> None:
    controller._future = _Future(value)

    controller._poll_analysis()

    assert controller._future is None
    assert status in controller.footer_status.text()
    if isinstance(value, AnalysisCancelledError):
        assert controller.model.report is None
    else:
        assert controller.model.report is value


def test_poll_failed_future_surfaces_error_and_clears_busy_state(controller, monkeypatch) -> None:
    errors: list[str] = []
    monkeypatch.setattr(controller, "_show_error", errors.append)
    controller._future = _Future(RuntimeError("service unavailable"))
    controller.model = controller.model.with_busy("Building canonical analysis report…")

    controller._poll_analysis()

    assert errors == ["Comparison failed: service unavailable"]
    assert controller.model.busy is False


def test_cancel_running_future_reports_waiting_instead_of_false_completion(controller) -> None:
    running = _Future(done=False, cancel=False)
    controller._future = running
    controller.model = controller.model.with_busy("Building canonical analysis report…")

    controller._cancel_analysis()

    assert running.cancel_called is True
    assert controller._future is running
    assert "Cancellation requested; waiting" in controller.footer_status.text()
    assert controller.model.report is None


def test_cancel_queued_future_publishes_explicit_cancel_state(controller) -> None:
    queued = _Future(done=False, cancel=True)
    controller._future = queued

    controller._cancel_analysis()

    assert controller._future is None
    assert controller.model.status == "Comparison cancelled."
    assert "no result was changed" in controller.footer_status.text()


def test_canonical_report_wins_over_stale_legacy_xlsx_records(controller, monkeypatch) -> None:
    report = _minimal_report()
    controller._populate_report(report)
    controller.set_v03_export_records(provenance=("stale",))
    routed: list[tuple[str, str]] = []
    errors: list[str] = []
    monkeypatch.setattr(controller, "_export", lambda suffix, _exporter, label: routed.append((suffix, label)))
    monkeypatch.setattr(controller, "_show_error", errors.append)

    controller._export_xlsx()

    assert routed == [("xlsx", "XLSX")]
    assert errors == []


def test_export_cancel_does_not_call_exporter(controller, monkeypatch) -> None:
    controller.model = controller.model.with_report(_minimal_report())
    calls: list[object] = []
    monkeypatch.setattr(controller.w.QFileDialog, "getSaveFileName", lambda *_: ("", ""))

    controller._export("json", lambda *_args, **_kwargs: calls.append(True), "JSON")

    assert calls == []


def test_exporter_failure_is_reported_without_claiming_success(controller, monkeypatch, tmp_path: Path) -> None:
    """A render-only report has no verified evidence provenance to export."""

    controller.model = controller.model.with_report(_minimal_report())
    output = tmp_path / "failed.json"
    monkeypatch.setattr(controller.w.QFileDialog, "getSaveFileName", lambda *_: (str(output), "JSON"))
    errors: list[str] = []
    monkeypatch.setattr(controller, "_show_error", errors.append)
    calls: list[object] = []

    def failing_exporter(*_args, **_kwargs):
        calls.append(True)
        raise OSError("disk full")

    controller._export("json", failing_exporter, "JSON")

    assert errors == ["Could not export JSON: canonical report has no verified source binding"]
    assert calls == []
    assert not output.exists()
    assert "JSON export written" not in controller.footer_status.text()


def test_verified_binding_exporter_failure_is_reported_without_claiming_success(
    controller, monkeypatch, tmp_path: Path
) -> None:
    """Once provenance is verified, exporter failures still remain truthful."""

    request, report = rich_report(tmp_path)
    assert controller._load_source("reference", Path(request.reference_selection.path))
    assert controller._load_source("target", Path(request.target_selection.path))
    controller._report_request = request
    controller._analysis_finished(report)
    output = tmp_path / "failed-verified.json"
    monkeypatch.setattr(controller.w.QFileDialog, "getSaveFileName", lambda *_: (str(output), "JSON"))
    errors: list[str] = []
    monkeypatch.setattr(controller, "_show_error", errors.append)
    calls: list[object] = []

    def failing_exporter(*_args, **_kwargs):
        calls.append(True)
        raise OSError("disk full")

    controller._export("json", failing_exporter, "JSON")

    assert errors == ["Could not export JSON: disk full"]
    assert calls == [True]
    assert not output.exists()
    assert "JSON export written" not in controller.footer_status.text()


def test_save_project_failure_is_user_visible(controller, monkeypatch, tmp_path: Path) -> None:
    output = tmp_path / "project.json"
    monkeypatch.setattr(controller.w.QFileDialog, "getSaveFileName", lambda *_: (str(output), "JSON"))
    monkeypatch.setattr(
        "structlens.plugin.gui.qt_exports_mixin.ProjectState.save",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(OSError("read-only directory")),
    )
    errors: list[str] = []
    monkeypatch.setattr(controller, "_show_error", errors.append)

    controller._save_project()

    assert errors == ["Could not save project: read-only directory"]
    assert not output.exists()


def test_open_project_failure_does_not_replace_current_state(controller, monkeypatch, tmp_path: Path) -> None:
    current = _legacy_result()
    controller.model = controller.model.with_analysis(current)
    missing = tmp_path / "missing.structlens.json"
    monkeypatch.setattr(controller.w.QFileDialog, "getOpenFileName", lambda *_: (str(missing), "JSON"))
    errors: list[str] = []
    monkeypatch.setattr(controller, "_show_error", errors.append)

    controller._open_project()

    assert errors and "Could not open project" in errors[-1]
    assert controller.model.analysis is current


def test_pymol_bundle_requires_report_and_both_coordinate_sources(controller, monkeypatch) -> None:
    errors: list[str] = []
    monkeypatch.setattr(controller, "_show_error", errors.append)
    controller._export_pymol_bundle()
    assert errors[-1] == "Run a comparison before exporting a PyMOL bundle."

    controller._populate_report(_minimal_report())
    controller._export_pymol_bundle()
    assert errors[-1] == "Run a comparison before exporting a PyMOL bundle."


def test_chart_export_cancel_keeps_authoritative_dataset_untouched(controller, monkeypatch) -> None:
    from structlens.application.chart_data import ChartDataset

    dataset = ChartDataset("msa", "MSA", "column", "conservation", "fraction", (), "descriptive")
    controller.model = controller.model.with_report(_minimal_report())
    controller.set_chart_datasets({"MSA conservation profile": dataset})
    controller.chart_combo.setCurrentText("MSA conservation profile")
    monkeypatch.setattr(controller.w.QFileDialog, "getSaveFileName", lambda *_: ("", ""))
    calls: list[object] = []

    monkeypatch.setattr(
        "structlens.plugin.gui.qt_panel.export_chart_xlsx",
        lambda *_args, **_kwargs: calls.append(True),
    )
    controller._export_chart_xlsx()

    assert calls == []
    assert controller._chart_datasets["MSA conservation profile"] is dataset


def test_chart_export_reports_unavailable_profile_without_opening_dialog(controller, monkeypatch) -> None:
    controller.model = controller.model.with_report(_minimal_report())
    controller.chart_combo.setCurrentText("Pairwise similarity heatmap")
    errors: list[str] = []
    monkeypatch.setattr(controller, "_show_error", errors.append)

    controller._export_chart_xlsx()

    assert errors == [
        "The Pairwise similarity heatmap dataset is unavailable. Run its scientific service before exporting."
    ]


def test_visualization_requires_report_and_pymol_objects(controller, monkeypatch) -> None:
    errors: list[str] = []
    monkeypatch.setattr(controller, "_show_error", errors.append)
    controller._apply_visualization()
    assert errors[-1] == "Run a comparison before applying a PyMOL view."

    controller.model = controller.model.with_analysis(_legacy_result())
    controller.command = object()
    controller._apply_visualization()
    assert "not loaded as PyMOL objects" in errors[-1]

    controller.command = None
    controller._apply_visualization()
    assert errors[-1] == "Apply a visualization from inside a PyMOL session."


def test_visualization_empty_state_has_explicit_legend_and_reset_status(controller) -> None:
    controller._update_legend()
    assert controller.legend_label.text() == "Legend appears after a comparison."
    assert controller.visualization_count.text() == "0 rows selected"

    controller._reset_visualization()
    assert controller.footer_status.text() == "StructLens-owned PyMOL selections removed"


def test_visualization_success_uses_the_report_compatible_adapter(controller, monkeypatch) -> None:
    controller.model = controller.model.with_analysis(_legacy_result())
    controller.command = object()
    controller.reference_object_name = "reference_obj"
    controller.target_object_name = "target_obj"
    calls: list[dict[str, object]] = []
    monkeypatch.setattr(controller._pymol, "apply", lambda result, **kwargs: calls.append({"result": result, **kwargs}))

    controller._apply_visualization()

    assert calls and calls[0]["result"] is controller.model.analysis
    assert calls[0]["reference_object"] == "reference_obj"
    assert calls[0]["target_object"] == "target_obj"
    assert "PyMOL view applied" in controller.footer_status.text()


def test_visualization_preset_is_reflected_in_immutable_panel_state(controller) -> None:
    before = controller.model.visualization_state

    controller._apply_preset("Structural deviation")

    assert controller.model.visualization_state != before
    assert controller.model.visualization_state.preset == "Structural deviation"
    assert controller.filter_combo.currentData() == "displacement"
    assert controller.color_combo.currentData() == "ca_displacement"
