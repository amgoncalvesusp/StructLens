"""Focused edge coverage for the report-first GUI adapters.

These tests exercise error boundaries and compatibility seams directly.  They
do not assert implementation details of the scientific services; every
scientific value is obtained from the existing typed report fixtures.
"""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pytest
from PySide6.QtWidgets import QApplication, QTableWidget

from structlens.application.project_state import ProjectState
from structlens.core.evidence import Availability
from structlens.core.models import AnalysisResult, ResidueId
from structlens.plugin.gui import (
    alignment_page,
    mutations_page,
    project_page,
    residues_page,
    results_page,
    visualization_page,
    widgets,
)
from structlens.plugin.gui import model as model_module
from structlens.plugin.gui import project_transaction as transaction
from structlens.plugin.gui.main_panel import build_qt_panel
from structlens.plugin.gui.model import CanonicalReportBinding, StructLensPanelModel
from structlens.plugin.gui.qt_chart import load_chart_classes, matrix_image_kwargs
from structlens.plugin.gui.qt_legacy import legacy_site_backend
from structlens.plugin.gui.qt_pages import (
    backend_label,
    configure_table,
    fraction_number,
    fraction_spin,
    human,
    label,
    number,
    page_subtitle,
    residue_label,
    state_dict,
    stylesheet,
)
from structlens.plugin.gui.qt_reports import legacy_analysis, report_chart_datasets, report_site_metrics
from structlens.plugin.gui.qt_sources import combo_data, find_residue, selected_chain, structure_meta
from structlens.plugin.visualization.renderer import VisualizationState

from ._task13_round4_fixtures import rich_report


@pytest.fixture(scope="module")
def application() -> QApplication:
    app = QApplication.instance() or QApplication([])
    yield app
    app.quit()


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


def test_compatibility_page_descriptors_are_importable() -> None:
    assert alignment_page.PAGE.name == "Structures"
    assert mutations_page.PAGE.name == "Mutations"
    assert project_page.PAGE.name == "Project"
    assert residues_page.PAGE.name == "Residues"
    assert results_page.PAGE.name == "Results"
    assert visualization_page.PAGE.name == "Charts"
    assert widgets.PageDescriptor is type(alignment_page.PAGE)


def test_model_binding_rejects_bad_types_mismatched_inputs_and_projection(tmp_path: Path, monkeypatch) -> None:
    request, report = rich_report(tmp_path)

    with pytest.raises(TypeError, match="report"):
        CanonicalReportBinding.create(object(), request)  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="request"):
        CanonicalReportBinding.create(report, object())  # type: ignore[arg-type]

    monkeypatch.setattr(model_module, "verify_report_inputs", lambda *_args: (_ for _ in ()).throw(ValueError("stale")))
    with pytest.raises(ValueError, match="stale"):
        CanonicalReportBinding.create(report, request)

    monkeypatch.undo()
    other_selection = replace(request.reference_selection, model_id="different")
    stale_request = replace(request, reference_selection=other_selection)
    with pytest.raises(ValueError, match="reference request selection"):
        CanonicalReportBinding.create(report, stale_request)
    other_selection = replace(request.target_selection, model_id="different")
    stale_request = replace(request, target_selection=other_selection)
    with pytest.raises(ValueError, match="target request selection"):
        CanonicalReportBinding.create(report, stale_request)

    presentation = model_module.present_report(report)
    with pytest.raises(ValueError, match="presentation"):
        CanonicalReportBinding.create(report, request, replace(presentation, report_id="wrong"))


def test_panel_model_transitions_are_strict_and_replace_values(tmp_path: Path) -> None:
    request, report = rich_report(tmp_path)
    binding = CanonicalReportBinding.create(report, request)
    legacy = _legacy_result()
    state = StructLensPanelModel().with_analysis(legacy)
    assert state.analysis is legacy and state.report is None and state.binding is None
    assert state.with_status("ready").status == "ready"
    assert state.with_sources(reference_path="a", target_path="b").reference_path == "a"
    assert state.with_visualization(VisualizationState()).visualization_state == VisualizationState()
    assert state.with_busy("running").busy is True
    assert state.with_error("bad").error == "bad"
    with pytest.raises(TypeError, match="report"):
        state.with_report(object())  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="binding"):
        state.with_binding(object())  # type: ignore[arg-type]
    canonical = state.with_binding(binding)
    assert canonical.report is report and canonical.analysis is None and canonical.selection is None


def test_project_transaction_validation_edges(tmp_path: Path, monkeypatch) -> None:
    request, report = rich_report(tmp_path)
    binding = CanonicalReportBinding.create(report, request)
    output = tmp_path / "edge.structlens.json"
    transaction.save_canonical_project(binding, output, visualization_state={})
    artifact_root = Path(f"{output}.artifacts") / report.report_id

    # Existing verified generations are checked instead of overwritten.
    transaction.save_canonical_project(binding, output, visualization_state={})
    assert transaction.stage_project(output).binding.report.report_id == report.report_id

    with pytest.raises(TypeError, match="verified canonical"):
        transaction.save_canonical_project(object(), output, visualization_state={})  # type: ignore[arg-type]

    for payload in (
        {"report_hash": report.report_id, "report_verification": "unverified"},
        {"report_hash": report.report_id, "report_verification": "verified"},
    ):
        ProjectState(**payload).save(tmp_path / "incomplete.json")  # type: ignore[arg-type]
        with pytest.raises(Exception, match="canonical"):
            transaction.stage_project(tmp_path / "incomplete.json")

    project_payload = json.loads(output.read_text(encoding="utf-8"))
    references = project_payload["v04"]["report_references"]
    project_payload["v04"]["report_references"] = list(reversed(references)) + [references[0]]
    output.write_text(json.dumps(project_payload), encoding="utf-8")
    with pytest.raises(Exception, match="references"):
        transaction.stage_project(output)

    # Restore a valid project and check the remaining path guards directly.
    transaction.save_canonical_project(binding, output, visualization_state={})
    with pytest.raises(transaction.ProjectSchemaError, match="non-empty"):
        transaction._validated_reference(output, artifact_root, "")
    with pytest.raises(transaction.ProjectSchemaError, match="escapes"):
        transaction._validated_reference(output, artifact_root, "../report.json")
    with pytest.raises(transaction.ProjectSchemaError, match="escapes"):
        transaction._validated_reference(output, artifact_root, str(output))
    with pytest.raises(transaction.ProjectSchemaError, match="missing or unsafe"):
        transaction._validate_generation_root(tmp_path / "missing-generation")
    artifact_file = artifact_root / "report.json"
    with pytest.raises(transaction.ProjectSchemaError, match="missing or unsafe"):
        transaction._validate_regular_artifact(artifact_file.with_name("missing"))

    monkeypatch.setattr(Path, "lstat", lambda _path: (_ for _ in ()).throw(OSError("unreadable")))
    assert transaction._is_link_or_reparse(tmp_path / "anything") is True


def test_project_transaction_private_metadata_and_restore_failures(tmp_path: Path) -> None:
    request, report = rich_report(tmp_path)
    project = ProjectState(
        reference_selection=request.reference_selection,
        target_selection=request.target_selection,
        evidence_sources={"canonical_request": {}},
        source_hashes={"reference_raw_sha256": request.reference_snapshot.raw_sha256},
    )
    assert transaction._has_canonical_metadata(project)
    assert transaction._snapshot_directory(tmp_path, "reference") == tmp_path / "snapshots"
    assert transaction._snapshot_directory(tmp_path, "target") == tmp_path / "snapshots" / "target"
    with pytest.raises(transaction.ProjectSchemaError, match="target_raw"):
        transaction._source_hash(project, "target_raw_sha256")
    assert transaction._source_hash(project, "reference_raw_sha256") == request.reference_snapshot.raw_sha256

    with pytest.raises(transaction.ProjectSchemaError, match="selections"):
        transaction._restore_request(ProjectState(), (request.reference_snapshot, request.target_snapshot))
    invalid_metadata = replace(project, evidence_sources={"canonical_request": []})
    with pytest.raises(transaction.ProjectSchemaError, match="metadata"):
        transaction._restore_request(invalid_metadata, (request.reference_snapshot, request.target_snapshot))
    invalid_manual = replace(project, evidence_sources={"canonical_request": {"manual_pairs": "bad"}})
    with pytest.raises(transaction.ProjectSchemaError, match="manual pairs"):
        transaction._restore_request(invalid_manual, (request.reference_snapshot, request.target_snapshot))
    with pytest.raises(transaction.ProjectSchemaError, match="manual pair"):
        transaction._manual_pair({"reference": {}})
    with pytest.raises(transaction.ProjectSchemaError, match="manual pair"):
        transaction._manual_pair("bad")
    with pytest.raises(transaction.ProjectSchemaError, match="site definition"):
        transaction._site_from_payload("bad")
    with pytest.raises(transaction.ProjectSchemaError, match="site definition"):
        transaction._site_from_payload({"site_id": "bad"})


def test_qt_chart_pages_and_report_compatibility_helpers(application, tmp_path: Path, monkeypatch) -> None:
    request, report = rich_report(tmp_path)
    matrix = SimpleNamespace(
        chart_id="distance_difference", title="distance", row_label="r", column_label="c", interpretation="delta"
    )
    assert matrix_image_kwargs(matrix, []) == {"aspect": "auto", "interpolation": "nearest"}
    assert matrix_image_kwargs(matrix, (-2.0, None, 3.0))["vmin"] == -3.0
    plain_matrix = SimpleNamespace(
        chart_id="ordinary",
        title="ordinary",
        row_label="r",
        column_label="c",
        interpretation="value",
    )
    assert matrix_image_kwargs(plain_matrix, (1.0, 1.0)) == {
        "aspect": "auto",
        "interpolation": "nearest",
    }
    classes = load_chart_classes()
    assert classes is None or len(classes) == 2

    widgets = SimpleNamespace(
        QLabel=lambda text: SimpleNamespace(text=text, setObjectName=lambda name: None),
        QPushButton=lambda text: SimpleNamespace(text=text, setObjectName=lambda name: None),
    )
    assert label(widgets, "x", "y").text == "x"
    assert (button := widgets.QPushButton("ok"))
    assert button.text == "ok"
    assert page_subtitle("Charts") != page_subtitle("Charts", standalone=True)
    assert fraction_number(None) == "—" and fraction_number(0.5) == "0.500"
    assert number(None) == "—" and number(1.25) == "1.250" and number(2) == "2"
    assert residue_label(None) == "—"
    assert residue_label(ResidueId("id", "1", "A", "1", None, "ALA")) == "A:1 ALA"
    assert human("ca_rmsd") == "Cα Rmsd"
    assert stylesheet()
    assert state_dict(VisualizationState())["representation"]

    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    spin = controller.coverage_spin
    fraction_spin(spin, 0.4)
    assert spin.value() == pytest.approx(0.4)
    table = QTableWidget()
    configure_table(table, controller.w)
    assert table.minimumHeight() == 220
    assert backend_label(_legacy_result()) == "sequence-guided"
    assert backend_label(replace(_legacy_result(), provenance={"backend": "fixture"})) == "fixture"
    assert report_site_metrics(report) == tuple(report.sites)
    assert report_chart_datasets(report)
    assert legacy_analysis(report) is not None
    no_analysis = replace(
        report, analysis=None, availability=replace(report.availability, analysis=Availability.NOT_DETECTED)
    )
    assert legacy_analysis(no_analysis) is None
    assert legacy_analysis(None) is None
    with pytest.raises(RuntimeError, match="unavailable"):
        legacy_site_backend()
    controller.close()
    panel.deleteLater()
    monkeypatch.setattr("structlens.plugin.gui.qt_compat.load_qt", lambda: None)
    with pytest.raises(RuntimeError, match="requires PySide6"):
        build_qt_panel()


def test_qt_chart_and_report_optional_paths_are_explicit(tmp_path: Path, monkeypatch) -> None:
    _, report = rich_report(tmp_path)
    delta_zero = SimpleNamespace(
        chart_id="distance_difference",
        title="distance",
        row_label="r",
        column_label="c",
        interpretation="delta",
    )
    assert matrix_image_kwargs(delta_zero, (0.0, 0.0)) == {
        "aspect": "auto",
        "interpolation": "nearest",
    }

    import structlens.plugin.gui.qt_chart as chart_module

    def unavailable_import(_name: str) -> object:
        raise ImportError("matplotlib not installed")

    monkeypatch.setattr(chart_module.importlib, "import_module", unavailable_import)
    assert chart_module.load_chart_classes() is None

    no_msa = replace(report, msa=None, availability=replace(report.availability, msa=Availability.NOT_DETECTED))
    assert "MSA conservation profile" not in report_chart_datasets(no_msa)
    no_analysis = replace(
        report,
        analysis=None,
        availability=replace(report.availability, analysis=Availability.NOT_DETECTED),
    )
    no_analysis_datasets = report_chart_datasets(no_analysis)
    assert "Structural deviation profile" not in no_analysis_datasets
    assert "MSA conservation profile" in no_analysis_datasets
    no_sites = replace(report, sites=None, availability=replace(report.availability, sites=Availability.NOT_DETECTED))
    assert report_site_metrics(no_sites) == ()


def test_report_mixin_compatibility_and_focus_failure_paths(application, tmp_path: Path, monkeypatch) -> None:
    request, report = rich_report(tmp_path)
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    assert controller._load_source("reference", Path(request.reference_selection.path))
    assert controller._load_source("target", Path(request.target_selection.path))
    controller._report_request = request
    controller._analysis_finished(report)
    controller._fill_report_structure_result(report)
    controller._fill_report_mutations(report)
    controller._fill_report_residues(report)
    controller._fill_report_history()
    controller._fill_report_evidence_label(report)
    assert "Reference" in controller.evidence_card_label.text()

    # A report with no cards keeps its exact unavailable state visible.
    no_cards = replace(
        report,
        evidence_cards=(),
        availability=replace(report.availability, evidence_cards=Availability.NOT_DETECTED),
    )
    controller._fill_report_evidence_label(no_cards)
    assert "unavailable" in controller.evidence_card_label.text().lower()
    controller._focus_report_alignment_index(0)
    errors: list[str] = []
    monkeypatch.setattr(controller, "_show_error", errors.append)
    controller._focus_report_alignment_index(-1)
    assert errors

    # Empty typed report sections take the safe early-return paths.
    controller._fill_report_structure_result(no_cards)
    controller._fill_report_mutations(no_cards)
    controller._fill_report_residues(no_cards)
    controller._fill_report_history()
    controller._focus_mutation(999, 0)
    controller._focus_residue(999, 0)
    controller._analysis_finished(_legacy_result())
    controller._focus_mutation(999, 0)
    controller._focus_residue(999, 0)
    controller.close()
    panel.deleteLater()


def test_report_mixin_legacy_rows_and_canonical_missing_sections(application, tmp_path: Path, monkeypatch) -> None:
    request, report = rich_report(tmp_path)
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    legacy = legacy_analysis(report)
    assert legacy is not None and legacy.mutations and legacy.correspondences

    # The compatibility tables still render a legacy AnalysisResult, including
    # mutation and correspondence rows, when a pre-v0.4 caller supplies one.
    controller.model = controller.model.with_analysis(legacy)
    controller._fill_mutations(legacy)
    controller._fill_residues(legacy)
    assert controller.mutation_table.rowCount() == len(legacy.mutations)
    assert controller.residue_table.rowCount() == len(legacy.correspondences)
    controller._focus_mutation(0, 0)
    controller._focus_residue(0, 0)
    controller._focus_alignment_index(legacy, legacy.correspondences[0].alignment_index)
    assert "selected" in controller.footer_status.text().lower()

    no_analysis = replace(
        report,
        analysis=None,
        availability=replace(report.availability, analysis=Availability.NOT_DETECTED),
    )
    controller._report_request = request
    monkeypatch.setattr(
        "structlens.plugin.gui.qt_reports_mixin.CanonicalReportBinding.create",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(ValueError("incompatible binding")),
    )
    controller._populate_report(no_analysis)
    assert controller.model.report is no_analysis
    assert controller.model.binding is None
    controller._fill_report_structure_result(no_analysis)
    controller._report_history = (no_analysis,)
    controller._fill_report_history()
    controller._update_report_evidence_label(no_analysis)
    controller._focus_report_alignment_index(0)

    # A valid report without evidence cards still permits correspondence focus;
    # the missing card is explicitly treated as unavailable.
    no_cards = replace(
        report,
        evidence_cards=(),
        availability=replace(report.availability, evidence_cards=Availability.NOT_DETECTED),
    )
    controller._populate_report(no_cards)
    controller._focus_report_alignment_index(report.analysis.correspondences[0].alignment_index)
    assert "unavailable" in controller.evidence_card_label.text().lower()
    controller.close()
    panel.deleteLater()


def test_qt_source_helpers_cover_empty_and_fallback_paths() -> None:
    class Combo:
        def __init__(self, value=None):
            self.value = value

        def currentData(self):
            return self.value

    assert combo_data(Combo(None)) is None
    assert combo_data(Combo(4)) == "4"
    assert selected_chain(None, Combo("A")) is None
    assert find_residue(SimpleNamespace(residue_records=(), residues=()), "malformed") is None
    assert structure_meta(SimpleNamespace(chains=(), structure_id="empty")) == "0 chain(s) · 0 residues · empty"


def test_export_and_open_cancel_paths_are_noops(application, monkeypatch) -> None:
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    errors: list[str] = []
    monkeypatch.setattr(controller, "_show_error", errors.append)
    monkeypatch.setattr(controller.w.QFileDialog, "getSaveFileName", lambda *_: ("", ""))
    controller._save_project()
    monkeypatch.setattr(controller.w.QFileDialog, "getOpenFileName", lambda *_: ("", ""))
    controller._open_project()
    controller._export("json", lambda *_args, **_kwargs: pytest.fail("exporter must not run"), "JSON")
    assert errors and "comparison" in errors[-1].lower()
    controller.close()
    panel.deleteLater()


def test_legacy_project_commit_source_and_object_failures(application, monkeypatch, tmp_path: Path) -> None:
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    fixture = Path(__file__).parents[2] / "fixtures" / "parsing" / "numbering_altloc.pdb"
    legacy = ProjectState(reference_source=str(fixture), target_sources=(str(fixture),))
    controller._commit_legacy_project(legacy, tmp_path / "legacy.json")
    assert controller.reference_structure is not None and controller.target_structure is not None

    calls: list[tuple[str, str]] = []
    monkeypatch.setattr(controller, "_load_named_pymol_object", lambda role, name: calls.append((role, name)) or True)
    object_project = ProjectState(source_objects={"reference": "ref", "target": "target"})
    controller._commit_legacy_project(object_project, tmp_path / "objects.json")
    assert calls == [("reference", "ref"), ("target", "target")]

    monkeypatch.setattr(controller, "_load_source", lambda *_args, **_kwargs: False)
    with pytest.raises(ValueError, match="reference source"):
        controller._commit_legacy_project(legacy, tmp_path / "failed.json")
    monkeypatch.setattr(controller, "_load_named_pymol_object", lambda *_args, **_kwargs: False)
    with pytest.raises(ValueError, match="reference PyMOL"):
        controller._commit_legacy_project(object_project, tmp_path / "failed-objects.json")
    controller.close()
    panel.deleteLater()


def test_legacy_export_and_chart_failure_paths_remain_truthful(application, monkeypatch, tmp_path: Path) -> None:
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    legacy = _legacy_result()
    controller.model = controller.model.with_analysis(legacy)
    errors: list[str] = []
    monkeypatch.setattr(controller, "_show_error", errors.append)

    controller.set_v03_export_records(provenance=("stale",))
    controller.model = controller.model.with_analysis(legacy)
    monkeypatch.setattr(controller.w.QFileDialog, "getSaveFileName", lambda *_: ("", ""))
    controller._export_xlsx()
    monkeypatch.setattr(controller.w.QFileDialog, "getSaveFileName", lambda *_: (str(tmp_path / "legacy.json"), "JSON"))
    calls: list[object] = []
    controller._export("json", lambda *_args, **_kwargs: calls.append(True), "JSON")
    assert calls == [True]

    controller.model = controller.model.with_report(rich_report(tmp_path)[1])
    controller._report_request = None
    controller.set_v03_export_records(provenance=("stale",))
    controller._export_xlsx()
    assert errors and "verified source binding" in errors[-1]

    controller.model = controller.model.with_analysis(legacy)
    monkeypatch.setattr(controller, "_selected_chart_dataset", lambda _result: object())
    monkeypatch.setattr(controller.w.QFileDialog, "getSaveFileName", lambda *_: (str(tmp_path / "chart.xlsx"), "XLSX"))
    monkeypatch.setattr(
        "structlens.plugin.gui.qt_panel.export_chart_xlsx",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(OSError("disk full")),
    )
    controller._export_chart_xlsx()
    assert errors[-1] == "Could not export chart data: disk full"

    monkeypatch.setattr(controller.w.QFileDialog, "getSaveFileName", lambda *_: (str(tmp_path / "chart.png"), "PNG"))
    monkeypatch.setattr(
        "structlens.plugin.gui.qt_panel.export_chart_image",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("renderer failed")),
    )
    controller._export_chart_image("png", 300)
    assert errors[-1] == "Could not export chart image: renderer failed"
    controller.close()
    panel.deleteLater()


def test_legacy_pymol_export_and_launch_failure_boundaries(application, monkeypatch, tmp_path: Path) -> None:
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    legacy = _legacy_result()
    controller.model = controller.model.with_analysis(legacy)
    errors: list[str] = []
    monkeypatch.setattr(controller, "_show_error", errors.append)
    controller._export_pymol_bundle()
    assert "coordinate" in errors[-1]
    controller._open_in_pymol()
    assert "coordinate" in errors[-1]

    fixture = Path(__file__).parents[2] / "fixtures" / "parsing" / "numbering_altloc.pdb"
    assert controller._load_source("reference", fixture)
    assert controller._load_source("target", fixture)
    legacy = _legacy_result()
    controller.model = controller.model.with_analysis(legacy)
    output = tmp_path / "legacy.structlens-pymol"
    monkeypatch.setattr(
        controller.w.QFileDialog,
        "getSaveFileName",
        lambda *_: (str(output), "StructLens-PyMOL bundle (*.structlens-pymol)"),
    )
    monkeypatch.setattr(
        "structlens.plugin.gui.qt_panel.write_pymol_bundle",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(OSError("disk full")),
    )
    controller._export_pymol_bundle()
    assert errors[-1] == "Could not export PyMOL bundle: disk full"

    class Launcher:
        def __init__(self, *_args, **_kwargs):
            pass

        def launch_bundle(self, _path):
            return SimpleNamespace(process_id=None)

        def locate(self):
            raise OSError("not installed")

    monkeypatch.setattr("structlens.plugin.gui.qt_exports_mixin.PyMOLLauncher", Launcher)
    monkeypatch.setattr(
        "structlens.plugin.gui.qt_panel.write_pymol_bundle",
        lambda path, **_kwargs: Path(path),
    )
    controller._open_in_pymol()
    assert "Validated bundle prepared" in controller.footer_status.text()
    controller.close()
    panel.deleteLater()


def test_source_mixin_snapshot_materialization_and_empty_collection_paths(
    application, monkeypatch, tmp_path: Path
) -> None:
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    fixture = Path(__file__).parents[2] / "fixtures" / "parsing" / "numbering_altloc.pdb"

    class PyMol:
        def __init__(self) -> None:
            self.names: set[str] = set()
            self.loaded: list[tuple[str, str]] = []

        def get_names(self, _kind: str) -> list[str]:
            return sorted(self.names)

        def load(self, path: str, name: str) -> None:
            self.loaded.append((path, name))
            self.names.add(name)

    command = PyMol()
    controller.command = command
    assert controller._load_source("reference", fixture)
    assert command.loaded
    loaded_count = len(command.loaded)
    # The second role has a distinct namespace.
    assert controller._load_source("target", fixture)
    assert len(command.loaded) == loaded_count + 1
    target_name = command.loaded[-1][1]
    controller._load_file_into_pymol("target", fixture)
    assert len(command.loaded) == loaded_count + 2
    assert target_name in command.names
    assert command.loaded[-1][1] != target_name
    assert command.loaded[-1][1].startswith("structlens_panel_numbering_altloc_")

    controller.reference_edit.clear()
    controller.target_edit.clear()
    assert controller._load_sources_from_edits() is False
    empty_structure = SimpleNamespace(chains=(), structure_id="empty")
    controller._populate_models(controller.reference_model_combo, empty_structure)
    controller._populate_chains(controller.reference_chain_combo, empty_structure)
    controller._model_changed("reference")
    controller.close()
    panel.deleteLater()


def test_source_mixin_compare_guards_and_compatibility_boundaries(application, monkeypatch, tmp_path: Path) -> None:
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    fixture = Path(__file__).parents[2] / "fixtures" / "parsing" / "numbering_altloc.pdb"
    errors: list[str] = []
    monkeypatch.setattr(controller, "_show_error", errors.append)

    controller._start_analysis()
    assert errors
    assert controller.nav.currentRow() == 0
    assert controller._load_source("reference", fixture)
    assert controller._load_source("target", fixture)

    controller.reference_edit.setText(str(tmp_path / "edited.pdb"))
    controller._start_analysis()
    assert "edited" in errors[-1]
    controller.reference_edit.setText(str(controller.reference_loaded_source.path))
    controller.target_edit.setText(str(controller.target_loaded_source.path))
    assert controller.reference_edit.text() == str(controller.reference_loaded_source.path)
    assert controller.target_edit.text() == str(controller.target_loaded_source.path)
    selected_chain = controller._selected_chain
    monkeypatch.setattr(controller, "_selected_chain", lambda *_args: None)
    controller._start_analysis()
    assert "chain" in errors[-1].lower()
    monkeypatch.setattr(controller, "_selected_chain", selected_chain)
    assert controller._load_source("reference", fixture)
    assert controller._load_source("target", fixture)

    verify_current_source = controller._report_controller.verify_current_source
    monkeypatch.setattr(
        controller._report_controller,
        "verify_current_source",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(OSError("capture race")),
    )
    controller._start_analysis()
    assert "capture" in errors[-1].lower()
    monkeypatch.setattr(controller._report_controller, "verify_current_source", verify_current_source)
    assert controller._load_source("reference", fixture)
    assert controller._load_source("target", fixture)

    reference = controller.reference_structure.chains[0]
    target = controller.target_structure.chains[0]
    controller.manual_edit.clear()
    with pytest.raises(ValueError, match="at least one"):
        controller._manual_pairs(reference, target)
    controller._report_controller.__dict__["build_request"] = lambda *_args, **_kwargs: SimpleNamespace()
    built = controller._build_snapshot_request(
        controller.reference_loaded_source,
        controller.target_loaded_source,
        reference=reference,
        target=target,
        analysis_settings=controller._settings(),
        site_definitions=(),
        manual_pairs=(),
    )
    assert isinstance(built, SimpleNamespace)
    del controller._report_controller.__dict__["build_request"]

    _, report = rich_report(tmp_path)
    controller.show_report(report)
    controller.show_presentation(SimpleNamespace())
    controller.show_status("compat status")
    assert controller.footer_status.text() == "compat status"
    controller.close()
    panel.deleteLater()


def test_export_callbacks_and_canonical_compatibility_failures(application, monkeypatch, tmp_path: Path) -> None:
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    errors: list[str] = []
    monkeypatch.setattr(controller, "_show_error", errors.append)

    # A staged legacy payload without a legacy result cannot be exported.
    controller.set_v03_export_records(provenance=("stale",))
    controller._export_xlsx()
    assert errors[-1] == "Run a comparison before exporting results."

    legacy = _legacy_result()
    controller.model = controller.model.with_analysis(legacy)
    output = tmp_path / "legacy.xlsx"
    monkeypatch.setattr(controller.w.QFileDialog, "getSaveFileName", lambda *_: (str(output), "XLSX"))
    monkeypatch.setattr(
        "structlens.plugin.gui.qt_panel.export_v03_xlsx",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(OSError("disk full")),
    )
    controller._export_xlsx()
    assert errors[-1] == "Could not export v0.3 XLSX: disk full"

    # The generic XLSX/TSV/JSON routes remain callable for explicit legacy
    # integrations and receive exactly the legacy result.
    controller._v03_export_records = {}
    monkeypatch.setattr("structlens.plugin.gui.qt_exports_mixin.export_analysis_xlsx", lambda *_args, **_kwargs: None)
    controller._export_xlsx()
    captured: list[tuple[object, object]] = []
    monkeypatch.setattr(
        controller.w.QFileDialog,
        "getSaveFileName",
        lambda _widget, title, *_args: (str(tmp_path / f"{title}.out"), title),
    )
    monkeypatch.setattr(
        "structlens.plugin.gui.qt_panel.export_analysis_tsv",
        lambda value, path, **_kwargs: captured.append((value, path)),
    )
    monkeypatch.setattr(
        "structlens.plugin.gui.qt_exports_mixin.export_analysis_json",
        lambda value, path, **_kwargs: captured.append((value, path)),
    )
    controller._export_tsv()
    controller._export_json()
    assert len(captured) >= 2 and all(value is legacy for value, _path in captured)

    # A canonical report with stale request identity is refused before any
    # legacy PyMOL bundle path can be reached.
    canonical_root = tmp_path / "canonical"
    canonical_root.mkdir()
    request, report = rich_report(canonical_root)
    assert controller._load_source("reference", Path(request.reference_selection.path))
    assert controller._load_source("target", Path(request.target_selection.path))
    controller._report_request = request
    controller._analysis_finished(report)
    controller._report_request = replace(request, target_selection=replace(request.target_selection, model_id="wrong"))
    controller._export_pymol_bundle()
    assert errors and "stale" in errors[-1].lower()
    controller.close()
    panel.deleteLater()


def test_export_chart_and_pymol_launch_lifecycle_edges(application, monkeypatch, tmp_path: Path) -> None:
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    errors: list[str] = []
    monkeypatch.setattr(controller, "_show_error", errors.append)
    controller._open_in_pymol()
    assert errors[-1] == "Run a comparison before opening PyMOL."

    legacy = _legacy_result()
    controller.model = controller.model.with_analysis(legacy)
    controller._open_in_pymol()
    assert "coordinate" in errors[-1]
    controller._export_pymol_bundle()
    assert "coordinate" in errors[-1]

    fixture = Path(__file__).parents[2] / "fixtures" / "parsing" / "numbering_altloc.pdb"
    assert controller._load_source("reference", fixture)
    assert controller._load_source("target", fixture)
    controller.model = controller.model.with_analysis(legacy)
    monkeypatch.setattr(controller.w.QFileDialog, "getSaveFileName", lambda *_: ("", ""))
    before_cancel = errors[-1]
    controller._export_pymol_bundle()
    assert errors[-1] == before_cancel

    class Launcher:
        def __init__(self, *_args, **_kwargs):
            pass

        def launch_bundle(self, _path):
            return SimpleNamespace(process_id=1)

        def locate(self):
            return Path("pymol.exe")

    monkeypatch.setattr("structlens.plugin.gui.qt_exports_mixin.PyMOLLauncher", Launcher)
    monkeypatch.setattr(
        "structlens.plugin.gui.qt_panel.write_pymol_bundle",
        lambda path, **_kwargs: Path(path),
    )
    controller._open_in_pymol()
    assert "PyMOL launched" in controller.footer_status.text()

    monkeypatch.setattr(
        "structlens.plugin.gui.qt_panel.write_pymol_bundle",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(OSError("writer failed")),
    )
    controller._open_in_pymol()
    assert errors[-1] == "Could not open in PyMOL: writer failed"

    controller.model = StructLensPanelModel()
    controller._export_chart_image("png", 300)
    assert "comparison" in errors[-1].lower()
    controller.model = controller.model.with_analysis(legacy)
    assert controller.model.analysis is legacy
    controller.chart_combo.setCurrentText("Pairwise similarity heatmap")
    controller._export_chart_image("png", 300)
    assert "dataset" in errors[-1].lower()
    monkeypatch.setattr(controller, "_selected_chart_dataset", lambda _result: object())
    monkeypatch.setattr(controller.w.QFileDialog, "getSaveFileName", lambda *_: ("", ""))
    controller._export_chart_image("png", 300)
    monkeypatch.setattr(controller.w.QFileDialog, "getSaveFileName", lambda *_: (str(tmp_path / "chart.png"), "PNG"))
    monkeypatch.setattr("structlens.plugin.gui.qt_panel.export_chart_image", lambda *_args, **_kwargs: None)
    controller._export_chart_image("png", 300)
    assert "Chart image exported" in controller.footer_status.text()
    controller.close()
    panel.deleteLater()
