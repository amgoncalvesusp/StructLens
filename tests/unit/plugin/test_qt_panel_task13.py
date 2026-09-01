"""Task 13 GUI report-first integration and source-policy regressions.

The site/request boundaries and the real fixture workflow are kept separate
from the historical Qt panel smoke tests.  PySide6 is imported directly so a
missing supported GUI dependency is a visible release failure, never a skip.
"""

import ast
import inspect
import os
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from PySide6 import QtWidgets  # noqa: E402

qt_widgets = QtWidgets

from structlens.core.models import (  # noqa: E402
    AnalysisResult,
    AtomRecord,
    ProteinChain,
    ProteinStructure,
    ResidueId,
    ResidueNumbering,
    ResidueRecord,
)
from structlens.core.sites import SiteDefinition, SiteDefinitionMode, SiteMetrics  # noqa: E402
from structlens.plugin.gui import qt_panel  # noqa: E402
from structlens.plugin.gui.main_panel import build_qt_panel  # noqa: E402


@pytest.fixture(scope="module")
def application():
    app = qt_widgets.QApplication.instance() or qt_widgets.QApplication([])
    yield app
    app.quit()


def test_site_control_records_definition_for_next_report_without_service_call(application, monkeypatch) -> None:
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    fixture = Path("tests/fixtures/parsing/numbering_altloc.pdb")
    assert controller._load_source("reference", fixture) is True
    assert controller._load_source("target", fixture) is True
    reference = controller._selected_chain(controller.reference_structure, controller.reference_chain_combo)
    assert reference is not None
    residue = reference.residue_records[0].residue_id
    controller.site_mode_combo.setCurrentText("Key residues")
    controller.site_residues_edit.setText(f"{residue.chain_id}:{residue.auth_seq_id}{residue.insertion_code or ''}")
    # A legacy analysis is deliberately present so the old implementation
    # reaches its service call.  The report-first action must only retain the
    # definition as request input and must never calculate metrics here.
    controller.model = controller.model.with_analysis(
        AnalysisResult(
            reference_id="reference",
            target_id="target",
            correspondences=(),
            mutations=(),
            sequence_identity=1.0,
            sequence_coverage=1.0,
            alignment_decision="compatibility",
        )
    )

    class ForbiddenSiteService:
        def calculate_site_metrics(self, *_args: object, **_kwargs: object) -> object:
            raise AssertionError("site metrics must be calculated by the report service, not the Qt action")

    monkeypatch.setattr(qt_panel, "site_service", ForbiddenSiteService(), raising=False)
    errors: list[str] = []
    monkeypatch.setattr(controller, "_show_error", errors.append)
    controller._define_site_from_controls()

    captured: list[dict[str, object]] = []
    sentinel_request = object()
    monkeypatch.setattr(
        controller._report_controller,
        "build_request",
        lambda *_args, **kwargs: captured.append(kwargs) or sentinel_request,
    )
    monkeypatch.setattr(controller._report_controller, "execute", lambda _request: None)
    controller._start_analysis()

    assert not errors
    assert len(captured) == 1
    definitions = captured[0]["site_definitions"]
    assert len(definitions) == 1
    definition = definitions[0]
    assert isinstance(definition, SiteDefinition)
    assert definition.mode is SiteDefinitionMode.KEY_RESIDUES
    assert definition.reference_residues == (residue,)
    assert controller._site_metrics == ()

    controller.close()
    panel.deleteLater()


def test_legacy_synthetic_site_control_stays_display_only_and_skips_site_service(application, monkeypatch) -> None:
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller

    reference_id = ResidueId("reference", "1", "A", "1", None, "ALA")
    target_id = ResidueId("target", "1", "A", "1", None, "VAL")
    reference_record = ResidueRecord(
        reference_id,
        ResidueNumbering("1", "1", None),
        "ALA",
        "A",
        (AtomRecord("CA", "C", (0.0, 0.0, 0.0)),),
    )
    target_record = ResidueRecord(
        target_id,
        ResidueNumbering("1", "1", None),
        "VAL",
        "V",
        (AtomRecord("CA", "C", (1.0, 0.0, 0.0)),),
    )
    controller.reference_structure = ProteinStructure(
        "reference", (ProteinChain("reference", "1", "A", (reference_id,), "A", (reference_record,)),)
    )
    controller.target_structure = ProteinStructure(
        "target", (ProteinChain("target", "1", "A", (target_id,), "V", (target_record,)),)
    )
    controller._populate_models(controller.reference_model_combo, controller.reference_structure)
    controller._populate_models(controller.target_model_combo, controller.target_structure)
    controller._populate_chains(
        controller.reference_chain_combo,
        controller.reference_structure,
        model_id="1",
    )
    controller._populate_chains(controller.target_chain_combo, controller.target_structure, model_id="1")
    controller.site_mode_combo.setCurrentText("Key residues")
    controller.site_residues_edit.setText("A:1")
    controller.model = controller.model.with_analysis(
        AnalysisResult(
            reference_id="reference",
            target_id="target",
            correspondences=(),
            mutations=(),
            sequence_identity=1.0,
            sequence_coverage=1.0,
            alignment_decision="compatibility",
        )
    )

    class ForbiddenSiteService:
        def calculate_site_metrics(self, *_args: object, **_kwargs: object) -> object:
            raise AssertionError("legacy synthetic site control must not calculate metrics")

    backend_calls: list[str] = []
    monkeypatch.setattr(
        "structlens.plugin.gui.qt_legacy.legacy_site_backend",
        lambda: backend_calls.append("called") or ForbiddenSiteService(),
    )
    monkeypatch.setattr(controller, "_show_error", lambda _message: None)
    controller._define_site_from_controls()

    captured: list[dict[str, object]] = []
    sentinel_request = object()
    monkeypatch.setattr(controller, "_load_sources_from_edits", lambda: True)
    monkeypatch.setattr(
        controller._report_controller,
        "build_request",
        lambda *_args, **kwargs: captured.append(kwargs) or sentinel_request,
    )
    monkeypatch.setattr(controller._report_controller, "execute", lambda _request: None)
    controller._start_analysis()

    assert backend_calls == []
    assert len(captured) == 1
    definitions = captured[0]["site_definitions"]
    assert len(definitions) == 1
    definition = definitions[0]
    assert isinstance(definition, SiteDefinition)
    assert definition.reference_residues == (reference_id,)
    with pytest.raises((AttributeError, TypeError)):
        definition.name = "mutated"

    controller.set_site_metrics((SiteMetrics("display-only", "target", 1, 1.0),))
    assert controller.site_metrics_table.rowCount() == 1
    controller.close()
    panel.deleteLater()


def test_qt_source_has_no_dynamic_site_service_or_private_controller_mutation() -> None:
    source = inspect.getsource(qt_panel)
    tree = ast.parse(source)

    assert "site_service" not in source
    assert not any(
        isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr in {"import_module", "calculate_site_metrics"}
        for node in ast.walk(tree)
    )

    def touches_report_controller(node: ast.AST) -> bool:
        if isinstance(node, ast.Attribute):
            return node.attr == "_report_controller" or touches_report_controller(node.value)
        return False

    assignments: list[ast.AST] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            assignments.extend(node.targets)
        elif isinstance(node, ast.AnnAssign):
            assignments.append(node.target)
        elif isinstance(node, ast.AugAssign):
            assignments.append(node.target)
    assert not any(touches_report_controller(target) for target in assignments)


def test_real_altloc_gui_workflow_keeps_implicit_model_zero_and_export_identity(
    application, monkeypatch, tmp_path: Path
) -> None:
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    fixture = Path("tests/fixtures/parsing/numbering_altloc.pdb")
    assert controller._load_source("reference", fixture) is True
    assert controller._load_source("target", fixture) is True
    reference_chain = controller._selected_chain(controller.reference_structure, controller.reference_chain_combo)
    target_chain = controller._selected_chain(controller.target_structure, controller.target_chain_combo)
    assert reference_chain is not None
    assert target_chain is not None
    assert reference_chain.model_id == "0"
    assert target_chain.model_id == "0"

    controller._start_analysis()
    future = controller._future
    assert future is not None
    report = future.result(timeout=20)
    controller._poll_analysis()

    request = controller._report_request
    assert request is not None
    assert request.reference_selection.model_id == reference_chain.model_id == "0"
    assert request.target_selection.model_id == target_chain.model_id == "0"
    assert controller.model.report is report

    output = tmp_path / "altloc-report.json"
    monkeypatch.setattr(
        controller.w.QFileDialog,
        "getSaveFileName",
        lambda *_args: (str(output), "JSON (*.json)"),
    )
    calls: list[tuple[object, object, object]] = []

    def exporter(value: object, path: object, **kwargs: object) -> None:
        calls.append((value, path, kwargs.get("snapshots")))

    controller._export("json", exporter, "JSON")

    assert calls == [(report, str(output), (request.reference_snapshot, request.target_snapshot))]
    controller.close()
    panel.deleteLater()
