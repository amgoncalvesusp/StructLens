"""Task 14 RED regressions for canonical GUI and PyMOL boundaries.

These tests are intentionally written against the typed contracts used by the
application.  They protect the GUI from presenting or exporting data that was
produced with a different request, and protect a live PyMOL session from
name/temporary-file side effects.
"""

from __future__ import annotations

import os
from dataclasses import replace
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from structlens.application.dto import AnalysisReportRequest
from structlens.application.report_service import ReportService
from structlens.core.evidence import InteractionEvidence
from structlens.core.interactions import InteractionRecord, InteractionType
from structlens.core.models import AlignmentMode, AnalysisResult, AnalysisSettings, ResidueCorrespondence, ResidueId
from structlens.core.parsing import InputSelection, StructureFormat, capture_snapshot
from structlens.core.sites import SiteDefinition, SiteDefinitionMode
from structlens.integrations.pymol.adapter import PyMOLAdapter
from structlens.plugin.gui.model import CanonicalReportBinding
from structlens.plugin.gui.presentation import present_report

_FIXTURE = Path("tests/fixtures/parsing/enriched.pdb")


def _request(tmp_path: Path, *, manual: bool = False) -> AnalysisReportRequest:
    source = _FIXTURE.read_bytes()
    reference_path = tmp_path / "reference.pdb"
    target_path = tmp_path / "target.pdb"
    reference_path.write_bytes(source)
    target_path.write_bytes(source)
    reference = capture_snapshot(reference_path)
    target = capture_snapshot(target_path)
    reference_selection = InputSelection(
        reference.content_id,
        reference.display_name,
        StructureFormat.PDB,
        "1",
        author_chain_ids=("A",),
        path=str(reference_path),
    )
    target_selection = InputSelection(
        target.content_id,
        target.display_name,
        StructureFormat.PDB,
        "1",
        author_chain_ids=("A",),
        path=str(target_path),
    )
    reference_id = ResidueId("reference", "1", "A", "100", None, "ALA")
    target_id = ResidueId("target", "1", "A", "100", None, "ALA")
    analysis_settings = AnalysisSettings(
        alignment_mode=AlignmentMode.MANUAL if manual else AlignmentMode.SEQUENCE,
    )
    return AnalysisReportRequest(
        reference,
        target,
        reference_selection,
        target_selection,
        analysis_settings=analysis_settings,
        manual_pairs=((reference_id, target_id),) if manual else (),
        site_definitions=(
            SiteDefinition(
                "active",
                "Active site",
                SiteDefinitionMode.KEY_RESIDUES,
                (reference_id,),
            ),
        ),
    )


@pytest.mark.parametrize(
    ("name", "change"),
    (
        (
            "alignment settings",
            lambda request: replace(
                request,
                analysis_settings=replace(request.analysis_settings, refinement_cutoff_angstrom=3.0),
            ),
        ),
        (
            "MSA settings",
            lambda request: replace(
                request,
                msa_settings=replace(request.msa_settings, algorithm="fallback-test"),
            ),
        ),
        (
            "interaction distance threshold",
            lambda request: replace(
                request,
                interaction_thresholds=replace(request.interaction_thresholds, hbond_distance_angstrom=3.6),
            ),
        ),
        (
            "interaction angle threshold",
            lambda request: replace(
                request,
                interaction_thresholds=replace(
                    request.interaction_thresholds,
                    pi_parallel_angle_tolerance_degrees=31.0,
                ),
            ),
        ),
        (
            "site definition",
            lambda request: replace(
                request,
                site_definitions=(replace(request.site_definitions[0], name="Different site"),),
            ),
        ),
        (
            "minimum vector magnitude",
            lambda request: replace(
                request,
                minimum_vector_magnitude_angstrom=1.0,
            ),
        ),
        ("maximum vector count", lambda request: replace(request, maximum_vectors=101)),
    ),
)
def test_canonical_binding_rejects_every_normalized_request_mismatch(tmp_path: Path, name: str, change: object) -> None:
    """A report must never be accepted under a request with different provenance."""

    request = _request(tmp_path)
    report = ReportService().analyze(request)
    altered = change(request)  # type: ignore[operator]

    with pytest.raises(ValueError, match="request|provenance|parameter|match|coherent"):
        CanonicalReportBinding.create(report, altered)


def test_canonical_binding_accepts_exact_request_and_is_report_id_stable(tmp_path: Path) -> None:
    request = _request(tmp_path)
    report = ReportService().analyze(request)

    binding = CanonicalReportBinding.create(report, request)

    assert binding.report.report_id == report.report_id
    assert binding.request == request
    assert binding.presentation == present_report(report)


def test_canonical_binding_rejects_manual_pair_provenance_mismatch(tmp_path: Path) -> None:
    request = _request(tmp_path, manual=True)
    report = ReportService().analyze(request)
    altered = replace(
        request,
        manual_pairs=(
            (
                ResidueId("reference", "1", "A", "100A", None, "GLY"),
                ResidueId("target", "1", "A", "100", None, "ALA"),
            ),
        ),
    )

    with pytest.raises(ValueError, match="request|provenance|parameter|match|coherent"):
        CanonicalReportBinding.create(report, altered)


class _CollisionCommand:
    def __init__(self, existing: set[str]) -> None:
        self.existing = set(existing)
        self.calls: list[tuple[str, tuple[object, ...]]] = []

    def get_names(self, kind: str) -> list[str]:
        del kind
        return sorted(self.existing)

    def select(self, name: str, expression: str) -> None:
        self.calls.append(("select", (name, expression)))

    def zoom(self, name: str) -> None:
        self.calls.append(("zoom", (name,)))

    def create(self, name: str, expression: str, *states: int) -> None:
        self.calls.append(("create", (name, expression, *states)))

    def show(self, *args: object) -> None:
        self.calls.append(("show", args))

    def color(self, *args: object) -> None:
        self.calls.append(("color", args))

    def label(self, *args: object) -> None:
        self.calls.append(("label", args))

    def delete(self, name: str) -> None:
        self.calls.append(("delete", (name,)))


def _correspondence() -> ResidueCorrespondence:
    return ResidueCorrespondence(
        0,
        ResidueId("reference", "1", "A", "100", None, "ALA"),
        ResidueId("target", "1", "A", "100", None, "ALA"),
        "A",
        "A",
        "conserved",
    )


def test_pymol_adapter_uses_collision_safe_owned_names_and_never_deletes_user_names() -> None:
    existing = {
        "structlens_default_target_focus",
        "structlens_default_target_target_view",
        "structlens_default_reference_reference_view",
    }
    command = _CollisionCommand(existing)
    adapter = PyMOLAdapter(command, project_id="default")

    focused = adapter.focus_residue(_correspondence(), "target")
    adapter.apply(
        AnalysisResult(
            "reference",
            "target",
            (_correspondence(),),
            (),
            1.0,
            1.0,
            "accepted",
        )
    )
    adapter.reset()

    created_or_selected = [args[0] for method, args in command.calls if method in {"create", "select"}]
    assert focused not in existing
    assert len(created_or_selected) == len(set(created_or_selected))
    deleted = {args[0] for method, args in command.calls if method == "delete"}
    assert deleted.isdisjoint(existing)


def test_presentation_retains_raw_interactions_when_no_difference_exists(tmp_path: Path) -> None:
    request = _request(tmp_path)
    report = ReportService().analyze(request)
    residue_a = ResidueId("reference", "1", "A", "100", None, "ALA")
    record = InteractionRecord(
        "reference",
        InteractionType.HYDROPHOBIC,
        residue_a,
        None,
        "CA",
        None,
        3.2,
        angle_degrees=45.0,
    )
    evidence = InteractionEvidence((), (record,), (replace(record, structure_id="target"),))
    report = replace(report, interactions=evidence)

    presentation = present_report(report)
    rows = presentation.sections.interactions.rows

    assert len(rows) == 2
    assert {row["reference_distance_angstrom"] for row in rows} == {"3.2", "Unavailable — value not reported"}
    assert {row["target_distance_angstrom"] for row in rows} == {"3.2", "Unavailable — value not reported"}
    assert {row["reference_angle_degrees"] for row in rows} == {"45.0", "Unavailable — value not reported"}


def test_gui_snapshot_capture_removes_temporary_coordinate_file_after_success_and_failure() -> None:
    class Command(_CollisionCommand):
        def __init__(self, fail: bool) -> None:
            super().__init__(set())
            self.fail = fail

        def load(self, path: str, object_name: str) -> None:
            self.calls.append(("load", (path, object_name)))
            if self.fail:
                raise RuntimeError("simulated PyMOL load failure")

    from PySide6 import QtWidgets

    from structlens.plugin.gui.main_panel import build_qt_panel

    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    fixture = _FIXTURE
    for fail in (False, True):
        panel = build_qt_panel(command=Command(fail))
        controller = panel._structlens_controller
        assert controller._load_source("reference", fixture) is True
        assert all(not path.exists() for path in controller._temporary_paths)
        controller.close()
        panel.deleteLater()
    app.processEvents()


def test_canonical_tabular_export_rejects_stale_binding_before_save_dialog(tmp_path: Path, monkeypatch) -> None:
    from PySide6 import QtWidgets

    from structlens.plugin.gui.main_panel import build_qt_panel

    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    request = _request(tmp_path)
    report = ReportService().analyze(request)
    controller._populate_report(report, request=request)
    controller._report_request = replace(
        request,
        minimum_vector_magnitude_angstrom=request.minimum_vector_magnitude_angstrom + 0.25,
    )
    dialogs: list[object] = []
    monkeypatch.setattr(
        controller.w.QFileDialog,
        "getSaveFileName",
        lambda *args: dialogs.append(args) or (str(tmp_path / "should-not-open.csv"), "CSV (*.csv)"),
    )
    errors: list[str] = []
    controller._show_error = errors.append

    controller._export("csv", lambda *_args, **_kwargs: None, "CSV")

    assert dialogs == []
    assert errors and "canonical" in errors[0].casefold()
    controller.close()
    panel.deleteLater()
    app.processEvents()
