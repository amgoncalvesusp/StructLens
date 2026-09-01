from dataclasses import dataclass, replace
from pathlib import Path

import pytest

from structlens.core.evidence import Availability
from structlens.core.models import (
    AlignmentMode,
    AnalysisResult,
    AnalysisSelection,
    AnalysisSettings,
    CorrespondenceStatus,
    ResidueCorrespondence,
    ResidueId,
)
from structlens.core.parsing import InputSelection, StructureFormat
from structlens.core.quality import StructureQualityReport
from structlens.core.reports import AnalysisReport, AnalysisSnapshot, InputQualityBundle, SectionAvailability
from structlens.core.sites import SiteDefinition, SiteDefinitionMode
from structlens.plugin.gui.main_panel import (
    GUI_SECTIONS,
    WORKFLOW_HELP,
    AnalysisReportController,
    StructLensPanelModel,
)


def _minimal_report(*, analysis: AnalysisResult | None = None) -> AnalysisReport:
    selection = InputSelection("a" * 64, "reference.pdb", StructureFormat.PDB, "1")
    target = InputSelection("b" * 64, "target.pdb", StructureFormat.PDB, "1")
    return AnalysisReport(
        selection,
        target,
        InputQualityBundle(
            StructureQualityReport(Availability.AVAILABLE),
            StructureQualityReport(Availability.AVAILABLE),
        ),
        analysis=AnalysisSnapshot.from_result(analysis) if analysis is not None else None,
        availability=SectionAvailability(
            input_quality=Availability.AVAILABLE,
            analysis=Availability.AVAILABLE if analysis is not None else Availability.NOT_APPLICABLE,
        ),
    )


def _explicit_model_one(path: Path) -> Path:
    fixture = Path("tests/fixtures/parsing/numbering_altloc.pdb")
    lines = [line for line in fixture.read_text(encoding="utf-8").splitlines() if line != "END"]
    path.write_text("MODEL        1\n" + "\n".join(lines) + "\nENDMDL\nEND\n", encoding="utf-8")
    return path


def test_gui_defines_six_scientific_sections_and_context_help() -> None:
    assert GUI_SECTIONS == (
        "Project",
        "Alignment",
        "Mutations",
        "Residues",
        "Visualization",
        "Results",
    )
    for key in (
        "Auto",
        "Sequence",
        "Structure",
        "Manual",
        "strict RMSD",
        "Cα displacement",
    ):
        assert key in WORKFLOW_HELP
        assert "what" in WORKFLOW_HELP[key].lower()


def test_panel_model_keeps_ui_state_separate_from_analysis_state() -> None:
    model = StructLensPanelModel()
    updated = model.with_status("Ready")
    assert updated.status == "Ready"
    assert updated.analysis is None
    assert model.status == "Choose a reference structure and a target to begin."


def test_panel_status_transition_leaves_cancelled_model_idle() -> None:
    model = StructLensPanelModel().with_busy("Comparing structures…")
    cancelled = model.with_status("Comparison cancelled.")

    assert model.busy is True
    assert cancelled.busy is False
    assert cancelled.status == "Comparison cancelled."


def test_report_controller_delivers_one_canonical_artifact_to_fake_widget() -> None:
    selection = InputSelection("a" * 64, "reference.pdb", StructureFormat.PDB, "1")
    target = InputSelection("b" * 64, "target.pdb", StructureFormat.PDB, "1")
    report = AnalysisReport(
        selection,
        target,
        InputQualityBundle(
            StructureQualityReport(Availability.AVAILABLE),
            StructureQualityReport(Availability.AVAILABLE),
        ),
        availability=SectionAvailability(input_quality=Availability.AVAILABLE),
    )

    class Service:
        def analyze(self, request: object) -> AnalysisReport:
            assert request is sentinel
            return report

    @dataclass
    class FakeWidget:
        rendered: AnalysisReport | None = None
        status: str = ""

        def show_report(self, value: AnalysisReport) -> None:
            self.rendered = value

        def show_status(self, value: str) -> None:
            self.status = value

    sentinel = object()
    widget = FakeWidget()
    controller = AnalysisReportController(Service(), widget)

    result = controller.analyze(sentinel)

    assert result is report
    assert controller.model.report is report
    assert widget.rendered is report
    assert report.report_id in widget.status
    assert not any("payload" in name for name in StructLensPanelModel.__dataclass_fields__)


def test_report_controller_delivers_one_exact_report_and_its_presentation() -> None:
    selection = InputSelection("a" * 64, "reference.pdb", StructureFormat.PDB, "1")
    target = InputSelection("b" * 64, "target.pdb", StructureFormat.PDB, "1")
    report = AnalysisReport(
        selection,
        target,
        InputQualityBundle(
            StructureQualityReport(Availability.AVAILABLE),
            StructureQualityReport(Availability.AVAILABLE),
        ),
        availability=SectionAvailability(input_quality=Availability.AVAILABLE),
    )

    class Service:
        def analyze(self, request: object) -> AnalysisReport:
            assert request is sentinel
            return report

    @dataclass
    class FakeWidget:
        rendered: AnalysisReport | None = None
        presented: object | None = None
        status: str = ""

        def show_report(self, value: AnalysisReport) -> None:
            self.rendered = value

        def show_presentation(self, value: object) -> None:
            self.presented = value

        def show_status(self, value: str) -> None:
            self.status = value

    sentinel = object()
    widget = FakeWidget()
    controller = AnalysisReportController(Service(), widget)

    result = controller.analyze(sentinel)

    assert result is report
    assert widget.rendered is report
    assert widget.presented is not None
    assert widget.presented.report_id == report.report_id


def test_report_controller_build_request_keeps_pdb_model_chain_sites_and_manual_pairs(tmp_path: Path) -> None:
    controller = AnalysisReportController(service=object())
    implicit = Path("tests/fixtures/parsing/numbering_altloc.pdb")
    explicit = _explicit_model_one(tmp_path / "numbering_altloc_model_1.pdb")
    reference_residue = ResidueId("numbering_altloc", "0", "A", "100", None, "ALA")
    target_residue = ResidueId("numbering_altloc", "0", "A", "100", None, "ALA")
    site = SiteDefinition(
        "active",
        mode=SiteDefinitionMode.KEY_RESIDUES,
        reference_residues=(reference_residue,),
    )
    settings = AnalysisSettings(alignment_mode=AlignmentMode.MANUAL)
    manual_pairs = ((reference_residue, target_residue),)

    request = controller.build_request(
        implicit,
        implicit,
        reference_chain_id="A",
        target_chain_id="A",
        analysis_settings=settings,
        site_definitions=(site,),
        manual_pairs=manual_pairs,
    )

    # Biopython represents a PDB with no MODEL record as model 0.  The GUI
    # must not manufacture model 1 merely because that is the CLI spelling of
    # its default model selector.
    assert request.reference_selection.model_id == "0"
    assert request.target_selection.model_id == "0"
    assert request.reference_selection.author_chain_ids == ("A",)
    assert request.target_selection.author_chain_ids == ("A",)
    assert request.site_definitions == (site,)
    assert request.manual_pairs == manual_pairs

    explicit_request = controller.build_request(explicit, explicit, reference_chain_id="A", target_chain_id="A")
    assert explicit_request.reference_selection.model_id == "1"
    assert explicit_request.target_selection.model_id == "1"
    assert explicit_request.reference_selection.author_chain_ids == ("A",)


def test_report_controller_validates_service_outputs_and_optional_view_hooks() -> None:
    report = _minimal_report()

    class ValidService:
        def analyze(self, request: object) -> AnalysisReport:
            assert request is sentinel
            return report

    class InvalidService:
        def analyze(self, _request: object) -> object:
            return object()

    @dataclass
    class MinimalView:
        rendered: AnalysisReport | None = None
        status: str = ""

        def show_report(self, value: AnalysisReport) -> None:
            self.rendered = value

        def show_status(self, value: str) -> None:
            self.status = value

    sentinel = object()
    with pytest.raises(TypeError, match="AnalysisReport"):
        AnalysisReportController(InvalidService()).execute(sentinel)
    with pytest.raises(TypeError, match="AnalysisReport"):
        AnalysisReportController().publish(object())

    headless = AnalysisReportController(ValidService())
    assert headless.analyze(sentinel) is report
    assert headless.last_report is report
    assert headless.last_presentation is not None

    view = MinimalView()
    controller = AnalysisReportController(ValidService(), view)
    assert controller.analyze(sentinel) is report
    assert view.rendered is report
    assert report.report_id in view.status


def test_with_report_makes_report_primary_and_clears_legacy_analysis_and_selection(tmp_path) -> None:
    from structlens.application.report_service import ReportService
    from tests.unit.application.test_report_service import _request

    # The model still accepts legacy AnalysisResult integrations, but a
    # canonical report must replace that competing state and any row focus.
    legacy = AnalysisResult(
        reference_id="legacy-reference",
        target_id="legacy-target",
        correspondences=(),
        mutations=(),
        sequence_identity=1.0,
        sequence_coverage=1.0,
        alignment_decision="legacy",
    )
    request = _request(tmp_path)
    report = ReportService().analyze(request)
    model = replace(
        StructLensPanelModel().with_analysis(legacy),
        selection=AnalysisSelection(reference_residue=ResidueId("r", "1", "A", "1", None, "ALA")),
    )

    updated = model.with_report(report)

    assert updated.report is report
    assert updated.analysis is None
    assert updated.selection is None
    assert model.analysis is legacy


def test_selection_controller_transitions_are_immutable_and_report_bound(tmp_path) -> None:
    from structlens.application.report_service import ReportService
    from structlens.plugin.gui.selection_controller import AnalysisSelectionController
    from tests.unit.application.test_report_service import _request

    report = ReportService().analyze(_request(tmp_path))
    assert report.analysis is not None
    assert report.analysis.correspondences
    assert report.evidence_cards

    controller = AnalysisSelectionController()
    bound = controller.bind_report(report)
    selected = bound.select_alignment_index(report.analysis.correspondences[0].alignment_index)

    assert bound is not selected
    assert selected.report_id == report.report_id
    assert selected.selection.reference_residue == report.analysis.correspondences[0].reference
    assert selected.selection.target_residue == report.analysis.correspondences[0].target
    assert selected.lookup_correspondence() == report.analysis.correspondences[0]
    assert selected.lookup_evidence_card() in report.evidence_cards
    assert controller.report_id is None
    assert controller.selection == AnalysisSelection()

    cleared = selected.clear_selection()
    assert cleared is not selected
    assert cleared.report_id == report.report_id
    assert cleared.selection == AnalysisSelection()

    replacement = cleared.bind_report(report)
    assert replacement is not cleared
    assert replacement.selection == AnalysisSelection()


def test_selection_controller_rejects_invalid_or_foreign_report_focus(tmp_path) -> None:
    from structlens.application.report_service import ReportService
    from structlens.plugin.gui.selection_controller import AnalysisSelectionController
    from tests.unit.application.test_report_service import _request

    report = ReportService().analyze(_request(tmp_path))
    controller = AnalysisSelectionController().bind_report(report)

    with pytest.raises((IndexError, KeyError, ValueError), match="alignment|correspondence|report"):
        controller.select_alignment_index(999_999)
    with pytest.raises((IndexError, KeyError, ValueError), match="evidence|alignment|correspondence"):
        controller.lookup_evidence_card(alignment_index=999_999)

    other_directory = tmp_path / "other"
    other_directory.mkdir()
    other_request = _request(other_directory)
    other_request = replace(
        other_request,
        analysis_settings=replace(
            other_request.analysis_settings,
            minimum_sequence_identity=0.25,
        ),
    )
    other = ReportService().analyze(other_request)
    assert other.report_id != report.report_id

    rebound = controller.bind_report(other)

    assert rebound.report_id == other.report_id
    assert rebound.selection == AnalysisSelection()
    assert controller.report_id == report.report_id


def test_selection_controller_reports_unbound_empty_and_type_errors() -> None:
    from structlens.plugin.gui.selection_controller import AnalysisSelectionController

    controller = AnalysisSelectionController()
    with pytest.raises(ValueError, match="report|bound"):
        controller.lookup_correspondence(0)
    with pytest.raises(ValueError, match="report|bound"):
        controller.lookup_evidence_card(0)
    with pytest.raises(ValueError, match="report|bound"):
        controller.select_evidence_card(ResidueId("reference", "1", "A", "1", None, "ALA"))
    with pytest.raises(TypeError, match="report"):
        controller.bind_report(object())


def test_selection_controller_validates_correspondence_and_evidence_focus() -> None:
    from structlens.plugin.gui.selection_controller import AnalysisSelectionController

    reference = ResidueId("reference", "1", "A", "1", None, "ALA")
    target = ResidueId("target", "1", "A", "1", None, "ALA")
    analysis = AnalysisResult(
        reference_id="reference",
        target_id="target",
        correspondences=(
            ResidueCorrespondence(
                alignment_index=3,
                reference=reference,
                target=target,
                reference_one_letter="A",
                target_one_letter="A",
                status=CorrespondenceStatus.CONSERVED,
            ),
        ),
        mutations=(),
        sequence_identity=1.0,
        sequence_coverage=1.0,
        alignment_decision="sequence",
    )
    report = _minimal_report(analysis=analysis)
    bound = AnalysisSelectionController().bind_report(report)

    with pytest.raises(TypeError, match="alignment"):
        bound.select_alignment_index("3")
    with pytest.raises(TypeError, match="alignment"):
        bound.lookup_correspondence("3")
    with pytest.raises(TypeError, match="alignment"):
        bound.lookup_evidence_card("3")
    with pytest.raises(IndexError, match="alignment index"):
        bound.lookup_correspondence(99)
    with pytest.raises(IndexError, match="alignment index"):
        bound.select_alignment_index(99)
    with pytest.raises(IndexError, match="alignment index"):
        bound.lookup_evidence_card(99)
    with pytest.raises(KeyError, match="evidence card"):
        bound.select_evidence_card(reference)
    with pytest.raises(KeyError, match="evidence card"):
        bound.lookup_evidence_card(3)

    foreign = ResidueId("foreign", "1", "A", "1", None, "ALA")
    foreign_focus = replace(bound, selection=AnalysisSelection(reference_residue=foreign))
    with pytest.raises(ValueError, match="selected correspondence"):
        foreign_focus.lookup_correspondence()
    with pytest.raises(ValueError, match="no correspondence"):
        bound.clear_selection().lookup_correspondence()


def test_selection_controller_rejects_analysis_free_reports() -> None:
    from structlens.plugin.gui.selection_controller import AnalysisSelectionController

    bound = AnalysisSelectionController().bind_report(_minimal_report())
    with pytest.raises(ValueError, match="analysis correspondences"):
        bound.select_alignment_index(0)
    with pytest.raises(ValueError, match="analysis correspondences"):
        bound.lookup_correspondence(0)
    with pytest.raises(ValueError, match="analysis correspondences"):
        bound.lookup_evidence_card(0)
