from dataclasses import dataclass

from structlens.core.evidence import Availability
from structlens.core.parsing import InputSelection, StructureFormat
from structlens.core.quality import StructureQualityReport
from structlens.core.reports import AnalysisReport, InputQualityBundle, SectionAvailability
from structlens.plugin.gui.main_panel import (
    GUI_SECTIONS,
    WORKFLOW_HELP,
    AnalysisReportController,
    StructLensPanelModel,
)


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
