"""Dependency-free immutable state and workflow labels for the GUI."""

from __future__ import annotations

from dataclasses import dataclass, replace

from structlens.application.dto import AnalysisReportRequest
from structlens.application.report_snapshot_io import verify_report_inputs
from structlens.core.models import AnalysisResult, AnalysisSelection
from structlens.core.parsing import SourceSnapshot
from structlens.core.reports import AnalysisReport
from structlens.plugin.visualization.renderer import VisualizationState

from .presentation import ReportPresentation, present_report

GUI_SECTIONS = (
    "Project",
    "Alignment",
    "Mutations",
    "Residues",
    "Visualization",
    "Results",
)

SCIENTIFIC_SECTIONS = (
    "Project",
    "Sequences",
    "Structures",
    "Residues",
    "Sites",
    "Charts",
    "PyMOL",
    "Results",
    "Export",
)

WORKFLOW_HELP = {
    "Auto": "What it does: evaluates sequence identity and coverage, then chooses sequence or structure mapping. Use it for most comparisons; the decision is reported.",
    "Sequence": "What it does: builds global amino-acid correspondence before superposition. Use it for homologous proteins with meaningful sequence similarity.",
    "Structure": "What it does: proposes correspondence from structural similarity through US-align. Use it when sequence identity is low; pairs remain inspectable.",
    "Manual": "What it does: accepts explicit locked residue pairs. Use it when biologically validated correspondence must be enforced.",
    "sequence identity": "What it means: the fraction of aligned canonical residues with the same amino acid.",
    "alignment coverage": "What it means: the fraction of the reference canonical residues represented by mapped canonical pairs.",
    "TM-score": "What it means: a length-normalized structural similarity score when US-align produces it.",
    "strict RMSD": "What it means: RMSD over every eligible mapped Cα pair before any refinement exclusions.",
    "refined RMSD": "What it means: RMSD after deterministic cutoff-based refinement; excluded correspondences remain listed as outliers.",
    "Cα displacement": "What it means: the aligned distance between one reference and target Cα pair; it is not a residue RMSD.",
    "backbone RMSD": "What it means: RMSD over matched N, Cα, C and O atoms for one mapped residue.",
    "side-chain RMSD": "What it means: symmetry-aware RMSD over matched side-chain heavy atoms.",
    "BLOSUM62": "What it means: a substitution descriptor from the embedded BLOSUM62 matrix; it does not imply function.",
    "Grantham distance": "What it means: a physicochemical substitution-distance descriptor; it does not imply pathogenicity.",
}


@dataclass(frozen=True, slots=True)
class CanonicalReportBinding:
    """One verified authority for canonical report consumers."""

    report: AnalysisReport
    request: AnalysisReportRequest
    snapshots: tuple[SourceSnapshot, SourceSnapshot]
    presentation: ReportPresentation

    @classmethod
    def create(
        cls,
        report: AnalysisReport,
        request: AnalysisReportRequest,
        presentation: ReportPresentation | None = None,
    ) -> CanonicalReportBinding:
        if not isinstance(report, AnalysisReport):
            raise TypeError("report must be an AnalysisReport")
        if not isinstance(request, AnalysisReportRequest):
            raise TypeError("request must be an AnalysisReportRequest")
        snapshots = (request.reference_snapshot, request.target_snapshot)
        verify_report_inputs(report, snapshots, (), None)
        if report.reference_selection.selection_id != request.reference_selection.selection_id:
            raise ValueError("reference request selection does not match the report")
        if report.target_selection.selection_id != request.target_selection.selection_id:
            raise ValueError("target request selection does not match the report")
        projection = presentation if presentation is not None else present_report(report)
        if projection.report_id != report.report_id:
            raise ValueError("presentation does not match the report ID")
        return cls(report, request, snapshots, projection)


@dataclass(frozen=True, slots=True)
class StructLensPanelModel:
    """Immutable UI state kept separate from the authoritative analysis result."""

    status: str = "Choose a reference structure and a target to begin."
    report: AnalysisReport | None = None
    binding: CanonicalReportBinding | None = None
    analysis: AnalysisResult | None = None
    selection: AnalysisSelection | None = None
    reference_path: str | None = None
    target_path: str | None = None
    reference_chain_id: str | None = None
    target_chain_id: str | None = None
    visualization_state: VisualizationState = VisualizationState()
    busy: bool = False
    error: str | None = None

    def with_status(self, status: str) -> StructLensPanelModel:
        return replace(self, status=status, busy=False, error=None)

    def with_analysis(self, analysis: AnalysisResult) -> StructLensPanelModel:
        return replace(
            self,
            report=None,
            binding=None,
            analysis=analysis,
            selection=None,
            status="Analysis complete.",
            busy=False,
            error=None,
        )

    def with_report(self, report: AnalysisReport) -> StructLensPanelModel:
        if not isinstance(report, AnalysisReport):
            raise TypeError("report must be an AnalysisReport")
        return replace(
            self,
            report=report,
            binding=None,
            analysis=None,
            selection=None,
            status="Analysis report complete.",
            busy=False,
            error=None,
        )

    def with_binding(self, binding: CanonicalReportBinding) -> StructLensPanelModel:
        if not isinstance(binding, CanonicalReportBinding):
            raise TypeError("binding must be a CanonicalReportBinding")
        return replace(
            self,
            report=binding.report,
            binding=binding,
            analysis=None,
            selection=None,
            status="Analysis report complete.",
            busy=False,
            error=None,
        )

    def with_sources(
        self,
        *,
        reference_path: str | None = None,
        target_path: str | None = None,
        reference_chain_id: str | None = None,
        target_chain_id: str | None = None,
    ) -> StructLensPanelModel:
        return replace(
            self,
            reference_path=reference_path,
            target_path=target_path,
            reference_chain_id=reference_chain_id,
            target_chain_id=target_chain_id,
            error=None,
        )

    def with_visualization(self, visualization_state: VisualizationState) -> StructLensPanelModel:
        return replace(self, visualization_state=visualization_state)

    def with_busy(self, status: str) -> StructLensPanelModel:
        return replace(self, status=status, busy=True, error=None)

    def with_error(self, message: str) -> StructLensPanelModel:
        return replace(self, status=message, busy=False, error=message)


__all__ = [
    "CanonicalReportBinding",
    "GUI_SECTIONS",
    "SCIENTIFIC_SECTIONS",
    "WORKFLOW_HELP",
    "StructLensPanelModel",
]
