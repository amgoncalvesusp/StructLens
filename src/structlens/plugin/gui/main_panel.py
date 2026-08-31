"""Operate-mode Qt panel for the optional StructLens PyMOL plugin.

Impeccable direction contract — Evidence Bench
THESIS: make the residue correspondence table the primary workspace, refusing a
decorative dashboard of unexplained metrics.
OWN-WORLD: host-compatible graphite surfaces, cool blue action states, compact
system sans controls, and dense tabular data with explicit Å/fraction units.
STORY: the user loads two sources, chooses an alignment policy, runs one
reproducible comparison, then moves from evidence to a reversible PyMOL view.
FIRST VIEWPORT: a narrow workflow rail anchors a two-column Project page; the
header keeps status and Compare visible while the right canvas holds source
controls, chain choices, and the next action.
FORM: Operate-mode split workspace, chosen to keep eight scientific stages visible
without burying the task in a tab strip or card mosaic.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Protocol

from structlens.application.dto import AnalysisReportRequest
from structlens.core.models import AnalysisResult
from structlens.core.reports import AnalysisReport
from structlens.plugin.visualization.renderer import VisualizationState

GUI_SECTIONS = (
    "Project",
    "Alignment",
    "Mutations",
    "Residues",
    "Visualization",
    "Results",
)

# v0.2 scientific labels are kept separately so host integrations that relied
# on the v0.1 menu identifiers remain source-compatible.
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
class StructLensPanelModel:
    """Immutable UI state kept separate from the authoritative analysis result."""

    status: str = "Choose a reference structure and a target to begin."
    report: AnalysisReport | None = None
    analysis: AnalysisResult | None = None
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
        """Retain the v0.3 compatibility path for callers not yet report-driven."""

        return replace(
            self,
            report=None,
            analysis=analysis,
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
            analysis=None,
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

    def with_visualization(
        self, visualization_state: VisualizationState
    ) -> StructLensPanelModel:
        return replace(self, visualization_state=visualization_state)

    def with_busy(self, status: str) -> StructLensPanelModel:
        return replace(self, status=status, busy=True, error=None)

    def with_error(self, message: str) -> StructLensPanelModel:
        return replace(self, status=message, busy=False, error=message)


class _ReportRunner(Protocol):
    def analyze(self, request: AnalysisReportRequest) -> AnalysisReport: ...


class _ReportView(Protocol):
    def show_report(self, report: AnalysisReport) -> None: ...

    def show_status(self, status: str) -> None: ...


class AnalysisReportController:
    """Headless GUI boundary that moves one typed report into one view."""

    def __init__(self, service: _ReportRunner, view: _ReportView) -> None:
        self._service = service
        self._view = view
        self.model = StructLensPanelModel()

    def analyze(self, request: AnalysisReportRequest) -> AnalysisReport:
        self.model = self.model.with_busy("Building canonical analysis report…")
        report = self._service.analyze(request)
        self.model = self.model.with_report(report)
        self._view.show_report(report)
        self._view.show_status(f"Analysis report complete · {report.report_id}")
        return report


def build_qt_panel(
    parent: object | None = None,
    *,
    command: object | None = None,
) -> object:
    """Build the full panel when a supported Qt binding is present."""

    from .qt_compat import load_qt

    qt = load_qt()
    if qt is None:
        raise RuntimeError(
            "StructLens GUI requires PySide6 or PyQt5; install structlens[gui] "
            "or launch it inside a PyMOL environment that provides Qt"
        )
    from .qt_panel import build_panel

    return build_panel(qt, parent=parent, command=command)


__all__ = [
    "GUI_SECTIONS",
    "SCIENTIFIC_SECTIONS",
    "WORKFLOW_HELP",
    "AnalysisReportController",
    "StructLensPanelModel",
    "build_qt_panel",
]
