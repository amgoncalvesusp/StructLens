"""Qt composition and controller for the StructLens PyMOL workflow."""

# ruff: noqa: F403,F405

from __future__ import annotations

from .qt_context import *  # noqa: F401,F403
from .qt_exports_mixin import ExportMixin
from .qt_pages_mixin import PageMixin
from .qt_presentation_mixin import PresentationMixin
from .qt_reports_mixin import ReportMixin
from .qt_shell_mixin import ShellMixin
from .qt_sources_mixin import SourceMixin
from .qt_visualization_mixin import VisualizationMixin


def build_panel(
    qt: QtBindings,
    *,
    parent: object | None = None,
    command: object | None = None,
) -> object:
    """Create the panel and retain its controller on the QWidget for host access."""

    controller = PanelController(qt, parent=parent, command=command)
    controller.widget._structlens_controller = controller
    controller.widget.destroyed.connect(lambda: controller.close())
    return controller.widget


class PanelController(
    ShellMixin,
    PageMixin,
    PresentationMixin,
    SourceMixin,
    VisualizationMixin,
    ReportMixin,
    ExportMixin,
):
    """Own Qt widgets and side effects while keeping the model immutable."""

    _report_controller: AnalysisReportController

    def __init__(
        self,
        qt: QtBindings,
        *,
        parent: object | None,
        command: object | None,
    ) -> None:
        self.qt = qt
        self.w = qt.widgets
        self.c = qt.core
        self.g = qt.gui
        self.command = command
        self.widget = self.w.QWidget(parent)
        self.widget.setObjectName("structlensPanel")
        self.widget.setWindowTitle("StructLens · Evidence Bench")
        self.widget.setMinimumSize(840, 600)
        icon_path = Path(__file__).parents[1] / "assets" / "structlens_icon.png"
        if icon_path.exists():
            self.widget.setWindowIcon(self.g.QIcon(str(icon_path)))
        self.model = StructLensPanelModel()
        self.reference_structure: ProteinStructure | None = None
        self.target_structure: ProteinStructure | None = None
        self.reference_loaded_source: LoadedSource | None = None
        self.target_loaded_source: LoadedSource | None = None
        self.reference_object_name: str | None = None
        self.target_object_name: str | None = None
        self._temporary_paths: list[Path] = []
        object.__setattr__(self, "_report_controller", AnalysisReportController(view=self))
        # The accepted request remains authoritative while a new Compare is
        # pending.  Only a fully validated result promotes the pending value.
        self._report_request: AnalysisReportRequest | None = None
        self._pending_report_request: AnalysisReportRequest | None = None
        self._report_presentation: ReportPresentation | None = None
        self._report_history: tuple[AnalysisReport, ...] = ()
        self._presentation_history: tuple[ReportPresentation, ...] = ()
        self._selection_controller = AnalysisSelectionController()
        self._renderer = VisualizationRenderer()
        self._visualization_service = VisualizationService(self._renderer)
        self._pymol = PyMOLAdapter(command, project_id="panel")
        self._future: Any = None
        self._poll_timer: Any = None
        self._cancel_event = Event()
        self.chart_export_buttons: list[Any] = []
        # v0.3 scientific services own these calculations.  The GUI only
        # retains JSON-ready payloads supplied by the orchestration layer so
        # the same authoritative values can be handed to PyMOL.
        self._v03_bundle_payloads: dict[str, Mapping[str, Any] | None] = {
            "msa_summary": None,
            "conservation": None,
            "interactions": None,
            "sites": None,
            "evidence": None,
            "vectors": None,
        }
        self._v03_export_records: dict[str, Any] = {}
        self._chart_datasets: dict[str, ChartDataset | MatrixDataset] = {}
        self._analysis_history: tuple[AnalysisResult, ...] = ()
        self._site_metrics: tuple[SiteMetrics, ...] = ()
        self._site_definitions: tuple[SiteDefinition, ...] = ()
        self._sequence_chart_canvas: Any = None
        self._chart_canvas: Any = None
        self._msa_chart_dataset: ChartDataset | None = None
        self._build_shell()
        self._build_project_page()
        self._build_mutations_page()
        self._build_alignment_page()
        self._build_residues_page()
        self._build_sites_page()
        self._build_visualization_page()
        self._build_pymol_page()
        self._build_results_page()
        self._build_export_page()
        self._wire_navigation()
        self._set_status(self.model.status)
        self._update_mode_help(self.mode_combo.currentText())
        self._update_comparison_help(self.comparison_combo.currentText())
        self._update_chart_explanation(self.chart_combo.currentText())
        self.usalign_status.setText(
            "Bundled backend · Ready"
            if bundled_executable() is not None
            else "Bundled backend · not present in this source build; custom/PATH fallback available"
        )

    # ------------------------------------------------------------------ shell


__all__ = ["PanelController", "build_panel"]
