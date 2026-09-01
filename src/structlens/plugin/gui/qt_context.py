"""Shared, binding-neutral symbols for the focused Qt controller modules."""

# This module intentionally re-exports composition dependencies for mixins;
# they are consumed by the concrete PanelController at runtime.
# ruff: noqa: F401

from __future__ import annotations

import re
import sys
from collections.abc import Iterable, Mapping
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
from pathlib import Path
from tempfile import NamedTemporaryFile
from threading import Event
from typing import Any

from structlens.application.chart_data import (
    ChartDataset,
    ChartSeries,
    MatrixDataset,
    mutation_conservation_matrix,
    structural_deviation_profile,
)
from structlens.application.chart_export import export_chart_image, export_chart_xlsx
from structlens.application.dto import AnalysisReportRequest
from structlens.application.export_service import (
    export_analysis_csv,
    export_analysis_json,
    export_analysis_tsv,
    export_analysis_xlsx,
    export_v03_xlsx,
)
from structlens.application.project_state import ProjectState
from structlens.application.visualization_service import VisualizationService
from structlens.core.errors import AnalysisCancelledError, BundleValidationError, ProjectSchemaError
from structlens.core.models import (
    AlignmentMode,
    AnalysisResult,
    AnalysisSettings,
    ComparisonMode,
    ProteinChain,
    ProteinStructure,
    ResidueId,
)
from structlens.core.msa import MultipleSequenceAlignment
from structlens.core.parsing import InputSelection, SnapshotError, load_structure
from structlens.core.reports import AnalysisReport
from structlens.core.sites import SiteDefinition, SiteDefinitionMode, SiteMetrics
from structlens.integrations.pymol.adapter import PyMOLAdapter
from structlens.integrations.pymol.launcher import PyMOLLauncher
from structlens.integrations.pymol.selections import selection_name
from structlens.integrations.pymol_bundle import write_pymol_bundle
from structlens.integrations.usalign.executable import bundled_executable
from structlens.plugin.visualization.legends import backbone_rmsd_legend, displacement_legend
from structlens.plugin.visualization.renderer import (
    ColorMode,
    HighlightFilter,
    Representation,
    VisualizationRenderer,
    VisualizationState,
)
from structlens.resources.backends import backend_versions

from .loaded_source import LoadedSource
from .model import SCIENTIFIC_SECTIONS, WORKFLOW_HELP, CanonicalReportBinding, StructLensPanelModel
from .presentation import ReportPresentation, freeze_json, present_report, thaw_json
from .project_transaction import (
    CanonicalProjectCandidate,
    LegacyProjectCandidate,
    parse_visualization_state,
    save_canonical_project,
    stage_project,
)
from .qt_chart import load_chart_classes, matrix_image_kwargs
from .qt_compat import QtBindings
from .qt_pages import (
    backend_label,
    button,
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
from .qt_reports import legacy_analysis, report_chart_datasets, report_site_metrics
from .qt_sources import combo_data, find_residue, selected_chain, structure_meta
from .report_controller import AnalysisReportController
from .selection_controller import AnalysisSelectionController


class QtMixinContext:
    """Typing/runtime seam for mixins whose widgets are injected by Qt."""

    model: StructLensPanelModel
    reference_structure: ProteinStructure | None
    target_structure: ProteinStructure | None
    reference_loaded_source: LoadedSource | None
    target_loaded_source: LoadedSource | None
    reference_object_name: str | None
    target_object_name: str | None
    _temporary_paths: list[Path]
    _report_controller: AnalysisReportController
    _report_request: AnalysisReportRequest | None
    _pending_report_request: AnalysisReportRequest | None
    _report_presentation: ReportPresentation | None
    _report_history: tuple[AnalysisReport, ...]
    _presentation_history: tuple[ReportPresentation, ...]
    _selection_controller: AnalysisSelectionController
    _renderer: VisualizationRenderer
    _visualization_service: VisualizationService
    _pymol: PyMOLAdapter
    _future: Any
    _poll_timer: Any
    _cancel_event: Event
    chart_export_buttons: list[Any]
    _v03_bundle_payloads: dict[str, Mapping[str, Any] | None]
    _v03_export_records: dict[str, Any]
    _chart_datasets: dict[str, ChartDataset | MatrixDataset]
    _analysis_history: tuple[AnalysisResult, ...]
    _site_metrics: tuple[SiteMetrics, ...]
    _site_definitions: tuple[SiteDefinition, ...]
    _sequence_chart_canvas: Any
    _chart_canvas: Any
    _msa_chart_dataset: ChartDataset | None

    def __getattr__(self, name: str) -> Any:
        raise AttributeError(name)


_SEQUENCES_PAGE_INDEX = SCIENTIFIC_SECTIONS.index("Sequences")
_STRUCTURES_PAGE_INDEX = SCIENTIFIC_SECTIONS.index("Structures")
_CHARTS_PAGE_INDEX = SCIENTIFIC_SECTIONS.index("Charts")
_ANALYSIS_EXECUTOR = ThreadPoolExecutor(max_workers=1, thread_name_prefix="structlens-analysis")
_RESULTS_PAGE_INDEX = SCIENTIFIC_SECTIONS.index("Results")


def compat_symbol(name: str, fallback: Any) -> Any:
    """Resolve a patched legacy export through the stable facade when present."""

    facade = sys.modules.get(f"{__package__}.qt_panel")
    return getattr(facade, name, fallback) if facade is not None else fallback


# Method bodies in the focused mixins retain the historical private helper
# names so downstream callers that introspect the controller remain stable.
_compat_symbol = compat_symbol
_backend_label = backend_label
_button = button
_combo_data = combo_data
_configure_table = configure_table
_find_residue = find_residue
_fraction_number = fraction_number
_fraction_spin = fraction_spin
_human = human
_label = label
_legacy_analysis = legacy_analysis
_matrix_image_kwargs = matrix_image_kwargs
_number = number
_page_subtitle = page_subtitle
_report_chart_datasets = report_chart_datasets
_report_site_metrics = report_site_metrics
_residue_label = residue_label
_selected_chain_for_gui = selected_chain
_state_dict = state_dict
_structure_meta = structure_meta
_stylesheet = stylesheet


__all__ = [
    "CanonicalReportBinding",
    "CanonicalProjectCandidate",
    "freeze_json",
    "InputSelection",
    "LoadedSource",
    "LegacyProjectCandidate",
    "ReportPresentation",
    "SnapshotError",
    "save_canonical_project",
    "stage_project",
    "AnalysisCancelledError",
    "AnalysisReport",
    "AnalysisReportController",
    "AnalysisReportRequest",
    "AnalysisResult",
    "AnalysisSelectionController",
    "AnalysisSettings",
    "AlignmentMode",
    "Any",
    "backend_label",
    "backbone_rmsd_legend",
    "backend_versions",
    "BundleValidationError",
    "button",
    "ChartDataset",
    "ChartSeries",
    "ColorMode",
    "combo_data",
    "ComparisonMode",
    "compat_symbol",
    "configure_table",
    "displacement_legend",
    "Event",
    "export_analysis_csv",
    "export_analysis_json",
    "export_analysis_tsv",
    "export_analysis_xlsx",
    "export_chart_image",
    "export_chart_xlsx",
    "export_v03_xlsx",
    "find_residue",
    "fraction_number",
    "fraction_spin",
    "human",
    "HighlightFilter",
    "label",
    "load_chart_classes",
    "load_structure",
    "legacy_analysis",
    "Mapping",
    "matrix_image_kwargs",
    "MatrixDataset",
    "mutation_conservation_matrix",
    "MultipleSequenceAlignment",
    "NamedTemporaryFile",
    "number",
    "page_subtitle",
    "Path",
    "ProjectState",
    "ProjectSchemaError",
    "ProteinChain",
    "ProteinStructure",
    "PyMOLAdapter",
    "PyMOLLauncher",
    "QtBindings",
    "re",
    "replace",
    "Representation",
    "report_chart_datasets",
    "report_site_metrics",
    "residue_label",
    "ResidueId",
    "SCIENTIFIC_SECTIONS",
    "selected_chain",
    "selection_name",
    "SiteDefinition",
    "SiteDefinitionMode",
    "SiteMetrics",
    "state_dict",
    "stylesheet",
    "structural_deviation_profile",
    "structure_meta",
    "sys",
    "ThreadPoolExecutor",
    "VisualizationRenderer",
    "VisualizationService",
    "VisualizationState",
    "WORKFLOW_HELP",
    "write_pymol_bundle",
    "bundled_executable",
    "QtMixinContext",
    "_ANALYSIS_EXECUTOR",
    "_CHARTS_PAGE_INDEX",
    "_RESULTS_PAGE_INDEX",
    "_SEQUENCES_PAGE_INDEX",
    "_STRUCTURES_PAGE_INDEX",
    "_backend_label",
    "_button",
    "_combo_data",
    "_compat_symbol",
    "_configure_table",
    "_find_residue",
    "_fraction_number",
    "_fraction_spin",
    "_human",
    "_label",
    "_legacy_analysis",
    "_matrix_image_kwargs",
    "_number",
    "_page_subtitle",
    "_report_chart_datasets",
    "_report_site_metrics",
    "_residue_label",
    "_selected_chain_for_gui",
    "_state_dict",
    "_structure_meta",
    "_stylesheet",
    "StructLensPanelModel",
    "thaw_json",
]
