"""Small compatibility facade for the StructLens Qt panel.

The controller implementation lives in :mod:`panel_controller`; this module
keeps the historical import path stable for PyMOL and test integrations.
"""

from __future__ import annotations

from structlens.application.chart_export import export_chart_image, export_chart_xlsx
from structlens.application.export_service import (
    export_analysis_csv,
    export_analysis_json,
    export_analysis_tsv,
    export_analysis_xlsx,
    export_v03_xlsx,
)
from structlens.integrations.pymol_bundle import write_pymol_bundle

from .main_panel import SCIENTIFIC_SECTIONS, WORKFLOW_HELP
from .panel_controller import PanelController, build_panel
from .qt_chart import matrix_image_kwargs as _matrix_image_kwargs
from .qt_compat import QtBindings
from .report_controller import AnalysisReportController

__all__ = [
    "PanelController",
    "AnalysisReportController",
    "SCIENTIFIC_SECTIONS",
    "QtBindings",
    "WORKFLOW_HELP",
    "build_panel",
    "export_chart_image",
    "export_chart_xlsx",
    "export_v03_xlsx",
    "export_analysis_csv",
    "export_analysis_json",
    "export_analysis_tsv",
    "export_analysis_xlsx",
    "write_pymol_bundle",
    "_matrix_image_kwargs",
]
