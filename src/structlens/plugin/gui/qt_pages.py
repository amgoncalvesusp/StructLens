"""Presentation-only widget and page helpers used by the Qt controller."""

from __future__ import annotations

from typing import Any

from structlens.core.models import AnalysisResult, ResidueId
from structlens.plugin.visualization.renderer import VisualizationState


def page_subtitle(section: str, *, standalone: bool = False) -> str:
    return {
        "Project": "Choose sources and chains.",
        "Sequences": "Review sequence identity, mutations, and conservation.",
        "Structures": "Choose structural comparison and alignment settings.",
        "Residues": "Inspect mapped positions.",
        "Sites": "Define and compare active-site or ligand-centered regions.",
        "Charts": (
            "Inspect chart-ready scientific profiles and visual filters."
            if standalone
            else "Apply a reversible PyMOL view from linked chart selections."
        ),
        "PyMOL": "Prepare a validated interchange bundle for the companion plugin.",
        "Results": "Review and export metrics.",
        "Export": "Write evidence tables, chart data, and interchange bundles.",
    }[section]


def label(widgets: Any, text: str, object_name: str) -> Any:
    widget = widgets.QLabel(text)
    widget.setObjectName(object_name)
    if object_name in {"fieldMeta", "helpText", "inlineNote", "pagePurpose", "resultDecision", "legend"}:
        widget.setWordWrap(True)
    return widget


def button(widgets: Any, text: str, object_name: str) -> Any:
    widget = widgets.QPushButton(text)
    widget.setObjectName(object_name)
    return widget


def fraction_spin(spin: Any, value: float) -> None:
    spin.setRange(0.0, 1.0)
    spin.setSingleStep(0.05)
    spin.setDecimals(2)
    spin.setValue(value)
    spin.setSuffix(" fraction")


def configure_table(table: Any, widgets: Any) -> None:
    table.setAlternatingRowColors(True)
    table.setSortingEnabled(False)
    view = widgets.QAbstractItemView
    table.setSelectionBehavior(getattr(view, "SelectionBehavior", view).SelectRows)
    table.setSelectionMode(getattr(view, "SelectionMode", view).SingleSelection)
    table.setEditTriggers(getattr(view, "EditTrigger", view).NoEditTriggers)
    table.verticalHeader().setVisible(False)
    table.horizontalHeader().setStretchLastSection(True)
    table.setMinimumHeight(220)


def residue_label(residue: ResidueId | None) -> str:
    if residue is None:
        return "—"
    insertion = residue.insertion_code or ""
    return f"{residue.chain_id}:{residue.auth_seq_id}{insertion} {residue.residue_name}"


def number(value: float | int | None) -> str:
    return "—" if value is None else f"{value:.3f}" if isinstance(value, float) else str(value)


def fraction_number(value: float | None) -> str:
    return "—" if value is None else f"{value:.3f}"


def backend_label(result: AnalysisResult) -> str:
    backend = result.provenance.get("backend") or result.provenance.get("structural_backend")
    return str(backend) if backend else "sequence-guided"


def human(value: str) -> str:
    return value.replace("_", " ").replace("ca", "Cα").title()


def state_dict(state: VisualizationState) -> dict[str, Any]:
    return {
        "highlight_filter": state.highlight_filter.value,
        "color_mode": state.color_mode.value,
        "representation": state.representation.value,
        "show_labels": state.show_labels,
        "show_reference": state.show_reference,
        "show_target": state.show_target,
        "local_radius_angstrom": state.local_radius_angstrom,
        "preset": state.preset,
    }


def stylesheet() -> str:
    """Return the panel's shared visual language."""

    return """
    QWidget#structlensPanel { background: #111827; color: #dbe7f3; }
    QScrollArea, QWidget[pageCanvas="true"] { background: #111827; color: #dbe7f3; }
    QFrame#header { background: #0c1421; border-bottom: 1px solid #26364b; }
    QFrame#sidebar { background: #0d1726; border-right: 1px solid #26364b; }
    QFrame#footer { background: #0c1421; border-top: 1px solid #26364b; }
    QLabel#eyebrow, QLabel#sectionKicker, QLabel#fieldLabel { color: #7f9ab8; font-size: 10px; font-weight: 700; }
    QLabel#eyebrow, QLabel#sectionKicker, QLabel#fieldLabel { letter-spacing: 1px; }
    QLabel#windowTitle { color: #f5f8fc; font-size: 20px; font-weight: 700; }
    QLabel#pageTitle { color: #f5f8fc; font-size: 24px; font-weight: 700; }
    QLabel#pagePurpose { color: #9eb0c5; font-size: 13px; }
    QLabel#workflowContext { color: #b8c9dc; font-size: 12px; }
    QLabel#statusPill { background: #19304a; color: #9fc7ff; border: 1px solid #2a5a8a; border-radius: 12px; padding: 5px 11px; font-weight: 700; }
    QLabel#footerStatus { color: #b8c9dc; }
    QLabel#footerMeta, QLabel#fieldMeta, QLabel#sidebarNote { color: #71869e; }
    QLabel#sidebarNote { font-size: 11px; line-height: 1.3; }
    QLabel#helpText, QLabel#inlineNote { color: #a8bbd0; line-height: 1.35; }
    QLabel#legend { background: #172438; color: #c3d6ea; border: 1px solid #2e435d; padding: 11px; }
    QLabel#resultDecision { background: #14253a; color: #d9e8f8; border: 1px solid #2a537d; padding: 13px; }
    QLabel#metricValue { color: #f5f8fc; font-size: 17px; font-weight: 700; }
    QListWidget#workflowNav { background: transparent; border: none; outline: none; color: #aabbd0; }
    QListWidget#workflowNav::item { padding: 6px 8px; border-radius: 6px; }
    QListWidget#workflowNav::item:hover { background: #16283e; color: #eaf3ff; }
    QListWidget#workflowNav::item:selected { background: #1e4d78; color: #ffffff; font-weight: 700; }
    QGroupBox { background: #141f30; border: 1px solid #293c55; border-radius: 5px; margin-top: 8px; padding-top: 15px; font-weight: 700; color: #e0ebf7; }
    QGroupBox::title { subcontrol-origin: margin; left: 13px; padding: 0 6px; color: #c6d8eb; }
    QLineEdit, QComboBox, QDoubleSpinBox, QPlainTextEdit { background: #0f1a2b; color: #e4eef9; border: 1px solid #344b67; border-radius: 4px; padding: 8px; selection-background-color: #2f7af8; }
    QLineEdit:focus, QComboBox:focus, QDoubleSpinBox:focus, QPlainTextEdit:focus { border: 1px solid #6ba9f8; }
    QComboBox QAbstractItemView { background: #172438; color: #e4eef9; selection-background-color: #24578a; }
    QPushButton { min-height: 32px; padding: 0 15px; border-radius: 4px; font-weight: 700; }
    QPushButton#primaryButton { background: #2f7af8; color: #ffffff; border: 1px solid #6aa4ff; }
    QPushButton#primaryButton:hover { background: #4b8cff; }
    QPushButton#primaryButton:pressed { background: #1f60cf; }
    QPushButton#primaryButton:disabled { background: #29415d; color: #7990a9; border-color: #29415d; }
    QPushButton#secondaryButton { background: #1b2a3e; color: #d8e6f5; border: 1px solid #3a516d; }
    QPushButton#secondaryButton:hover { background: #243a54; }
    QCheckBox { color: #c3d4e7; spacing: 8px; }
    QTableWidget { background: #0f1a2b; alternate-background-color: #142236; color: #dce9f6; border: 1px solid #2c4059; gridline-color: #22354c; selection-background-color: #24578a; selection-color: #ffffff; }
    QHeaderView::section { background: #1a2b41; color: #a9c0d9; padding: 8px; border: none; border-right: 1px solid #2c4059; font-weight: 700; }
    QProgressBar { background: #172438; border: 1px solid #2d4662; height: 8px; text-align: center; color: transparent; }
    QProgressBar::chunk { background: #2f7af8; }
    QScrollBar:vertical { background: #0e1828; width: 11px; }
    QScrollBar::handle:vertical { background: #35516e; min-height: 28px; }
    QScrollBar:horizontal { background: #0e1828; height: 11px; }
    QScrollBar::handle:horizontal { background: #35516e; min-width: 28px; }
    """


__all__ = [
    "backend_label",
    "button",
    "configure_table",
    "fraction_number",
    "fraction_spin",
    "human",
    "label",
    "number",
    "page_subtitle",
    "residue_label",
    "state_dict",
    "stylesheet",
]
