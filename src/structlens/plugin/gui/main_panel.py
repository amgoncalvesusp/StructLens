"""Optional Qt panel with a table-first, six-section scientific workflow."""

from __future__ import annotations

import importlib
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

from structlens.core.models import AnalysisResult

GUI_SECTIONS = (
    "Project",
    "Alignment",
    "Mutations",
    "Residues",
    "Visualization",
    "Results",
)

WORKFLOW_HELP = {
    "Auto": "What it does: evaluates sequence identity and coverage, then chooses sequence or structure mapping. Use it for most comparisons; the decision is reported.",
    "Sequence": "What it does: builds global amino-acid correspondence before superposition. Use it for homologous proteins with meaningful sequence similarity.",
    "Structure": "What it does: proposes correspondence from structural similarity through US-align. Use it when sequence identity is low; pairs remain inspectable.",
    "Manual": "What it does: accepts explicit locked residue pairs. Use it when biologically validated correspondence must be enforced.",
    "sequence identity": "What it means: the fraction of aligned canonical residues with the same amino acid.",
    "alignment coverage": "What it means: the fraction of the reference canonical residues represented by mapped canonical pairs.",
    "TM-score": "What it means: a length-normalized structural similarity score when US-align produces it.",
    "strict RMSD": "What it means: RMSD over every eligible mapped atom before any refinement exclusions.",
    "refined RMSD": "What it means: RMSD after deterministic cutoff-based refinement; excluded correspondences remain listed as outliers.",
    "Cα displacement": "What it means: the aligned distance between one reference and target Cα pair; it is not a residue RMSD.",
    "backbone RMSD": "What it means: RMSD over matched N, Cα, C and O atoms for one mapped residue.",
    "side-chain RMSD": "What it means: symmetry-aware RMSD over matched side-chain heavy atoms.",
    "BLOSUM62": "What it means: a substitution descriptor from the embedded BLOSUM62 matrix; it does not imply function.",
    "Grantham distance": "What it means: a physicochemical substitution-distance descriptor; it does not imply pathogenicity.",
}


@dataclass(frozen=True, slots=True)
class StructLensPanelModel:
    status: str = "Choose a reference structure and a target to begin."
    analysis: AnalysisResult | None = None

    def with_status(self, status: str) -> StructLensPanelModel:
        return replace(self, status=status)

    def with_analysis(self, analysis: AnalysisResult) -> StructLensPanelModel:
        return replace(self, analysis=analysis, status="Analysis complete.")


def build_qt_panel(parent: object | None = None) -> object:
    """Build a minimal Qt panel when a supported Qt binding is present.

    Importing StructLens remains safe outside PyMOL; callers receive an
    actionable error if the host distribution does not provide Qt.
    """

    qt = _load_qt()
    if qt is None:
        raise RuntimeError("StructLens GUI requires the Qt binding shipped with PyMOL")
    QWidget, QVBoxLayout, QLabel, QTabWidget = qt
    panel = QWidget(parent)
    layout = QVBoxLayout(panel)
    title = QLabel("StructLens — reproducible protein-structure comparison")
    title.setObjectName("structlensTitle")
    layout.addWidget(title)
    tabs = QTabWidget(panel)
    for section in GUI_SECTIONS:
        page = QWidget(tabs)
        page_layout = QVBoxLayout(page)
        page_layout.addWidget(QLabel(_subtitle(section)))
        tabs.addTab(page, section)
    layout.addWidget(tabs)
    icon_path = Path(__file__).parents[1] / "assets" / "structlens_icon.png"
    if icon_path.exists() and hasattr(panel, "setWindowIcon"):
        QIcon = _load_qicon()
        if QIcon is not None:
            panel.setWindowIcon(QIcon(str(icon_path)))
    return panel


def _subtitle(section: str) -> str:
    return {
        "Project": "Choose the reference structure and the structures you want to compare.",
        "Alignment": "Choose how StructLens should determine equivalent residues before superposition.",
        "Mutations": "Review substitutions, insertions, deletions, and non-standard residues detected from the residue map.",
        "Residues": "Inspect residue-by-residue correspondence and structural differences.",
        "Visualization": "Choose what to highlight in PyMOL and how it should be represented.",
        "Results": "Review global alignment quality, structural metrics, and mutation counts.",
    }[section]


def _load_qt() -> tuple[Any, Any, Any, Any] | None:
    try:
        module = importlib.import_module("PySide6.QtWidgets")
        return module.QWidget, module.QVBoxLayout, module.QLabel, module.QTabWidget
    except ImportError:
        try:
            module = importlib.import_module("PyQt5.QtWidgets")
            return module.QWidget, module.QVBoxLayout, module.QLabel, module.QTabWidget
        except ImportError:
            return None


def _load_qicon() -> Any:
    try:
        return importlib.import_module("PySide6.QtGui").QIcon
    except ImportError:
        try:
            return importlib.import_module("PyQt5.QtGui").QIcon
        except ImportError:
            return None


__all__ = ["GUI_SECTIONS", "WORKFLOW_HELP", "StructLensPanelModel", "build_qt_panel"]
