"""Qt composition and controller for the StructLens PyMOL workflow."""

from __future__ import annotations

import re
from collections.abc import Mapping
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from tempfile import NamedTemporaryFile
from threading import Event
from typing import Any

from structlens.application.analysis_service import AnalysisService
from structlens.application.chart_data import structural_deviation_profile
from structlens.application.chart_export import export_chart_image, export_chart_xlsx
from structlens.application.export_service import (
    export_analysis_csv,
    export_analysis_json,
    export_analysis_xlsx,
)
from structlens.application.project_state import ProjectState
from structlens.core.errors import AnalysisCancelledError, BundleValidationError
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
from structlens.core.parsing import load_structure
from structlens.integrations.pymol.adapter import PyMOLAdapter
from structlens.integrations.pymol.launcher import PyMOLLauncher
from structlens.integrations.pymol.selections import selection_name
from structlens.integrations.pymol_bundle import write_pymol_bundle
from structlens.integrations.usalign.executable import bundled_executable
from structlens.plugin.visualization.legends import (
    backbone_rmsd_legend,
    displacement_legend,
)
from structlens.plugin.visualization.renderer import (
    ColorMode,
    HighlightFilter,
    Representation,
    VisualizationRenderer,
    VisualizationState,
)
from structlens.resources.backends import backend_versions

from .main_panel import SCIENTIFIC_SECTIONS, WORKFLOW_HELP, StructLensPanelModel
from .qt_compat import QtBindings

_RESIDUE_TOKEN = re.compile(r"^(?P<chain>[^:]+):(?P<number>-?\d+)(?P<insertion>[A-Za-z]?)$")
_STRUCTURES_PAGE_INDEX = SCIENTIFIC_SECTIONS.index("Structures")
_ANALYSIS_EXECUTOR = ThreadPoolExecutor(
    max_workers=1, thread_name_prefix="structlens-analysis"
)
_RESULTS_PAGE_INDEX = SCIENTIFIC_SECTIONS.index("Results")


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


class PanelController:
    """Own Qt widgets and side effects while keeping the model immutable."""

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
        self.reference_object_name: str | None = None
        self.target_object_name: str | None = None
        self._temporary_paths: list[Path] = []
        self._analysis_service = AnalysisService()
        self._renderer = VisualizationRenderer()
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
        self._build_shell()
        self._build_project_page()
        self._build_mutations_page()
        self._build_alignment_page()
        self._build_residues_page()
        self._build_sites_page()
        self._build_visualization_page()
        self._build_pymol_page()
        self._build_export_page()
        self._build_results_page()
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
    def _build_shell(self) -> None:
        self.widget.setStyleSheet(_stylesheet())
        root = self.w.QVBoxLayout(self.widget)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        header = self.w.QFrame(self.widget)
        header.setObjectName("header")
        header_layout = self.w.QHBoxLayout(header)
        header_layout.setContentsMargins(24, 18, 24, 18)
        header_layout.setSpacing(16)
        brand = self.w.QVBoxLayout()
        brand.setSpacing(2)
        brand.addWidget(_label(self.w, "STRUCTLENS / EVIDENCE BENCH", "eyebrow"))
        brand.addWidget(_label(self.w, "Integrated sequence and structure analysis", "windowTitle"))
        header_layout.addLayout(brand)
        header_layout.addStretch(1)
        self.header_status = _label(self.w, "Ready", "statusPill")
        header_layout.addWidget(self.header_status)
        self.cancel_button = _button(self.w, "Cancel", "secondaryButton")
        self.cancel_button.setVisible(False)
        self.cancel_button.clicked.connect(self._cancel_analysis)
        header_layout.addWidget(self.cancel_button)
        self.compare_button = _button(self.w, "Compare", "primaryButton")
        header_layout.addWidget(self.compare_button)
        root.addWidget(header)

        body = self.w.QWidget(self.widget)
        body_layout = self.w.QHBoxLayout(body)
        body_layout.setContentsMargins(0, 0, 0, 0)
        body_layout.setSpacing(0)
        self.sidebar = self.w.QFrame(body)
        self.sidebar.setObjectName("sidebar")
        self.sidebar.setMinimumWidth(200)
        self.sidebar.setMaximumWidth(250)
        side_layout = self.w.QVBoxLayout(self.sidebar)
        side_layout.setContentsMargins(14, 18, 14, 18)
        side_layout.setSpacing(12)
        side_layout.addWidget(_label(self.w, "WORKFLOW", "sectionKicker"))
        self.nav = self.w.QListWidget(self.sidebar)
        self.nav.setObjectName("workflowNav")
        self.nav.setSpacing(4)
        for section in SCIENTIFIC_SECTIONS:
            item = self.w.QListWidgetItem(section)
            item.setToolTip(_page_subtitle(section, standalone=self.command is None))
            self.nav.addItem(item)
        self.nav.setCurrentRow(0)
        side_layout.addWidget(self.nav, 1)
        side_layout.addWidget(
            _label(
                self.w,
                (
                    "The correspondence table is the source of truth. PyMOL only renders reversible selections."
                    if self.command is not None
                    else "The correspondence table is the source of truth. Filters and exports remain available here."
                ),
                "sidebarNote",
            )
        )
        body_layout.addWidget(self.sidebar)

        self.pages = self.w.QStackedWidget(body)
        self.pages.setObjectName("pageStack")
        body_layout.addWidget(self.pages, 1)
        root.addWidget(body, 1)

        footer = self.w.QFrame(self.widget)
        footer.setObjectName("footer")
        footer_layout = self.w.QHBoxLayout(footer)
        footer_layout.setContentsMargins(20, 8, 20, 8)
        self.footer_status = _label(self.w, "Choose two structures to begin.", "footerStatus")
        footer_layout.addWidget(self.footer_status, 1)
        self.progress = self.w.QProgressBar(footer)
        self.progress.setObjectName("analysisProgress")
        self.progress.setRange(0, 0)
        self.progress.setFixedWidth(180)
        self.progress.setVisible(False)
        footer_layout.addWidget(self.progress)
        footer_layout.addWidget(_label(self.w, "v0.3.0 · units are explicit · descriptive evidence only", "footerMeta"))
        root.addWidget(footer)

        self.compare_button.clicked.connect(self._start_analysis)

    def _add_page(self, title: str, purpose: str) -> tuple[Any, Any]:
        page = self.w.QWidget(self.pages)
        page.setObjectName(f"page{title}")
        outer = self.w.QVBoxLayout(page)
        outer.setContentsMargins(30, 28, 34, 32)
        outer.setSpacing(18)
        outer.addWidget(_label(self.w, title, "pageTitle"))
        outer.addWidget(_label(self.w, purpose, "pagePurpose"))
        content = self.w.QVBoxLayout()
        content.setSpacing(16)
        outer.addLayout(content)
        outer.addStretch(1)
        self.pages.addWidget(page)
        return page, content

    # --------------------------------------------------------------- Project
    def _build_project_page(self) -> None:
        _, content = self._add_page(
            "Project",
            "Load two coordinate sources, select chains, and keep the evidence trail reproducible.",
        )
        source_group = self.w.QGroupBox("Sources", self.widget)
        source_layout = self.w.QGridLayout(source_group)
        source_layout.setContentsMargins(18, 20, 18, 18)
        source_layout.setHorizontalSpacing(12)
        source_layout.setVerticalSpacing(12)
        source_layout.addWidget(_label(self.w, "REFERENCE", "fieldLabel"), 0, 0)
        self.reference_edit = self.w.QLineEdit(source_group)
        self.reference_edit.setPlaceholderText("PDB or mmCIF file")
        source_layout.addWidget(self.reference_edit, 0, 1)
        reference_browse = _button(self.w, "Browse…", "secondaryButton")
        reference_browse.clicked.connect(lambda: self._browse_source("reference"))
        source_layout.addWidget(reference_browse, 0, 2)
        self.reference_meta = _label(self.w, "Not loaded", "fieldMeta")
        source_layout.addWidget(self.reference_meta, 0, 3)
        reference_object = _button(self.w, "Use object…", "secondaryButton")
        reference_object.clicked.connect(lambda: self._use_pymol_object("reference"))
        reference_object.setVisible(self.command is not None)
        source_layout.addWidget(reference_object, 0, 4)
        source_layout.addWidget(_label(self.w, "TARGET", "fieldLabel"), 1, 0)
        self.target_edit = self.w.QLineEdit(source_group)
        self.target_edit.setPlaceholderText("PDB or mmCIF file")
        source_layout.addWidget(self.target_edit, 1, 1)
        target_browse = _button(self.w, "Browse…", "secondaryButton")
        target_browse.clicked.connect(lambda: self._browse_source("target"))
        source_layout.addWidget(target_browse, 1, 2)
        self.target_meta = _label(self.w, "Not loaded", "fieldMeta")
        source_layout.addWidget(self.target_meta, 1, 3)
        target_object = _button(self.w, "Use object…", "secondaryButton")
        target_object.clicked.connect(lambda: self._use_pymol_object("target"))
        target_object.setVisible(self.command is not None)
        source_layout.addWidget(target_object, 1, 4)
        source_layout.setColumnStretch(1, 1)
        source_layout.setColumnStretch(3, 1)
        if self.command is None:
            source_layout.addWidget(
                _label(
                    self.w,
                    "Standalone mode · load coordinate files here; PyMOL object sources are available in the plugin.",
                    "fieldMeta",
                ),
                2,
                0,
                1,
                5,
            )
        content.addWidget(source_group)

        chain_group = self.w.QGroupBox("Chain selection", self.widget)
        chain_layout = self.w.QGridLayout(chain_group)
        chain_layout.setContentsMargins(18, 20, 18, 18)
        chain_layout.setHorizontalSpacing(14)
        chain_layout.setVerticalSpacing(8)
        chain_layout.addWidget(_label(self.w, "REFERENCE CHAIN", "fieldLabel"), 0, 0)
        self.reference_chain_combo = self.w.QComboBox(chain_group)
        chain_layout.addWidget(self.reference_chain_combo, 1, 0)
        chain_layout.addWidget(_label(self.w, "TARGET CHAIN", "fieldLabel"), 0, 1)
        self.target_chain_combo = self.w.QComboBox(chain_group)
        chain_layout.addWidget(self.target_chain_combo, 1, 1)
        self.reference_chain_combo.currentIndexChanged.connect(self._chain_changed)
        self.target_chain_combo.currentIndexChanged.connect(self._chain_changed)
        chain_layout.setColumnStretch(0, 1)
        chain_layout.setColumnStretch(1, 1)
        content.addWidget(chain_group)

        note = _label(
            self.w,
            "StructLens reads coordinates without modifying source files. Save Project captures source paths, settings, analysis, and SHA-256 hashes.",
            "inlineNote",
        )
        note.setWordWrap(True)
        content.addWidget(note)
        backend_group = self.w.QGroupBox("About · Scientific Backends", self.widget)
        backend_layout = self.w.QVBoxLayout(backend_group)
        backend_layout.setContentsMargins(18, 14, 18, 14)
        versions = backend_versions()
        backend_layout.addWidget(_label(self.w, " · ".join(f"{key}: {value}" for key, value in versions.items()), "fieldMeta"))
        content.addWidget(backend_group)
        actions = self.w.QHBoxLayout()
        load_button = _button(self.w, "Load sources", "secondaryButton")
        load_button.clicked.connect(self._load_sources_from_edits)
        actions.addWidget(load_button)
        save_button = _button(self.w, "Save Project…", "secondaryButton")
        save_button.clicked.connect(self._save_project)
        actions.addWidget(save_button)
        open_button = _button(self.w, "Open Project…", "secondaryButton")
        open_button.clicked.connect(self._open_project)
        actions.addWidget(open_button)
        actions.addStretch(1)
        content.addLayout(actions)

    # -------------------------------------------------------------- Alignment
    def _build_alignment_page(self) -> None:
        _, content = self._add_page(
            "Structures",
            "Compare folds and choose the correspondence policy before superposition. The bundled US-align status is recorded in provenance.",
        )
        policy_group = self.w.QGroupBox("Correspondence policy", self.widget)
        policy_layout = self.w.QGridLayout(policy_group)
        policy_layout.setContentsMargins(18, 20, 18, 18)
        policy_layout.setHorizontalSpacing(14)
        policy_layout.setVerticalSpacing(10)
        policy_layout.addWidget(_label(self.w, "MODE", "fieldLabel"), 0, 0)
        self.mode_combo = self.w.QComboBox(policy_group)
        for label, value in (
            ("Auto", AlignmentMode.AUTO.value),
            ("Sequence", AlignmentMode.SEQUENCE.value),
            ("Structure · US-align", AlignmentMode.STRUCTURE.value),
            ("Manual · locked pairs", AlignmentMode.MANUAL.value),
        ):
            self.mode_combo.addItem(label, value)
        policy_layout.addWidget(self.mode_combo, 1, 0)
        self.mode_help = _label(self.w, "", "helpText")
        self.mode_help.setWordWrap(True)
        policy_layout.addWidget(self.mode_help, 1, 1, 1, 2)
        policy_layout.setColumnStretch(1, 1)
        self.mode_combo.currentTextChanged.connect(self._update_mode_help)
        content.addWidget(policy_group)

        comparison_group = self.w.QGroupBox("Comparison mode", self.widget)
        comparison_layout = self.w.QGridLayout(comparison_group)
        comparison_layout.setContentsMargins(18, 20, 18, 18)
        comparison_layout.setHorizontalSpacing(14)
        comparison_layout.setVerticalSpacing(8)
        comparison_layout.addWidget(_label(self.w, "TOPOLOGY", "fieldLabel"), 0, 0)
        self.comparison_combo = self.w.QComboBox(comparison_group)
        # The desktop panel currently has one reference and one target source
        # picker. Keep this selector honest until a multi-target source
        # collection workflow is promoted into the GUI; the application service
        # already exposes the v0.2 multi-structure APIs for scripted use.
        for label, value in (("Pairwise · one reference + one target", ComparisonMode.PAIRWISE.value),):
            self.comparison_combo.addItem(label, value)
        comparison_layout.addWidget(self.comparison_combo, 1, 0)
        self.comparison_help = _label(self.w, "", "helpText")
        self.comparison_help.setWordWrap(True)
        comparison_layout.addWidget(self.comparison_help, 1, 1)
        comparison_layout.setColumnStretch(1, 1)
        self.comparison_combo.currentTextChanged.connect(self._update_comparison_help)
        content.addWidget(comparison_group)

        thresholds = self.w.QGroupBox("Auto thresholds and refinement", self.widget)
        threshold_layout = self.w.QGridLayout(thresholds)
        threshold_layout.setContentsMargins(18, 20, 18, 18)
        threshold_layout.setHorizontalSpacing(14)
        threshold_layout.setVerticalSpacing(10)
        threshold_layout.addWidget(_label(self.w, "MINIMUM IDENTITY", "fieldLabel"), 0, 0)
        self.identity_spin = self.w.QDoubleSpinBox(thresholds)
        _fraction_spin(self.identity_spin, 0.30)
        threshold_layout.addWidget(self.identity_spin, 1, 0)
        threshold_layout.addWidget(_label(self.w, "MINIMUM COVERAGE", "fieldLabel"), 0, 1)
        self.coverage_spin = self.w.QDoubleSpinBox(thresholds)
        _fraction_spin(self.coverage_spin, 0.70)
        threshold_layout.addWidget(self.coverage_spin, 1, 1)
        threshold_layout.addWidget(_label(self.w, "US-ALIGN", "fieldLabel"), 0, 2)
        self.usalign_status = _label(self.w, "Bundled backend · Ready", "fieldMeta")
        threshold_layout.addWidget(self.usalign_status, 1, 2)
        self.usalign_edit = self.w.QLineEdit(thresholds)
        self.usalign_edit.setPlaceholderText("Optional custom executable (Advanced Settings)")
        threshold_layout.addWidget(_label(self.w, "CUSTOM EXECUTABLE (ADVANCED)", "fieldLabel"), 4, 0)
        threshold_layout.addWidget(self.usalign_edit, 5, 0, 1, 3)
        self.refined_check = self.w.QCheckBox("Refine outliers after strict fit", thresholds)
        threshold_layout.addWidget(self.refined_check, 2, 0, 1, 2)
        threshold_layout.addWidget(_label(self.w, "CUTOFF (Å)", "fieldLabel"), 2, 2)
        self.cutoff_spin = self.w.QDoubleSpinBox(thresholds)
        self.cutoff_spin.setRange(0.1, 50.0)
        self.cutoff_spin.setValue(2.0)
        self.cutoff_spin.setDecimals(1)
        self.cutoff_spin.setSuffix(" Å")
        threshold_layout.addWidget(self.cutoff_spin, 3, 2)
        threshold_layout.setColumnStretch(2, 1)
        content.addWidget(thresholds)

        self.manual_group = self.w.QGroupBox("Manual pairs", self.widget)
        manual_layout = self.w.QVBoxLayout(self.manual_group)
        manual_layout.setContentsMargins(18, 20, 18, 18)
        manual_layout.addWidget(
            _label(self.w, "One pair per line: reference_chain:auth_resi -> target_chain:auth_resi", "fieldMeta")
        )
        self.manual_edit = self.w.QPlainTextEdit(self.manual_group)
        self.manual_edit.setPlaceholderText("A:42 -> A:42\nA:43 -> A:43")
        self.manual_edit.setFixedHeight(86)
        manual_layout.addWidget(self.manual_edit)
        content.addWidget(self.manual_group)
        run_layout = self.w.QHBoxLayout()
        run_layout.addStretch(1)
        self.run_button = _button(self.w, "Run comparison", "primaryButton")
        self.run_button.clicked.connect(self._start_analysis)
        run_layout.addWidget(self.run_button)
        content.addLayout(run_layout)

    # -------------------------------------------------------------- Mutations
    def _build_mutations_page(self) -> None:
        _, content = self._add_page(
            "Sequences",
            "Compare amino-acid sequences, map equivalent positions, and inspect mutations and conservation.",
        )
        self.mutation_summary = _label(self.w, "No comparison yet.", "inlineNote")
        content.addWidget(self.mutation_summary)
        self.mutation_table = self.w.QTableWidget(0, 8, self.widget)
        self.mutation_table.setHorizontalHeaderLabels(
            ["Index", "Kind", "Reference", "Target", "Notation", "BLOSUM62", "Grantham", "Class"]
        )
        _configure_table(self.mutation_table, self.w)
        self.mutation_table.cellDoubleClicked.connect(self._focus_mutation)
        content.addWidget(self.mutation_table, 1)
        msa_group = self.w.QGroupBox("Multiple Sequence Alignment", self.widget)
        msa_layout = self.w.QVBoxLayout(msa_group)
        msa_layout.setContentsMargins(18, 16, 18, 16)
        self.msa_summary_label = _label(
            self.w,
            "No MSA result loaded. Run the MSA service and pass its immutable result to set_msa_result().",
            "inlineNote",
        )
        self.msa_summary_label.setWordWrap(True)
        msa_layout.addWidget(self.msa_summary_label)
        self.msa_table = self.w.QTableWidget(0, 3, msa_group)
        self.msa_table.setHorizontalHeaderLabels(["Structure", "Aligned sequence", "Source"])
        _configure_table(self.msa_table, self.w)
        self.msa_table.setMinimumHeight(150)
        msa_layout.addWidget(self.msa_table, 1)
        msa_layout.addWidget(
            _label(
                self.w,
                "Alignment columns retain reference-relative insertion labels; conservation excludes gaps and ambiguous residues from entropy.",
                "helpText",
            )
        )
        content.addWidget(msa_group, 1)

    def set_msa_result(self, alignment: MultipleSequenceAlignment | None) -> None:
        """Display an authoritative MSA without recalculating its values."""

        self.msa_table.setRowCount(0)
        if alignment is None:
            self.msa_summary_label.setText("MSA unavailable.")
            return
        self.msa_table.setRowCount(len(alignment.aligned_rows))
        for row_index, (structure_id, row) in enumerate(alignment.aligned_rows):
            source = next((sequence.source for sequence in alignment.sequences if sequence.structure_id == structure_id), "unknown")
            self.msa_table.setItem(row_index, 0, self.w.QTableWidgetItem(structure_id))
            self.msa_table.setItem(row_index, 1, self.w.QTableWidgetItem(row))
            self.msa_table.setItem(row_index, 2, self.w.QTableWidgetItem(source))
        self.msa_summary_label.setText(
            f"{len(alignment.aligned_rows)} sequences · {len(alignment.columns)} alignment columns · "
            f"{sum(column.reference_residue is None for column in alignment.columns)} reference-relative insertion columns."
        )

    # --------------------------------------------------------------- Residues
    def _build_residues_page(self) -> None:
        _, content = self._add_page(
            "Residues",
            (
                "Inspect the authoritative correspondence table; double-click a row to select it."
                if self.command is None
                else "Inspect the authoritative correspondence table; double-click a row to focus it in PyMOL."
            ),
        )
        self.residue_summary = _label(self.w, "No comparison yet.", "inlineNote")
        content.addWidget(self.residue_summary)
        self.residue_table = self.w.QTableWidget(0, 10, self.widget)
        self.residue_table.setHorizontalHeaderLabels(
            [
                "Index",
                "Reference",
                "Target",
                "Status",
                "Cα Δ (Å)",
                "Backbone (Å)",
                "Side-chain (Å)",
                "Heavy (Å)",
                "Outlier",
                "Key",
            ]
        )
        _configure_table(self.residue_table, self.w)
        self.residue_table.cellDoubleClicked.connect(self._focus_residue)
        content.addWidget(self.residue_table, 1)
        evidence_group = self.w.QGroupBox("Residue Evidence Card", self.widget)
        evidence_layout = self.w.QVBoxLayout(evidence_group)
        evidence_layout.setContentsMargins(18, 16, 18, 16)
        self.evidence_card_label = _label(
            self.w,
            "Select a correspondence row to inspect sequence, structure, interaction, site, and quality evidence. Missing values remain unavailable.",
            "helpText",
        )
        self.evidence_card_label.setWordWrap(True)
        evidence_layout.addWidget(self.evidence_card_label)
        content.addWidget(evidence_group)

    # ---------------------------------------------------------- Visualization
    def _build_sites_page(self) -> None:
        _, content = self._add_page(
            "Sites",
            "Define active sites or ligand-centered regions and compare coverage, geometry, exposure, and interaction fingerprints.",
        )
        site_group = self.w.QGroupBox("Site definition", self.widget)
        site_layout = self.w.QGridLayout(site_group)
        site_layout.setContentsMargins(18, 16, 18, 16)
        site_layout.setHorizontalSpacing(12)
        site_layout.setVerticalSpacing(8)
        self.site_mode_combo = self.w.QComboBox(site_group)
        for label, value in (
            ("Key residues", "key_residues"),
            ("Ligand radius", "ligand_radius"),
            ("Residue radius", "residue_radius"),
        ):
            self.site_mode_combo.addItem(label, value)
        site_layout.addWidget(_label(self.w, "DEFINITION", "fieldLabel"), 0, 0)
        site_layout.addWidget(self.site_mode_combo, 1, 0)
        self.site_residues_edit = self.w.QLineEdit(site_group)
        self.site_residues_edit.setPlaceholderText("A:70, A:73, A:166")
        site_layout.addWidget(_label(self.w, "REFERENCE POSITIONS", "fieldLabel"), 0, 1)
        site_layout.addWidget(self.site_residues_edit, 1, 1)
        self.site_radius_spin = self.w.QDoubleSpinBox(site_group)
        self.site_radius_spin.setRange(0.1, 20.0)
        self.site_radius_spin.setValue(5.0)
        self.site_radius_spin.setSuffix(" Å")
        site_layout.addWidget(_label(self.w, "RADIUS", "fieldLabel"), 0, 2)
        site_layout.addWidget(self.site_radius_spin, 1, 2)
        self.site_define_button = _button(self.w, "Define site", "secondaryButton")
        self.site_define_button.clicked.connect(self._define_site_from_controls)
        site_layout.addWidget(self.site_define_button, 1, 3)
        self.site_status_label = _label(
            self.w,
            "Site metrics appear after an analysis service run; this panel never estimates them locally.",
            "inlineNote",
        )
        self.site_status_label.setWordWrap(True)
        site_layout.addWidget(self.site_status_label, 2, 0, 1, 4)
        content.addWidget(site_group)

        metrics = self.w.QGroupBox("Authoritative site metrics", self.widget)
        metrics_layout = self.w.QVBoxLayout(metrics)
        metrics_layout.setContentsMargins(18, 16, 18, 16)
        metrics_layout.addWidget(
            _label(
                self.w,
                "Coverage, global-frame/site-fitted RMSD, SASA (Å²), atomic envelope volume (Å³), and interaction fingerprints are read from the site service result.",
                "helpText",
            )
        )
        content.addWidget(metrics)

    # ---------------------------------------------------------- Visualization
    def _build_visualization_page(self) -> None:
        _, content = self._add_page(
            "Charts",
            "Use authoritative sequence and structure datasets, inspect what each chart means, and route exports through a reproducible workflow.",
        )
        chart_group = self.w.QGroupBox("Scientific chart", self.widget)
        chart_layout = self.w.QGridLayout(chart_group)
        chart_layout.setContentsMargins(18, 20, 18, 18)
        chart_layout.setHorizontalSpacing(14)
        chart_layout.setVerticalSpacing(10)
        chart_layout.addWidget(_label(self.w, "PROFILE", "fieldLabel"), 0, 0)
        self.chart_combo = self.w.QComboBox(chart_group)
        for label in (
            "Structural deviation profile",
            "Mutation / conservation matrix",
            "Pairwise similarity heatmap",
            "Sequence–structure relationship",
            "Structural conservation profile",
            "Key-residue comparison",
            "MSA conservation profile",
            "Sequence logo",
            "Interaction difference matrix",
            "Site comparison",
            "Distance-difference heatmap",
            "Evidence Card figure",
        ):
            self.chart_combo.addItem(label)
        chart_layout.addWidget(self.chart_combo, 1, 0)
        self.chart_explanation = _label(
            self.w,
            "Charts consume the authoritative analysis state. Values include units and remain exportable as data.",
            "helpText",
        )
        self.chart_explanation.setWordWrap(True)
        chart_layout.addWidget(self.chart_explanation, 1, 1)
        chart_layout.setColumnStretch(1, 1)
        self.chart_combo.currentTextChanged.connect(self._update_chart_explanation)
        content.addWidget(chart_group)
        note = _label(
            self.w,
            "Use the Export page for XLSX, JPEG, or TIFF output. PyMOL-specific filtering and reversible rendering controls remain on the PyMOL page.",
            "inlineNote",
        )
        note.setWordWrap(True)
        content.addWidget(note)

    # --------------------------------------------------------------- PyMOL
    def _build_pymol_page(self) -> None:
        _, content = self._add_page(
            "PyMOL",
            "Prepare the current StructLens analysis for interactive 3D exploration in PyMOL.",
        )
        controls = self.w.QGroupBox("Visualization controls", self.widget)
        layout = self.w.QGridLayout(controls)
        layout.setContentsMargins(18, 20, 18, 18)
        layout.setHorizontalSpacing(14)
        layout.setVerticalSpacing(10)
        self.preset_combo = self.w.QComboBox(controls)
        for name in ("Minimal", "Publication", "Mutation focus", "Structural deviation", "Active site", "Presentation"):
            self.preset_combo.addItem(name)
        self.preset_combo.currentTextChanged.connect(self._apply_preset)
        layout.addWidget(_label(self.w, "PRESET", "fieldLabel"), 0, 0)
        layout.addWidget(self.preset_combo, 1, 0)
        self.filter_combo = self.w.QComboBox(controls)
        for enum in HighlightFilter:
            self.filter_combo.addItem(_human(enum.value), enum.value)
        layout.addWidget(_label(self.w, "HIGHLIGHT", "fieldLabel"), 0, 1)
        layout.addWidget(self.filter_combo, 1, 1)
        self.color_combo = self.w.QComboBox(controls)
        for color_mode in ColorMode:
            self.color_combo.addItem(_human(color_mode.value), color_mode.value)
        layout.addWidget(_label(self.w, "COLOR BY", "fieldLabel"), 0, 2)
        layout.addWidget(self.color_combo, 1, 2)
        self.representation_combo = self.w.QComboBox(controls)
        for representation in Representation:
            self.representation_combo.addItem(_human(representation.value), representation.value)
        layout.addWidget(_label(self.w, "REPRESENTATION", "fieldLabel"), 2, 0)
        layout.addWidget(self.representation_combo, 3, 0)
        self.radius_spin = self.w.QDoubleSpinBox(controls)
        self.radius_spin.setRange(0.5, 30.0)
        self.radius_spin.setValue(5.0)
        self.radius_spin.setSuffix(" Å")
        self.radius_spin.setDecimals(1)
        layout.addWidget(_label(self.w, "LOCAL RADIUS", "fieldLabel"), 2, 1)
        layout.addWidget(self.radius_spin, 3, 1)
        self.labels_check = self.w.QCheckBox("Show residue labels", controls)
        self.reference_check = self.w.QCheckBox("Reference", controls)
        self.reference_check.setChecked(True)
        self.target_check = self.w.QCheckBox("Target", controls)
        self.target_check.setChecked(True)
        layout.addWidget(self.labels_check, 2, 2)
        visible_layout = self.w.QHBoxLayout()
        visible_layout.addWidget(self.reference_check)
        visible_layout.addWidget(self.target_check)
        visible_layout.addStretch(1)
        layout.addLayout(visible_layout, 3, 2)
        content.addWidget(controls)
        for control in (
            self.filter_combo,
            self.color_combo,
            self.representation_combo,
            self.radius_spin,
            self.labels_check,
            self.reference_check,
            self.target_check,
        ):
            signal = getattr(control, "currentIndexChanged", None) or getattr(control, "valueChanged", None) or getattr(control, "toggled", None)
            if signal is not None:
                signal.connect(self._visualization_changed)
        self.legend_label = _label(self.w, "Legend appears after a comparison.", "legend")
        self.legend_label.setWordWrap(True)
        content.addWidget(self.legend_label)
        self.visualization_count = _label(self.w, "0 rows selected", "inlineNote")
        content.addWidget(self.visualization_count)
        integration = self.w.QGroupBox("PyMOL integration", self.widget)
        integration_layout = self.w.QGridLayout(integration)
        integration_layout.setContentsMargins(18, 20, 18, 18)
        integration_layout.setHorizontalSpacing(14)
        integration_layout.setVerticalSpacing(10)
        integration_layout.addWidget(_label(self.w, "STATUS", "fieldLabel"), 0, 0)
        self.pymol_status = _label(self.w, "Ready to export", "helpText")
        self.pymol_status.setWordWrap(True)
        integration_layout.addWidget(self.pymol_status, 1, 0, 1, 3)
        integration_layout.addWidget(_label(self.w, "EXECUTABLE (OPTIONAL)", "fieldLabel"), 2, 0)
        self.pymol_edit = self.w.QLineEdit(integration)
        self.pymol_edit.setPlaceholderText("Configured PyMOL executable or PATH-resolved command")
        self.pymol_edit.textChanged.connect(lambda _text: self._refresh_pymol_status())
        integration_layout.addWidget(self.pymol_edit, 3, 0, 1, 3)
        actions = self.w.QHBoxLayout()
        open_button = _button(self.w, "Open in PyMOL", "primaryButton")
        open_button.clicked.connect(self._open_in_pymol)
        actions.addWidget(open_button)
        export_bundle = _button(self.w, "Export for PyMOL…", "secondaryButton")
        export_bundle.clicked.connect(self._export_pymol_bundle)
        actions.addWidget(export_bundle)
        plugin_help = _button(self.w, "Plugin installation instructions", "secondaryButton")
        plugin_help.clicked.connect(
            lambda: self._set_status(
                "Install StructLens-PyMOL from the amgoncalvesusp/pymol-plugins GitHub release, then open the bundle in PyMOL."
            )
        )
        actions.addWidget(plugin_help)
        actions.addStretch(1)
        integration_layout.addLayout(actions, 4, 0, 1, 3)
        content.addWidget(integration)
        host_actions = self.w.QHBoxLayout()
        host_actions.addStretch(1)
        if self.command is None:
            host_actions.addWidget(
                _label(
                    self.w,
                    "Standalone mode · Open in PyMOL creates a validated bundle and launches the external application when configured.",
                    "fieldMeta",
                )
            )
        else:
            reset = _button(self.w, "Reset StructLens view", "secondaryButton")
            reset.clicked.connect(self._reset_visualization)
            host_actions.addWidget(reset)
            apply_button = _button(self.w, "Apply to PyMOL", "primaryButton")
            apply_button.clicked.connect(self._apply_visualization)
            host_actions.addWidget(apply_button)
        content.addLayout(host_actions)

    # --------------------------------------------------------------- Export
    def _build_export_page(self) -> None:
        _, content = self._add_page(
            "Export",
            "Write the current evidence state as tabular data, chart data, or a portable PyMOL interchange bundle.",
        )
        exports = self.w.QGroupBox("Evidence exports", self.widget)
        export_layout = self.w.QHBoxLayout(exports)
        export_layout.setContentsMargins(18, 20, 18, 18)
        for label, callback in (
            ("XLSX", self._export_xlsx),
            ("CSV", self._export_csv),
            ("JSON", self._export_json),
            ("Chart XLSX", self._export_chart_xlsx),
            ("Chart JPEG", lambda: self._export_chart_image("jpeg", 300)),
            ("Chart TIFF", lambda: self._export_chart_image("tiff", 600)),
        ):
            button = _button(self.w, f"Export {label}…", "secondaryButton")
            button.clicked.connect(callback)
            if label.startswith("Chart"):
                self.chart_export_buttons.append(button)
            export_layout.addWidget(button)
        export_layout.addStretch(1)
        content.addWidget(exports)
        self._update_chart_export_state(self.chart_combo.currentText())
        note = _label(
            self.w,
            "The selected chart profile controls chart exports. Non-deviation profiles remain available as structured data through the application API.",
            "inlineNote",
        )
        note.setWordWrap(True)
        content.addWidget(note)

    # ---------------------------------------------------------------- Results
    def _build_results_page(self) -> None:
        _, content = self._add_page(
            "Results",
            (
                "Review global metrics, branch choice, and export the same result used by the evidence tables."
                if self.command is None
                else "Review global metrics, branch choice, and export the same result used by the PyMOL view."
            ),
        )
        self.result_decision = _label(self.w, "No comparison yet.", "resultDecision")
        self.result_decision.setWordWrap(True)
        content.addWidget(self.result_decision)
        metrics = self.w.QGroupBox("Global metrics", self.widget)
        metric_layout = self.w.QGridLayout(metrics)
        metric_layout.setContentsMargins(18, 20, 18, 18)
        metric_layout.setHorizontalSpacing(32)
        metric_layout.setVerticalSpacing(14)
        self.result_labels: dict[str, Any] = {}
        for row, (key, title, unit) in enumerate(
            (
                ("sequence_identity", "Sequence identity", "fraction"),
                ("sequence_similarity", "Sequence similarity", "fraction"),
                ("sequence_coverage", "Sequence coverage", "fraction"),
                ("strict_rmsd_angstrom", "Strict Cα RMSD", "Å"),
                ("refined_rmsd_angstrom", "Refined Cα RMSD", "Å"),
                ("tm_score", "TM-score", "score"),
                ("mapped_residue_count", "Mapped residues", "residues"),
                ("mutation_count", "Mutation count", "events"),
            )
        ):
            metric_layout.addWidget(_label(self.w, title, "fieldLabel"), row, 0)
            value = _label(self.w, "—", "metricValue")
            metric_layout.addWidget(value, row, 1)
            metric_layout.addWidget(_label(self.w, unit, "fieldMeta"), row, 2)
            self.result_labels[key] = value
        metric_layout.setColumnStretch(1, 1)
        content.addWidget(metrics)

    # --------------------------------------------------------------- bindings
    def _wire_navigation(self) -> None:
        self.nav.currentRowChanged.connect(self.pages.setCurrentIndex)

    def _browse_source(self, role: str) -> None:
        path, _ = self.w.QFileDialog.getOpenFileName(
            self.widget,
            f"Choose {role} structure",
            "",
            "Structure files (*.pdb *.ent *.cif *.mmcif *.gz);;All files (*)",
        )
        if path:
            edit = self.reference_edit if role == "reference" else self.target_edit
            edit.setText(path)
            self._load_source(role, Path(path))

    def _load_sources_from_edits(self) -> bool:
        loaded = True
        if self.reference_edit.text().strip():
            loaded = self._load_source("reference", Path(self.reference_edit.text().strip())) and loaded
        else:
            self._clear_source("reference")
            loaded = False
        if self.target_edit.text().strip():
            loaded = self._load_source("target", Path(self.target_edit.text().strip())) and loaded
        else:
            self._clear_source("target")
            loaded = False
        return loaded

    def _load_source(
        self,
        role: str,
        path: Path,
        *,
        pymol_object_name: str | None = None,
    ) -> bool:
        self._clear_source(role)
        try:
            structure = load_structure(path)
            if role == "reference":
                self.reference_structure = structure
                self.reference_object_name = pymol_object_name or self._load_file_into_pymol(role, path)
                self.reference_edit.setText(str(path))
                self._populate_chains(self.reference_chain_combo, structure)
                self.reference_meta.setText(_structure_meta(structure))
            else:
                self.target_structure = structure
                self.target_object_name = pymol_object_name or self._load_file_into_pymol(role, path)
                self.target_edit.setText(str(path))
                self._populate_chains(self.target_chain_combo, structure)
                self.target_meta.setText(_structure_meta(structure))
            self._sync_source_model()
            self._set_status(f"{role.title()} loaded · choose chains or run comparison")
            return True
        except Exception as exc:  # parser errors are user-facing recovery states
            self._show_error(f"Could not load {role}: {exc}")
            return False

    def _clear_source(self, role: str) -> None:
        if role == "reference":
            self.reference_structure = None
            self.reference_object_name = None
            self.reference_chain_combo.clear()
            self.reference_meta.setText("Not loaded")
        else:
            self.target_structure = None
            self.target_object_name = None
            self.target_chain_combo.clear()
            self.target_meta.setText("Not loaded")
        self._sync_source_model()

    def _use_pymol_object(self, role: str) -> None:
        if self.command is None:
            self._show_error("PyMOL object sources are available only inside a PyMOL session.")
            return
        get_names = getattr(self.command, "get_names", None)
        save = getattr(self.command, "save", None)
        if get_names is None or save is None:
            self._show_error("The PyMOL command proxy cannot list and save objects.")
            return
        names = [str(name) for name in get_names("objects")]
        if not names:
            self._show_error("No PyMOL objects are available. Load a structure first.")
            return
        selected, accepted = self.w.QInputDialog.getItem(
            self.widget,
            f"Choose {role} object",
            "PyMOL object",
            names,
            0,
            False,
        )
        if not accepted or not selected:
            return
        with NamedTemporaryFile(prefix="structlens_", suffix=".pdb", delete=False) as handle:
            temporary = Path(handle.name)
        try:
            save(str(temporary), selected)
            self._temporary_paths.append(temporary)
            self._load_source(role, temporary, pymol_object_name=str(selected))
        except Exception as exc:
            temporary.unlink(missing_ok=True)
            self._show_error(f"Could not read PyMOL object {selected}: {exc}")

    def _load_named_pymol_object(self, role: str, object_name: str) -> bool:
        save = getattr(self.command, "save", None)
        if save is None:
            return False
        with NamedTemporaryFile(prefix="structlens_", suffix=".pdb", delete=False) as handle:
            temporary = Path(handle.name)
        try:
            save(str(temporary), object_name)
            self._temporary_paths.append(temporary)
            return self._load_source(
                role, temporary, pymol_object_name=object_name
            )
        except Exception:
            temporary.unlink(missing_ok=True)
            return False

    def _load_file_into_pymol(self, role: str, path: Path) -> str | None:
        if self.command is None:
            return None
        load = getattr(self.command, "load", None)
        if load is None:
            return None
        object_name = selection_name("panel", path.stem, role)
        try:
            get_names = getattr(self.command, "get_names", None)
            names = {str(name) for name in get_names("objects")} if get_names is not None else set()
            if object_name not in names:
                load(str(path), object_name)
            return object_name
        except Exception:
            return None

    def _populate_chains(self, combo: Any, structure: ProteinStructure) -> None:
        combo.blockSignals(True)
        combo.clear()
        for chain in structure.chains:
            count = len(chain.residue_records or chain.residues)
            combo.addItem(f"Chain {chain.chain_id} · {count} residues", chain.chain_id)
        combo.blockSignals(False)
        if combo.count():
            combo.setCurrentIndex(0)

    def _chain_changed(self, _: int) -> None:
        self._sync_source_model()

    def _sync_source_model(self) -> None:
        self.model = self.model.with_sources(
            reference_path=self.reference_edit.text().strip() or None,
            target_path=self.target_edit.text().strip() or None,
            reference_chain_id=self._combo_data(self.reference_chain_combo),
            target_chain_id=self._combo_data(self.target_chain_combo),
        )

    def _combo_data(self, combo: Any) -> str | None:
        value = combo.currentData()
        return str(value) if value is not None else None

    # --------------------------------------------------------------- analysis
    def _settings(self) -> AnalysisSettings:
        mode = str(self.mode_combo.currentData() or AlignmentMode.AUTO.value)
        return AnalysisSettings(
            alignment_mode=AlignmentMode(mode),
            minimum_sequence_identity=float(self.identity_spin.value()),
            minimum_sequence_coverage=float(self.coverage_spin.value()),
            refined_rmsd=bool(self.refined_check.isChecked()),
            refinement_cutoff_angstrom=float(self.cutoff_spin.value()),
            usalign_executable=self.usalign_edit.text().strip() or None,
        )

    def _selected_chain(self, structure: ProteinStructure | None, combo: Any) -> ProteinChain | None:
        if structure is None:
            return None
        chain_id = self._combo_data(combo)
        for chain in structure.chains:
            if chain_id is None or chain.chain_id == chain_id:
                return chain
        return None

    def _start_analysis(self) -> None:
        if not self._load_sources_from_edits():
            self.nav.setCurrentRow(0)
            return
        reference = self._selected_chain(self.reference_structure, self.reference_chain_combo)
        target = self._selected_chain(self.target_structure, self.target_chain_combo)
        if reference is None or target is None:
            self._show_error("Load one reference and one target chain before comparing.")
            self.nav.setCurrentRow(0)
            return
        reference_chain = reference
        target_chain = target
        try:
            settings = self._settings()
            manual_pairs = self._manual_pairs(reference_chain, target_chain) if settings.alignment_mode is AlignmentMode.MANUAL else None
            service_manual_pairs: list[tuple[object, object]] | None = (
                [(reference_id, target_id) for reference_id, target_id in manual_pairs]
                if manual_pairs is not None
                else None
            )
        except ValueError as exc:
            self._show_error(str(exc))
            self.nav.setCurrentRow(_STRUCTURES_PAGE_INDEX)
            return
        self.model = self.model.with_busy("Comparing structures…")
        self._set_busy(True)
        self._cancel_event = Event()

        self._future = _ANALYSIS_EXECUTOR.submit(
            self._analysis_service.analyze,
            reference_chain,
            target_chain,
            settings,
            manual_pairs=service_manual_pairs,
            cancel_event=self._cancel_event,
        )
        self._poll_timer = self.c.QTimer(self.widget)
        self._poll_timer.setInterval(50)
        self._poll_timer.timeout.connect(self._poll_analysis)
        self._poll_timer.start()

    def _poll_analysis(self) -> None:
        if self._future is None or not self._future.done():
            return
        future = self._future
        self._future = None
        if self._poll_timer is not None:
            self._poll_timer.stop()
        self._poll_timer = None
        try:
            self._analysis_finished(future.result())
        except AnalysisCancelledError:
            self._analysis_cancelled()
        except Exception as exc:
            self._analysis_failed(str(exc))

    def _manual_pairs(self, reference: ProteinChain, target: ProteinChain) -> list[tuple[ResidueId, ResidueId]]:
        pairs: list[tuple[ResidueId, ResidueId]] = []
        for line_number, raw_line in enumerate(self.manual_edit.toPlainText().splitlines(), 1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            pieces = re.split(r"\s*(?:->|=|,)\s*", line, maxsplit=1)
            if len(pieces) != 2:
                raise ValueError(f"Manual pair line {line_number} must use ref -> target syntax")
            reference_id = _find_residue(reference, pieces[0].strip())
            target_id = _find_residue(target, pieces[1].strip())
            if reference_id is None or target_id is None:
                raise ValueError(f"Manual pair line {line_number} names a residue not in the selected chains")
            pairs.append((reference_id, target_id))
        if not pairs:
            raise ValueError("Manual mode needs at least one residue pair")
        return pairs

    def _analysis_finished(self, result: AnalysisResult) -> None:
        self.model = self.model.with_analysis(result)
        self._set_busy(False)
        self._populate_result(result)
        self._set_status(f"Analysis complete · {len(result.correspondences)} aligned positions")
        self.nav.setCurrentRow(_RESULTS_PAGE_INDEX)

    def _analysis_failed(self, message: str) -> None:
        self.model = self.model.with_error(f"Comparison failed: {message}")
        self._set_busy(False)
        self._show_error(f"Comparison failed: {message}")

    def _analysis_cancelled(self) -> None:
        self.model = self.model.with_status("Comparison cancelled.")
        self._set_busy(False)
        self._set_status("Comparison cancelled; no result was changed")

    def _set_busy(self, busy: bool) -> None:
        self.progress.setVisible(busy)
        self.cancel_button.setVisible(busy)
        self.compare_button.setEnabled(not busy)
        self.run_button.setEnabled(not busy)
        self.compare_button.setText("Comparing…" if busy else "Compare")
        self.run_button.setText("Comparing…" if busy else "Run comparison")

    def _cancel_analysis(self) -> None:
        self._cancel_event.set()
        if self._future is not None:
            self._future.cancel()
            self._future = None
        if self._poll_timer is not None:
            self._poll_timer.stop()
            self._poll_timer = None
        self._analysis_cancelled()

    # -------------------------------------------------------------- rendering
    def _state_from_controls(self) -> VisualizationState:
        return VisualizationState(
            highlight_filter=HighlightFilter(str(self.filter_combo.currentData())),
            color_mode=ColorMode(str(self.color_combo.currentData())),
            representation=Representation(str(self.representation_combo.currentData())),
            show_labels=self.labels_check.isChecked(),
            show_reference=self.reference_check.isChecked(),
            show_target=self.target_check.isChecked(),
            local_radius_angstrom=float(self.radius_spin.value()),
            preset=self.preset_combo.currentText(),
        )

    def _visualization_changed(self, *_: object) -> None:
        state = self._state_from_controls()
        self.model = self.model.with_visualization(state)
        self._update_legend()

    def _apply_preset(self, name: str) -> None:
        if not name:
            return
        state = self._renderer.apply_preset(name)
        self._set_combo_value(self.filter_combo, state.highlight_filter.value)
        self._set_combo_value(self.color_combo, state.color_mode.value)
        self._set_combo_value(self.representation_combo, state.representation.value)
        self.labels_check.setChecked(state.show_labels)
        self.reference_check.setChecked(state.show_reference)
        self.target_check.setChecked(state.show_target)
        self.model = self.model.with_visualization(state)
        self._update_legend()

    def _apply_visualization(self) -> None:
        result = self.model.analysis
        if result is None:
            self._show_error("Run a comparison before applying a PyMOL view.")
            return
        if self.command is None:
            self._show_error("Apply a visualization from inside a PyMOL session.")
            return
        state = self._state_from_controls()
        if self.target_object_name is None or (
            state.show_reference and self.reference_object_name is None
        ):
            self._show_error(
                "The selected sources are not loaded as PyMOL objects. Use Browse in PyMOL or choose Use object… first."
            )
            return
        self.model = self.model.with_visualization(state)
        self._pymol.apply(
            result,
            state=state,
            reference_object=self.reference_object_name,
            target_object=self.target_object_name,
        )
        selected = len(self._renderer.filtered_correspondences(result.correspondences, state))
        self._set_status(f"PyMOL view applied · {selected} correspondence rows selected")

    def _reset_visualization(self) -> None:
        self._pymol.reset()
        self._set_status("StructLens-owned PyMOL selections removed")

    def _update_legend(self) -> None:
        result = self.model.analysis
        if result is None:
            self.legend_label.setText("Legend appears after a comparison.")
            self.visualization_count.setText("0 rows selected")
            return
        state = self._state_from_controls()
        selected = self._renderer.filtered_correspondences(result.correspondences, state)
        self.visualization_count.setText(f"{len(selected)} of {len(result.correspondences)} rows selected")
        if state.color_mode is ColorMode.CA_DISPLACEMENT:
            values = [item.ca_displacement_angstrom for item in result.correspondences if item.ca_displacement_angstrom is not None]
            if values:
                legend = displacement_legend(min(values), max(values))
                self.legend_label.setText(f"Legend · {legend.title} · {legend.unit} · {legend.minimum:.2f}–{legend.maximum:.2f}")
                return
        if state.color_mode is ColorMode.BACKBONE_RMSD:
            values = [item.backbone_rmsd_angstrom for item in result.correspondences if item.backbone_rmsd_angstrom is not None]
            if values:
                legend = backbone_rmsd_legend(min(values), max(values))
                self.legend_label.setText(f"Legend · {legend.title} · {legend.unit} · {legend.minimum:.2f}–{legend.maximum:.2f}")
                return
        self.legend_label.setText("Legend · status colors are paired with text labels; no scientific meaning is inferred from color alone.")

    # -------------------------------------------------------------- result UI
    def _populate_result(self, result: AnalysisResult) -> None:
        self.result_decision.setText(f"<b>{result.alignment_decision}</b><br>Reference: {result.reference_id} · Target: {result.target_id}")
        values = {
            "sequence_identity": f"{result.sequence_identity:.3f}",
            "sequence_similarity": _number(result.sequence_similarity),
            "sequence_coverage": f"{result.sequence_coverage:.3f}",
            "strict_rmsd_angstrom": _number(result.strict_rmsd_angstrom),
            "refined_rmsd_angstrom": _number(result.refined_rmsd_angstrom),
            "tm_score": _number(result.tm_score),
            "mapped_residue_count": str(result.mapped_residue_count),
            "mutation_count": str(result.mutation_count),
        }
        for key, value in values.items():
            self.result_labels[key].setText(value)
        focus_hint = "double-click a row to focus it" if self.command is not None else "double-click a row to select it"
        self.mutation_summary.setText(f"{result.mutation_count} mutation event(s) · {focus_hint}")
        self.residue_summary.setText(f"{len(result.correspondences)} aligned positions · {result.mapped_residue_count} mapped Cα pairs")
        self._fill_mutations(result)
        self._fill_residues(result)
        self._update_legend()

    def _fill_mutations(self, result: AnalysisResult) -> None:
        self.mutation_table.setRowCount(0)
        for event in result.mutations:
            row = self.mutation_table.rowCount()
            self.mutation_table.insertRow(row)
            values = (
                str(event.alignment_index),
                event.kind.value,
                _residue_label(event.reference),
                _residue_label(event.target),
                event.canonical_notation,
                _number(event.blosum62_score),
                _number(event.grantham_distance),
                event.physicochemical_class or "—",
            )
            for column, value in enumerate(values):
                self.mutation_table.setItem(row, column, self.w.QTableWidgetItem(value))

    def _fill_residues(self, result: AnalysisResult) -> None:
        self.residue_table.setRowCount(0)
        for item in result.correspondences:
            row = self.residue_table.rowCount()
            self.residue_table.insertRow(row)
            values = (
                str(item.alignment_index),
                _residue_label(item.reference),
                _residue_label(item.target),
                item.status.value,
                _number(item.ca_displacement_angstrom),
                _number(item.backbone_rmsd_angstrom),
                _number(item.sidechain_rmsd_angstrom),
                _number(item.all_heavy_atom_rmsd_angstrom),
                "yes" if item.is_outlier else "no",
                "yes" if item.is_key_residue else "no",
            )
            for column, value in enumerate(values):
                self.residue_table.setItem(row, column, self.w.QTableWidgetItem(value))

    def _focus_mutation(self, row: int, _: int) -> None:
        result = self.model.analysis
        if result is None or row >= len(result.mutations):
            return
        alignment_index = result.mutations[row].alignment_index
        self._focus_alignment_index(result, alignment_index)

    def _focus_residue(self, row: int, _: int) -> None:
        result = self.model.analysis
        if result is None or row >= len(result.correspondences):
            return
        self._focus_alignment_index(result, result.correspondences[row].alignment_index)

    def _focus_alignment_index(self, result: AnalysisResult, alignment_index: int) -> None:
        item = next((entry for entry in result.correspondences if entry.alignment_index == alignment_index), None)
        if item is None:
            return
        self.evidence_card_label.setText(
            f"Reference {_residue_label(item.reference) or 'unavailable'} → Target {_residue_label(item.target) or 'unavailable'}\n"
            f"Status: {item.status.value} · Cα displacement: {_number(item.ca_displacement_angstrom)} Å · "
            f"Backbone RMSD: {_number(item.backbone_rmsd_angstrom)} Å · "
            "Interactions/site evidence: unavailable until the corresponding v0.3 service result is present."
        )
        selection = self._pymol.focus_residue(item, result.target_id)
        if selection:
            self._set_status(f"Focused {selection} in PyMOL")
        elif self.command is None:
            self._set_status(f"Residue {alignment_index} selected in the evidence table")
        else:
            self._set_status(f"Residue {alignment_index} selected; open inside PyMOL to focus it")

    # -------------------------------------------------------------- persistence
    def _save_project(self) -> None:
        path, _ = self.w.QFileDialog.getSaveFileName(self.widget, "Save StructLens project", "structlens_project.json", "StructLens project (*.json)")
        if not path:
            return
        result = self.model.analysis
        reference_source = self.reference_edit.text().strip() or None
        target_source = self.target_edit.text().strip() or None
        project = ProjectState(
            reference_source=(None if self.reference_object_name else reference_source),
            target_sources=(
                () if self.target_object_name is not None else ((target_source,) if target_source else ())
            ),
            settings=self._settings(),
            analysis_results=(result,) if result is not None else (),
            visualization_state=_state_dict(self._state_from_controls()),
            source_objects={
                key: value
                for key, value in (
                    ("reference", self.reference_object_name),
                    ("target", self.target_object_name),
                )
                if value is not None
            },
            comparison_mode=ComparisonMode(str(self.comparison_combo.currentData() or ComparisonMode.PAIRWISE.value)),
        ).with_source_hashes()
        try:
            project.save(path)
            self._set_status(f"Project saved · {Path(path).name}")
        except OSError as exc:
            self._show_error(f"Could not save project: {exc}")

    def _open_project(self) -> None:
        path, _ = self.w.QFileDialog.getOpenFileName(self.widget, "Open StructLens project", "", "StructLens project (*.json)")
        if not path:
            return
        try:
            project = ProjectState.load(path)
            self._apply_project_settings(project.settings)
            self._set_combo_value(self.comparison_combo, project.comparison_mode.value)
            self._apply_visualization_payload(project.visualization_state)
            if project.reference_source:
                self._load_source("reference", Path(project.reference_source))
            elif project.source_objects.get("reference"):
                if not self._load_named_pymol_object(
                    "reference", project.source_objects["reference"]
                ):
                    self._show_error(
                        "The saved reference PyMOL object is not available in this session."
                    )
            if project.target_sources:
                self._load_source("target", Path(project.target_sources[0]))
            elif project.source_objects.get("target"):
                if not self._load_named_pymol_object(
                    "target", project.source_objects["target"]
                ):
                    self._show_error(
                        "The saved target PyMOL object is not available in this session."
                    )
            if project.analysis_results:
                self.model = self.model.with_analysis(project.analysis_results[-1])
                self._populate_result(project.analysis_results[-1])
            self._set_status(f"Project opened · {Path(path).name}")
        except (OSError, ValueError, BundleValidationError) as exc:
            self._show_error(f"Could not open project: {exc}")

    def _apply_project_settings(self, settings: AnalysisSettings) -> None:
        self._set_combo_value(self.mode_combo, settings.alignment_mode.value)
        self.identity_spin.setValue(settings.minimum_sequence_identity)
        self.coverage_spin.setValue(settings.minimum_sequence_coverage)
        self.refined_check.setChecked(settings.refined_rmsd)
        self.cutoff_spin.setValue(settings.refinement_cutoff_angstrom)
        self.usalign_edit.setText(settings.usalign_executable or "")

    def _apply_visualization_payload(self, payload: Any) -> None:
        if not payload:
            return
        state = VisualizationState(
            highlight_filter=HighlightFilter(str(payload.get("highlight_filter", HighlightFilter.ALL.value))),
            color_mode=ColorMode(str(payload.get("color_mode", ColorMode.MUTATION_STATUS.value))),
            representation=Representation(str(payload.get("representation", Representation.STICKS.value))),
            show_labels=bool(payload.get("show_labels", False)),
            show_reference=bool(payload.get("show_reference", True)),
            show_target=bool(payload.get("show_target", True)),
            local_radius_angstrom=float(payload.get("local_radius_angstrom", 5.0)),
            preset=str(payload.get("preset", "Minimal")),
        )
        self._set_combo_value(self.preset_combo, state.preset)
        self._set_combo_value(self.filter_combo, state.highlight_filter.value)
        self._set_combo_value(self.color_combo, state.color_mode.value)
        self._set_combo_value(self.representation_combo, state.representation.value)
        self.radius_spin.setValue(state.local_radius_angstrom)
        self.labels_check.setChecked(state.show_labels)
        self.reference_check.setChecked(state.show_reference)
        self.target_check.setChecked(state.show_target)
        self.model = self.model.with_visualization(state)
        self._update_legend()

    # ---------------------------------------------------------------- exports
    def set_v03_bundle_payloads(
        self,
        *,
        msa_summary: Mapping[str, Any] | None = None,
        conservation: Mapping[str, Any] | None = None,
        interactions: Mapping[str, Any] | None = None,
        sites: Mapping[str, Any] | None = None,
        evidence: Mapping[str, Any] | None = None,
        vectors: Mapping[str, Any] | None = None,
    ) -> None:
        """Stage authoritative v0.3 payloads for the next PyMOL export.

        The GUI deliberately accepts already calculated, JSON-ready mappings
        rather than invoking a scientific service.  ``None`` means that the
        corresponding analysis is unavailable and the bundle writer will
        omit that optional entry instead of fabricating an empty result.
        Shallow copies prevent later top-level mutations by a producer from
        changing the staged export.
        """

        payloads = {
            "msa_summary": msa_summary,
            "conservation": conservation,
            "interactions": interactions,
            "sites": sites,
            "evidence": evidence,
            "vectors": vectors,
        }
        self._v03_bundle_payloads = {
            name: None if payload is None else dict(payload)
            for name, payload in payloads.items()
        }

    def _v03_bundle_kwargs(self) -> dict[str, Mapping[str, Any] | None]:
        """Return optional v0.3 payloads without calculating or normalizing them."""

        return dict(self._v03_bundle_payloads)

    def _export_xlsx(self) -> None:
        self._export("xlsx", export_analysis_xlsx, "XLSX")

    def _export_csv(self) -> None:
        self._export("csv", export_analysis_csv, "CSV")

    def _export_json(self) -> None:
        self._export("json", export_analysis_json, "JSON")

    def _export_pymol_bundle(self) -> None:
        result = self.model.analysis
        if result is None:
            self._show_error("Run a comparison before exporting a PyMOL bundle.")
            return
        if self.reference_structure is None or self.target_structure is None:
            self._show_error("Load reference and target coordinate files before exporting a PyMOL bundle.")
            return
        default_name = f"StructLens_{result.reference_id}_pairwise.structlens-pymol"
        path, _ = self.w.QFileDialog.getSaveFileName(
            self.widget,
            "Export StructLens-PyMOL bundle",
            default_name,
            "StructLens-PyMOL bundle (*.structlens-pymol)",
        )
        if not path:
            return
        try:
            write_pymol_bundle(
                path,
                reference=self.reference_structure,
                targets={result.target_id: self.target_structure},
                analysis=result,
                provenance=dict(result.provenance),
                **self._v03_bundle_kwargs(),
            )
            self._set_status(f"Validated PyMOL bundle written · {Path(path).name}")
        except (OSError, ValueError, BundleValidationError) as exc:
            self._show_error(f"Could not export PyMOL bundle: {exc}")

    def _open_in_pymol(self) -> None:
        result = self.model.analysis
        if result is None:
            self._show_error("Run a comparison before opening PyMOL.")
            return
        if self.reference_structure is None or self.target_structure is None:
            self._show_error("Load reference and target coordinate files before opening PyMOL.")
            return
        with NamedTemporaryFile(prefix="structlens_", suffix=".structlens-pymol", delete=False) as handle:
            bundle_path = Path(handle.name)
        self._temporary_paths.append(bundle_path)
        try:
            write_pymol_bundle(
                bundle_path,
                reference=self.reference_structure,
                targets={result.target_id: self.target_structure},
                analysis=result,
                provenance=dict(result.provenance),
                **self._v03_bundle_kwargs(),
            )
            launcher = PyMOLLauncher(self.pymol_edit.text().strip() or None)
            launch = launcher.launch_bundle(bundle_path)
            self._set_status(
                f"PyMOL launched with validated bundle · {bundle_path.name}"
                if launch.process_id is not None
                else f"Validated bundle prepared · {bundle_path.name}"
            )
            self._refresh_pymol_status()
        except Exception as exc:
            bundle_path.unlink(missing_ok=True)
            self._temporary_paths = [path for path in self._temporary_paths if path != bundle_path]
            self._show_error(f"Could not open in PyMOL: {exc}")

    def _export_chart_xlsx(self) -> None:
        result = self.model.analysis
        if result is None:
            self._show_error("Run a comparison before exporting chart data.")
            return
        if self.chart_combo.currentText() != "Structural deviation profile":
            self._show_error("Select Structural deviation profile before exporting chart data.")
            return
        path, _ = self.w.QFileDialog.getSaveFileName(
            self.widget, "Export chart data", "structlens_chart.xlsx", "XLSX (*.xlsx)"
        )
        if not path:
            return
        try:
            export_chart_xlsx(structural_deviation_profile(result), path)
            self._set_status(f"Chart data exported · {Path(path).name}")
        except (OSError, ValueError) as exc:
            self._show_error(f"Could not export chart data: {exc}")

    def _export_chart_image(self, suffix: str, dpi: int) -> None:
        result = self.model.analysis
        if result is None:
            self._show_error("Run a comparison before exporting a chart image.")
            return
        if self.chart_combo.currentText() != "Structural deviation profile":
            self._show_error("Select Structural deviation profile before exporting a chart image.")
            return
        path, _ = self.w.QFileDialog.getSaveFileName(
            self.widget,
            f"Export {suffix.upper()} chart ({dpi} dpi)",
            f"structlens_chart.{suffix}",
            f"{suffix.upper()} (*.{suffix})",
        )
        if not path:
            return
        try:
            export_chart_image(structural_deviation_profile(result), path, dpi=dpi)
            self._set_status(f"Chart image exported · {Path(path).name} · {dpi} dpi")
        except (OSError, RuntimeError, ValueError) as exc:
            self._show_error(f"Could not export chart image: {exc}")

    def _export(self, suffix: str, exporter: Any, label: str) -> None:
        result = self.model.analysis
        if result is None:
            self._show_error("Run a comparison before exporting results.")
            return
        path, _ = self.w.QFileDialog.getSaveFileName(self.widget, f"Export {label}", f"structlens_result.{suffix}", f"{label} (*.{suffix})")
        if not path:
            return
        try:
            exporter(result, path)
            self._set_status(f"{label} export written · {Path(path).name}")
        except (OSError, ValueError) as exc:
            self._show_error(f"Could not export {label}: {exc}")

    # --------------------------------------------------------------- feedback
    def _set_status(self, message: str) -> None:
        self.header_status.setText("Busy" if self.model.busy else "Ready")
        self.footer_status.setText(message)

    def _show_error(self, message: str) -> None:
        self.model = self.model.with_error(message)
        self._set_busy(False)
        self._set_status(message)
        self.header_status.setText("Needs attention")
        if hasattr(self.w.QMessageBox, "warning"):
            self.w.QMessageBox.warning(self.widget, "StructLens", message)

    def _update_mode_help(self, label: str) -> None:
        key = label.split(" ·", 1)[0]
        self.mode_help.setText(WORKFLOW_HELP.get(key, WORKFLOW_HELP["Auto"]))
        self.manual_group.setVisible(key == "Manual")

    def _update_comparison_help(self, label: str) -> None:
        key = label.split(" ·", 1)[0]
        self.comparison_help.setText(
            {
                "Pairwise": "One reference and one target. Use it for detailed residue inspection and focused WT/mutant analysis.",
            }.get(key, "Choose a comparison topology to see its outputs and computational implications.")
        )

    def _update_chart_explanation(self, label: str) -> None:
        explanations = {
            "Structural deviation profile": "Reference position on X; selectable Cα displacement, backbone RMSD, side-chain RMSD, or local RMSD on Y (Å).",
            "Mutation / conservation matrix": "Rows are structures and columns are reference-aligned positions; cell text preserves mutation identity.",
            "Pairwise similarity heatmap": "Mirrors one stored value per pair for sequence identity, TM-score, RMSD, or key-site RMSD.",
            "Sequence–structure relationship": "Sequence identity (%) versus TM-score by default; points link back to a selected pair.",
            "Structural conservation profile": "Cα positional variability (Å) with a separate position-coverage track.",
            "Key-residue comparison": "Explicit key reference residues compared across targets with units and missing mappings visible.",
            "MSA conservation profile": "Alignment-column conservation with gap and ambiguous-residue tracks; gaps are not amino acids.",
            "Sequence logo": "Letter heights are pᵢ × information content, with canonical amino acids only.",
            "Interaction difference matrix": "Reference-normalized conserved, gained, lost, and target-only-unmapped interactions.",
            "Site comparison": "Global-frame and site-fitted RMSD remain separate; envelope volume is not a cavity volume.",
            "Distance-difference heatmap": "Target internal distance minus reference internal distance (Å), with missing pairs masked.",
            "Evidence Card figure": "Publication-ready descriptive evidence; no impact or functional score is inferred.",
        }
        self.chart_explanation.setText(explanations.get(label, "Charts consume the authoritative analysis state."))
        self._update_chart_export_state(label)

    def _define_site_from_controls(self) -> None:
        positions = tuple(item.strip() for item in self.site_residues_edit.text().split(",") if item.strip())
        mode = str(self.site_mode_combo.currentData())
        if mode == "key_residues" and not positions:
            self._show_error("Key-residue sites need at least one reference position.")
            return
        self.site_status_label.setText(
            f"Site definition staged · mode={mode} · {len(positions)} reference position(s) · radius {self.site_radius_spin.value():.1f} Å. "
            "Run the site service to calculate coverage, RMSDs, SASA, and interaction fingerprints."
        )

    def _update_chart_export_state(self, label: str) -> None:
        """Keep export controls aligned with the chart profile they produce."""

        supported = label == "Structural deviation profile"
        for button in self.chart_export_buttons:
            button.setEnabled(supported)
            if supported:
                button.setToolTip("Export the selected structural deviation profile.")
            else:
                button.setToolTip(
                    "This chart profile is available as structured data in the application API; "
                    "GUI image/XLSX export is currently implemented for Structural deviation profile."
                )

    def _set_combo_value(self, combo: Any, value: str) -> None:
        for index in range(combo.count()):
            if str(combo.itemData(index)) == value or combo.itemText(index) == value:
                combo.setCurrentIndex(index)
                return

    def _refresh_pymol_status(self) -> None:
        try:
            executable = PyMOLLauncher(self.pymol_edit.text().strip() or None).locate()
        except Exception:
            self.pymol_status.setText(
                "Ready to export · PyMOL executable not configured or not found. Export remains available even without PyMOL."
            )
            return
        self.pymol_status.setText(
            f"Ready to open in PyMOL · executable resolved to {executable.name}. Scientific calculations stay in StructLens."
        )

    def close(self) -> None:
        if self._future is not None:
            self._cancel_event.set()
            self._future.cancel()
            self._future = None
        if self._poll_timer is not None:
            self._poll_timer.stop()
            self._poll_timer = None
        self._pymol.reset()
        for path in self._temporary_paths:
            path.unlink(missing_ok=True)


def _stylesheet() -> str:
    return """
    QWidget#structlensPanel { background: #111827; color: #dbe7f3; }
    QFrame#header { background: #0c1421; border-bottom: 1px solid #26364b; }
    QFrame#sidebar { background: #0d1726; border-right: 1px solid #26364b; }
    QFrame#footer { background: #0c1421; border-top: 1px solid #26364b; }
    QLabel#eyebrow, QLabel#sectionKicker, QLabel#fieldLabel { color: #7f9ab8; font-size: 10px; font-weight: 700; }
    QLabel#eyebrow, QLabel#sectionKicker, QLabel#fieldLabel { letter-spacing: 1px; }
    QLabel#windowTitle { color: #f5f8fc; font-size: 20px; font-weight: 700; }
    QLabel#pageTitle { color: #f5f8fc; font-size: 24px; font-weight: 700; }
    QLabel#pagePurpose { color: #9eb0c5; font-size: 13px; }
    QLabel#statusPill { background: #19304a; color: #9fc7ff; border: 1px solid #2a5a8a; border-radius: 12px; padding: 5px 11px; font-weight: 700; }
    QLabel#footerStatus { color: #b8c9dc; }
    QLabel#footerMeta, QLabel#fieldMeta, QLabel#sidebarNote { color: #71869e; }
    QLabel#sidebarNote { font-size: 11px; line-height: 1.3; }
    QLabel#helpText, QLabel#inlineNote { color: #a8bbd0; line-height: 1.35; }
    QLabel#legend { background: #172438; color: #c3d6ea; border: 1px solid #2e435d; padding: 11px; }
    QLabel#resultDecision { background: #14253a; color: #d9e8f8; border: 1px solid #2a537d; padding: 13px; }
    QLabel#metricValue { color: #f5f8fc; font-size: 17px; font-weight: 700; }
    QListWidget#workflowNav { background: transparent; border: none; outline: none; color: #aabbd0; }
    QListWidget#workflowNav::item { padding: 11px 12px; border-radius: 6px; }
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
    """


def _label(widgets: Any, text: str, object_name: str) -> Any:
    label = widgets.QLabel(text)
    label.setObjectName(object_name)
    return label


def _button(widgets: Any, text: str, object_name: str) -> Any:
    button = widgets.QPushButton(text)
    button.setObjectName(object_name)
    return button


def _fraction_spin(spin: Any, value: float) -> None:
    spin.setRange(0.0, 1.0)
    spin.setSingleStep(0.05)
    spin.setDecimals(2)
    spin.setValue(value)
    spin.setSuffix(" fraction")


def _configure_table(table: Any, widgets: Any) -> None:
    table.setAlternatingRowColors(True)
    table.setSortingEnabled(False)
    view = widgets.QAbstractItemView
    table.setSelectionBehavior(getattr(view, "SelectionBehavior", view).SelectRows)
    table.setSelectionMode(getattr(view, "SelectionMode", view).SingleSelection)
    table.setEditTriggers(getattr(view, "EditTrigger", view).NoEditTriggers)
    table.verticalHeader().setVisible(False)
    table.horizontalHeader().setStretchLastSection(True)
    table.setMinimumHeight(220)


def _page_subtitle(section: str, *, standalone: bool = False) -> str:
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


def _structure_meta(structure: ProteinStructure) -> str:
    count = sum(len(chain.residue_records or chain.residues) for chain in structure.chains)
    return f"{len(structure.chains)} chain(s) · {count} residues · {structure.structure_id}"


def _find_residue(chain: ProteinChain, token: str) -> ResidueId | None:
    match = _RESIDUE_TOKEN.match(token)
    if match is None:
        return None
    chain_id = match.group("chain")
    auth_seq_id = match.group("number")
    insertion_code = match.group("insertion") or None
    for residue in chain.residue_records:
        residue_id = residue.residue_id
        if (
            residue_id.chain_id == chain_id
            and residue_id.auth_seq_id == auth_seq_id
            and residue_id.insertion_code == insertion_code
        ):
            return residue_id
    for residue_id in chain.residues:
        if (
            residue_id.chain_id == chain_id
            and residue_id.auth_seq_id == auth_seq_id
            and residue_id.insertion_code == insertion_code
        ):
            return residue_id
    return None


def _residue_label(residue: ResidueId | None) -> str:
    if residue is None:
        return "—"
    insertion = residue.insertion_code or ""
    return f"{residue.chain_id}:{residue.auth_seq_id}{insertion} {residue.residue_name}"


def _number(value: float | int | None) -> str:
    return "—" if value is None else f"{value:.3f}" if isinstance(value, float) else str(value)


def _human(value: str) -> str:
    return value.replace("_", " ").replace("ca", "Cα").title()


def _state_dict(state: VisualizationState) -> dict[str, Any]:
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


__all__ = ["PanelController", "build_panel"]
