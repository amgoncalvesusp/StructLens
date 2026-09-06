"""Focused pages presentation/routing mixin for the Qt panel."""

# ruff: noqa: F403,F405

from __future__ import annotations

from .qt_context import *  # noqa: F401,F403


class PageMixin(QtMixinContext):
    """Cohesive GUI-only pages behavior composed into PanelController."""

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
        chain_layout.addWidget(_label(self.w, "REFERENCE MODEL", "fieldLabel"), 0, 0)
        self.reference_model_combo = self.w.QComboBox(chain_group)
        chain_layout.addWidget(self.reference_model_combo, 1, 0)
        chain_layout.addWidget(_label(self.w, "TARGET MODEL", "fieldLabel"), 0, 1)
        self.target_model_combo = self.w.QComboBox(chain_group)
        chain_layout.addWidget(self.target_model_combo, 1, 1)
        chain_layout.addWidget(_label(self.w, "REFERENCE CHAIN", "fieldLabel"), 2, 0)
        self.reference_chain_combo = self.w.QComboBox(chain_group)
        chain_layout.addWidget(self.reference_chain_combo, 3, 0)
        chain_layout.addWidget(_label(self.w, "TARGET CHAIN", "fieldLabel"), 2, 1)
        self.target_chain_combo = self.w.QComboBox(chain_group)
        chain_layout.addWidget(self.target_chain_combo, 3, 1)
        self.reference_model_combo.currentIndexChanged.connect(lambda _: self._model_changed("reference"))
        self.target_model_combo.currentIndexChanged.connect(lambda _: self._model_changed("target"))
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
        backend_layout.addWidget(
            _label(self.w, " · ".join(f"{key}: {value}" for key, value in versions.items()), "fieldMeta")
        )
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
        # The header Compare action is the single primary submit control.
        # Keep the old attribute as an adapter for integrations that inspect
        # busy state, without creating a second in-page action.
        self.run_button = self.compare_button
        content.addLayout(run_layout)

        structure_results = self.w.QGroupBox("Latest structure comparison", self.widget)
        structure_results_layout = self.w.QVBoxLayout(structure_results)
        structure_results_layout.setContentsMargins(18, 16, 18, 16)
        self.structure_result_summary = _label(
            self.w,
            "No structure comparison result yet. Run comparison to populate this tab.",
            "inlineNote",
        )
        self.structure_result_summary.setWordWrap(True)
        structure_results_layout.addWidget(self.structure_result_summary)
        self.structure_result_table = self.w.QTableWidget(0, 9, structure_results)
        self.structure_result_table.setHorizontalHeaderLabels(
            [
                "Reference",
                "Target",
                "Decision",
                "Strict RMSD (Å)",
                "Refined RMSD (Å)",
                "TM-score",
                "Mapped",
                "Excluded",
                "Backend",
            ]
        )
        _configure_table(self.structure_result_table, self.w)
        self.structure_result_table.setMinimumHeight(92)
        structure_results_layout.addWidget(self.structure_result_table)
        content.addWidget(structure_results)

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
        self.sequence_chart_status = _label(
            self.w,
            "Sequence conservation chart unavailable until an authoritative MSA or comparison result is present.",
            "inlineNote",
        )
        self.sequence_chart_status.setWordWrap(True)
        msa_layout.addWidget(self.sequence_chart_status)
        self.sequence_chart_layout = self.w.QVBoxLayout()
        self.sequence_chart_layout.setContentsMargins(0, 0, 0, 0)
        msa_layout.addLayout(self.sequence_chart_layout, 1)
        msa_layout.addWidget(
            _label(
                self.w,
                "Alignment columns retain reference-relative insertion labels; conservation excludes gaps and ambiguous residues from entropy.",
                "helpText",
            )
        )
        content.addWidget(msa_group, 1)

    def _render_msa_result(self, alignment: MultipleSequenceAlignment | None) -> None:
        """Render an authoritative MSA without recalculating its values."""

        self.msa_table.setRowCount(0)
        if alignment is None:
            self.msa_summary_label.setText("MSA unavailable.")
            self._msa_chart_dataset = None
            self._clear_chart_layout(self.sequence_chart_layout)
            self.sequence_chart_status.setText(
                "Sequence conservation chart unavailable until an authoritative MSA or comparison result is present."
            )
            self._update_chart_export_state(self.chart_combo.currentText())
            return
        self.msa_table.setRowCount(len(alignment.aligned_rows))
        for row_index, (structure_id, row) in enumerate(alignment.aligned_rows):
            source = next(
                (sequence.source for sequence in alignment.sequences if sequence.structure_id == structure_id),
                "unknown",
            )
            self.msa_table.setItem(row_index, 0, self.w.QTableWidgetItem(structure_id))
            self.msa_table.setItem(row_index, 1, self.w.QTableWidgetItem(row))
            self.msa_table.setItem(row_index, 2, self.w.QTableWidgetItem(source))
        self.msa_summary_label.setText(
            f"{len(alignment.aligned_rows)} sequences · {len(alignment.columns)} alignment columns · "
            f"{sum(column.reference_residue is None for column in alignment.columns)} reference-relative insertion columns."
        )
        self._render_msa_chart(alignment)
        self._update_chart_export_state(self.chart_combo.currentText())
        self._render_selected_chart()
        self.nav.setCurrentRow(_SEQUENCES_PAGE_INDEX)

    def set_msa_result(self, alignment: MultipleSequenceAlignment | None) -> None:
        """Compatibility adapter for callers that provide an MSA directly."""

        self._render_msa_result(alignment)

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
        self.site_ligand_edit = self.w.QLineEdit(site_group)
        self.site_ligand_edit.setPlaceholderText("Ligand ID (for ligand radius)")
        site_layout.addWidget(_label(self.w, "LIGAND ID", "fieldLabel"), 0, 3)
        site_layout.addWidget(self.site_ligand_edit, 1, 3)
        self.site_radius_spin = self.w.QDoubleSpinBox(site_group)
        self.site_radius_spin.setRange(0.1, 20.0)
        self.site_radius_spin.setValue(5.0)
        self.site_radius_spin.setSuffix(" Å")
        site_layout.addWidget(_label(self.w, "RADIUS", "fieldLabel"), 0, 2)
        site_layout.addWidget(self.site_radius_spin, 1, 2)
        self.site_define_button = _button(self.w, "Define site", "secondaryButton")
        self.site_define_button.clicked.connect(self._define_site_from_controls)
        site_layout.addWidget(self.site_define_button, 1, 4)
        self.site_status_label = _label(
            self.w,
            "Site metrics appear after an analysis service run; this panel never estimates them locally.",
            "inlineNote",
        )
        self.site_status_label.setWordWrap(True)
        site_layout.addWidget(self.site_status_label, 2, 0, 1, 5)
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
        self.site_metrics_table = self.w.QTableWidget(0, 12, metrics)
        self.site_metrics_table.setHorizontalHeaderLabels(
            [
                "Site",
                "Structure",
                "Mapped",
                "Coverage",
                "Global RMSD (Å)",
                "Site-fitted RMSD (Å)",
                "Centroid Δ (Å)",
                "R gyration (Å)",
                "Envelope (Å³)",
                "SASA (Å²)",
                "Polar",
                "Charged",
            ]
        )
        _configure_table(self.site_metrics_table, self.w)
        self.site_metrics_table.setMinimumHeight(108)
        metrics_layout.addWidget(self.site_metrics_table)
        content.addWidget(metrics)
        self._build_pocket_panel(content)

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
        chart_preview = self.w.QGroupBox("Authoritative chart preview", self.widget)
        chart_preview_layout = self.w.QVBoxLayout(chart_preview)
        chart_preview_layout.setContentsMargins(18, 16, 18, 16)
        self.chart_preview_status = _label(
            self.w,
            "Chart unavailable. Run the corresponding scientific service first.",
            "inlineNote",
        )
        self.chart_preview_status.setWordWrap(True)
        chart_preview_layout.addWidget(self.chart_preview_status)
        self.chart_preview_layout = self.w.QVBoxLayout()
        self.chart_preview_layout.setContentsMargins(0, 0, 0, 0)
        chart_preview_layout.addLayout(self.chart_preview_layout, 1)
        content.addWidget(chart_preview, 1)
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
            signal = (
                getattr(control, "currentIndexChanged", None)
                or getattr(control, "valueChanged", None)
                or getattr(control, "toggled", None)
            )
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
            ("TSV", self._export_tsv),
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
        self._build_presentation_sections(content)
        self._build_quality_panel(content)
        history = self.w.QGroupBox("Compiled analysis history", self.widget)
        history_layout = self.w.QVBoxLayout(history)
        history_layout.setContentsMargins(18, 16, 18, 16)
        self.results_history_status = _label(
            self.w,
            "No completed analyses yet. Results from each comparison will be retained here.",
            "inlineNote",
        )
        self.results_history_status.setWordWrap(True)
        history_layout.addWidget(self.results_history_status)
        self.results_table = self.w.QTableWidget(0, 10, history)
        self.results_table.setHorizontalHeaderLabels(
            [
                "Reference",
                "Target",
                "Decision",
                "Identity",
                "Coverage",
                "Strict RMSD (Å)",
                "Refined RMSD (Å)",
                "TM-score",
                "Mapped",
                "Mutations",
            ]
        )
        _configure_table(self.results_table, self.w)
        self.results_table.setMinimumHeight(130)
        history_layout.addWidget(self.results_table, 1)
        content.addWidget(history, 1)

    # --------------------------------------------------------------- bindings


__all__ = ["PageMixin"]
