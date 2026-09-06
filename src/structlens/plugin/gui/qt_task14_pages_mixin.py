"""Qt construction for the report-backed Task 14 evidence panels."""

# ruff: noqa: F403,F405

from __future__ import annotations

from .qt_context import *  # noqa: F401,F403


class Task14PageMixin(QtMixinContext):
    """Build only the widgets that project canonical quality/pocket evidence."""

    def _build_pocket_panel(self, content: Any) -> None:
        """Build the report-backed Sites & Pockets evidence panel.

        Detection and volume measurement are application services owned by the
        canonical Compare report.  This panel renders that report and never
        starts a second GUI-local calculation.
        """

        group = self.w.QGroupBox("Sites & Pockets", self.widget)
        layout = self.w.QVBoxLayout(group)
        layout.setContentsMargins(18, 16, 18, 16)
        layout.setSpacing(10)
        layout.addWidget(_label(self.w, "Pocket evidence", "sectionKicker"))
        layout.addWidget(
            _label(
                self.w,
                "Pocket candidates, matching states, and coarse/fine free-volume estimates are read from the canonical Compare report.",
                "helpText",
            )
        )
        self.pocket_status_label = _label(self.w, "Pocket evidence appears after Compare.", "inlineNote")
        self.pocket_status_label.setWordWrap(True)
        layout.addWidget(self.pocket_status_label)
        self.pocket_capability_label = _label(self.w, "No canonical pocket report loaded.", "inlineNote")
        self.pocket_capability_label.setWordWrap(True)
        layout.addWidget(self.pocket_capability_label)
        self.pocket_table = self.w.QTableWidget(0, 7, group)
        self.pocket_table.setHorizontalHeaderLabels(
            [
                "Candidate",
                "Structure",
                "Match state",
                "α spheres",
                "Lining residues",
                "Fine volume (Å³)",
                "Sensitivity (Å³)",
            ]
        )
        _configure_table(self.pocket_table, self.w)
        self.pocket_table.cellClicked.connect(self._select_pocket_candidate)
        layout.addWidget(self.pocket_table)
        self.pocket_detail_label = _label(self.w, "Select a pocket candidate to inspect volume evidence.", "inlineNote")
        self.pocket_detail_label.setWordWrap(True)
        layout.addWidget(self.pocket_detail_label)
        actions = self.w.QHBoxLayout()
        self.pocket_detect_button = _button(self.w, "Detect pockets (rerun Compare)", "secondaryButton")
        self.pocket_detect_button.setToolTip(
            "Pocket detection is computed by Compare and is read-only in this report view."
        )
        self.pocket_detect_button.setEnabled(False)
        self.pocket_detect_button.clicked.connect(self._request_pocket_detection)
        actions.addWidget(self.pocket_detect_button)
        self.pocket_measure_button = _button(self.w, "Measure volume (rerun Compare)", "secondaryButton")
        self.pocket_measure_button.setToolTip(
            "Volume measurement is computed by Compare and is read-only in this report view."
        )
        self.pocket_measure_button.setEnabled(False)
        self.pocket_measure_button.clicked.connect(self._request_pocket_measurement)
        actions.addWidget(self.pocket_measure_button)
        actions.addStretch(1)
        layout.addLayout(actions)
        content.addWidget(group)

    def _build_quality_panel(self, content: Any) -> None:
        """Build the report-backed coordinate quality evidence panel."""

        group = self.w.QGroupBox("Quality", self.widget)
        layout = self.w.QVBoxLayout(group)
        layout.setContentsMargins(18, 16, 18, 16)
        layout.setSpacing(10)
        layout.addWidget(_label(self.w, "Quality", "sectionKicker"))
        controls = self.w.QHBoxLayout()
        controls.addWidget(_label(self.w, "FILTER", "fieldLabel"))
        self.quality_filter = self.w.QComboBox(group)
        for label, value in (
            ("All diagnostics", QualityFilter.ALL.value),
            ("Errors", QualityFilter.ERRORS.value),
            ("Warnings", QualityFilter.WARNINGS.value),
            ("Info", QualityFilter.INFO.value),
        ):
            self.quality_filter.addItem(label, value)
        self.quality_filter.currentIndexChanged.connect(lambda _index: self._render_quality())
        controls.addWidget(self.quality_filter)
        controls.addStretch(1)
        layout.addLayout(controls)
        self.quality_summary_label = _label(self.w, "No QC report loaded.", "inlineNote")
        self.quality_summary_label.setWordWrap(True)
        layout.addWidget(self.quality_summary_label)
        self.quality_table = self.w.QTableWidget(0, 6, group)
        self.quality_table.setHorizontalHeaderLabels(
            ["Structure", "Severity", "Code", "Location", "Message", "Remediation"]
        )
        _configure_table(self.quality_table, self.w)
        layout.addWidget(self.quality_table)
        layout.addWidget(
            _label(
                self.w,
                "Diagnostics are descriptive evidence from the quality service; no score or causal interpretation is inferred here.",
                "helpText",
            )
        )
        content.addWidget(group)


__all__ = ["Task14PageMixin"]
