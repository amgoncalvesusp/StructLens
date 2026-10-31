"""Residue table, evidence inspector and navigation from mutation counts."""

from __future__ import annotations

from collections import Counter
from dataclasses import replace
from typing import Any

from structlens.core.models.mutation import SEQUENCE_CHANGE_KINDS, MutationEvent, MutationKind

from .main_panel import SCIENTIFIC_SECTIONS
from .qt_context import QtMixinContext, _button, _configure_table, _label


class ResiduePageMixin(QtMixinContext):
    def _build_residues_page(self) -> None:
        _, content = self._add_page(
            "Residues", "Select a position to inspect the evidence for its reference and target residues."
        )
        self.residue_summary = _label(self.w, "No comparison yet.", "inlineNote")
        content.addWidget(self.residue_summary)
        controls = self.w.QHBoxLayout()
        filter_label = _label(self.w, "Show positions", "fieldLabel")
        self.residue_filter_combo = self.w.QComboBox()
        filter_label.setBuddy(self.residue_filter_combo)
        self.residue_filter_combo.setAccessibleName("Residue filter")
        self.residue_filter_combo.addItem("All positions", "all")
        self.residue_filter_combo.addItem("Only changes", "changes")
        for kind in SEQUENCE_CHANGE_KINDS:
            self.residue_filter_combo.addItem(kind.value.capitalize() + "s", kind.value)
        controls.addWidget(filter_label)
        controls.addWidget(self.residue_filter_combo)
        controls.addStretch()
        content.addLayout(controls)

        details = self.w.QHBoxLayout()
        self.residue_table = self.w.QTableWidget(0, 10)
        self.residue_table.setAccessibleName("Residue correspondences")
        self.residue_table.setHorizontalHeaderLabels(
            ["Index", "Reference", "Target", "Status", "Cα Δ (Å)", "Backbone (Å)", "Side-chain (Å)",
             "Heavy (Å)", "Outlier", "Key"]
        )
        _configure_table(self.residue_table, self.w)
        self.residue_table.currentCellChanged.connect(
            lambda row, column, *_: self._focus_residue(row, column, focus_pymol=False)
        )
        self.residue_table.cellDoubleClicked.connect(self._focus_residue)
        details.addWidget(self.residue_table, 2)

        card = self.w.QGroupBox("Residue Evidence Card")
        card.setObjectName("residueEvidence")
        card.setMinimumWidth(180)
        card.setMinimumHeight(300)
        card_layout = self.w.QVBoxLayout(card)
        self.evidence_card_label = _label(self.w, "Select a residue to inspect its evidence.", "helpText")
        policy = self.evidence_card_label.sizePolicy()
        policy.setHorizontalPolicy(getattr(self.w.QSizePolicy, "Policy", self.w.QSizePolicy).Ignored)
        self.evidence_card_label.setSizePolicy(policy)
        self.evidence_card_label.setTextFormat(getattr(self.c.Qt, "TextFormat", self.c.Qt).PlainText)
        alignment = getattr(self.c.Qt, "AlignmentFlag", self.c.Qt)
        self.evidence_card_label.setAlignment(alignment.AlignLeft | alignment.AlignTop)
        card_layout.addWidget(self.evidence_card_label)
        details.addWidget(card, 1)
        content.addLayout(details, 1)
        self.residue_filter_combo.currentIndexChanged.connect(self._apply_residue_filter)

    def _add_mutation_breakdown(self, layout: Any) -> None:
        actions = self.w.QHBoxLayout()
        self.mutation_kind_buttons: dict[MutationKind, Any] = {}
        for kind in SEQUENCE_CHANGE_KINDS:
            button = _button(self.w, "", "secondaryButton")
            button.clicked.connect(lambda _checked=False, key=kind.value: self._open_residue_filter(key))
            actions.addWidget(button)
            self.mutation_kind_buttons[kind] = button
        layout.addLayout(actions, layout.rowCount(), 0, 1, 3)
        self._set_mutation_counts(None)

    def _set_mutation_counts(self, events: tuple[MutationEvent, ...] | None) -> None:
        counts = Counter(event.kind for event in events or ())
        total = sum(counts[kind] for kind in SEQUENCE_CHANGE_KINDS)
        total_button = self.result_labels["mutation_count"]
        total_button.setText(str(total) if events is not None else "—")
        total_button.setEnabled(bool(total))
        for kind, button in self.mutation_kind_buttons.items():
            count = counts[kind]
            button.setText(f"{kind.value.capitalize()}s ({count if events is not None else '—'})")
            button.setEnabled(bool(count))

    def _open_sequence_change(self, row: int, _column: int) -> None:
        analysis = self.model.report.analysis if self.model.report is not None else self.model.analysis
        if analysis is not None and 0 <= row < len(analysis.mutations):
            event = analysis.mutations[row]
            if event.kind in SEQUENCE_CHANGE_KINDS:
                self._open_residue_filter(event.kind.value, event.alignment_index)

    def _open_residue_filter(self, key: str, alignment_index: int | None = None) -> None:
        self.nav.setCurrentRow(SCIENTIFIC_SECTIONS.index("Residues"))
        self._set_residue_filter(key, alignment_index)
        self.residue_table.setFocus()

    def _set_residue_filter(self, key: str, alignment_index: int | None) -> None:
        combo = self.residue_filter_combo
        blocked = combo.blockSignals(True)
        try:
            combo.setCurrentIndex(combo.findData(key))
        finally:
            combo.blockSignals(blocked)
        self._apply_residue_filter(alignment_index=alignment_index)

    def _apply_residue_filter(self, _index: int = -1, *, alignment_index: int | None = None) -> None:
        table = self.residue_table
        mode = self.residue_filter_combo.currentData()
        report = self.model.report
        analysis = report.analysis if report is not None else self.model.analysis
        kinds = {e.alignment_index: e.kind for e in analysis.mutations} if analysis is not None else {}
        first = selected = None
        visible_count = 0
        for row in range(table.rowCount()):
            index = int(table.item(row, 0).text())
            kind = kinds.get(index)
            visible = mode == "all" or (kind in SEQUENCE_CHANGE_KINDS if mode == "changes" else kind == mode)
            table.setRowHidden(row, not visible)
            visible_count += int(visible)
            if visible and first is None:
                first = row
            if visible and index == alignment_index:
                selected = row
        self.residue_summary.setText(f"{visible_count} of {table.rowCount()} aligned positions shown")
        if first is None:
            table.setCurrentCell(-1, -1)
            self._selection_controller = self._selection_controller.clear_selection()
            self.model = replace(self.model, selection=None)
            self.evidence_card_label.setText("No positions match this filter.")
            return
        row = selected if selected is not None else first
        same_cell = table.currentRow() == row and table.currentColumn() == 0
        table.setCurrentCell(row, 0)
        if same_cell:
            self._focus_residue(row, 0, focus_pymol=False)
        table.scrollToItem(table.item(row, 0))
