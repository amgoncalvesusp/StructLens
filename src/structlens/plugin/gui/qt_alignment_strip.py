"""Windowed, report-indexed sequence alignment for linked residue inspection."""

from __future__ import annotations

from html import escape
from typing import Any

from structlens.core.models import ResidueCorrespondence
from structlens.core.reports import CorrespondenceSnapshot

from .qt_context import QtMixinContext, _label, _residue_label

_COLUMNS = 24


class AlignmentStripMixin(QtMixinContext):
    def _add_alignment_strip(self, content: Any) -> None:
        group = self.w.QGroupBox("Aligned positions")
        layout = self.w.QVBoxLayout(group)
        controls = self.w.QHBoxLayout()
        controls.addWidget(_label(self.w, "Window", "fieldLabel"))
        self.alignment_window = self.w.QSpinBox()
        self.alignment_window.setAccessibleName("Alignment window")
        self.alignment_window.setRange(1, 1)
        controls.addWidget(self.alignment_window)
        self.alignment_range = _label(self.w, "", "inlineNote")
        controls.addWidget(self.alignment_range)
        controls.addStretch()
        layout.addLayout(controls)
        self.alignment_strip = self.w.QTableWidget(2, 0)
        self.alignment_strip.setAccessibleName("Pairwise sequence alignment")
        self.alignment_strip.setVerticalHeaderLabels(["Reference", "Target"])
        view = self.w.QAbstractItemView
        self.alignment_strip.setEditTriggers(getattr(view, "EditTrigger", view).NoEditTriggers)
        self.alignment_strip.setSelectionMode(getattr(view, "SelectionMode", view).SingleSelection)
        self.alignment_strip.setSelectionBehavior(getattr(view, "SelectionBehavior", view).SelectColumns)
        self.alignment_strip.setFixedHeight(110)
        self.alignment_strip.cellClicked.connect(self._inspect_alignment_column)
        self.alignment_strip.cellActivated.connect(self._inspect_alignment_column)
        self.alignment_window.valueChanged.connect(self._draw_alignment_window)
        layout.addWidget(self.alignment_strip)
        layout.addWidget(_label(self.w, "Numbers are alignment columns. Select a column to inspect both residue identifiers. A dash marks a gap.", "helpText"))
        content.addWidget(group)
        self._set_alignment_rows(())

    def _set_alignment_rows(self, rows: tuple[CorrespondenceSnapshot | ResidueCorrespondence, ...]) -> None:
        self._alignment_rows = rows
        self._alignment_focus_index: int | None = None
        pages = max(1, (len(rows) + _COLUMNS - 1) // _COLUMNS)
        self.alignment_window.setMaximum(pages)
        self.alignment_window.setEnabled(pages > 1)
        self.alignment_window.setValue(1)
        self._draw_alignment_window()

    def _draw_alignment_window(self, *_: object) -> None:
        start = (self.alignment_window.value() - 1) * _COLUMNS
        rows = self._alignment_rows[start:start + _COLUMNS]
        self.alignment_strip.setCurrentCell(-1, -1)
        self.alignment_strip.setColumnCount(len(rows))
        self.alignment_strip.setHorizontalHeaderLabels([str(start + column + 1) for column in range(len(rows))])
        self.alignment_range.setText(f"{start + 1}–{start + len(rows)} of {len(self._alignment_rows)}" if rows else "No alignment yet")
        for column, row in enumerate(rows):
            tooltip = f"{_residue_label(row.reference)} → {_residue_label(row.target)} · {row.status.value}"
            for side, (residue, letter) in enumerate(((row.reference, row.reference_one_letter), (row.target, row.target_one_letter))):
                item = self.w.QTableWidgetItem((letter or "X") if residue is not None else "–")
                item.setToolTip(f"<span>{escape(tooltip)}</span>")
                if row.status.value in {"substitution", "insertion", "deletion"}:
                    font = item.font()
                    font.setBold(True)
                    item.setFont(font)
                self.alignment_strip.setItem(side, column, item)
            self.alignment_strip.setColumnWidth(column, 36)
            if row.alignment_index == self._alignment_focus_index:
                self.alignment_strip.setCurrentCell(0, column)

    def _inspect_alignment_column(self, _side: int, column: int) -> None:
        index = (self.alignment_window.value() - 1) * _COLUMNS + column
        if 0 <= index < len(self._alignment_rows):
            self._open_residue_filter("all", self._alignment_rows[index].alignment_index)

    def _sync_alignment_selection(self, alignment_index: int | None) -> None:
        self._alignment_focus_index = alignment_index
        if alignment_index is None:
            self.alignment_strip.setCurrentCell(-1, -1)
            return
        for position, row in enumerate(self._alignment_rows):
            if row.alignment_index == alignment_index:
                page, column = divmod(position, _COLUMNS)
                self.alignment_window.setValue(page + 1)
                self.alignment_strip.setCurrentCell(0, column)
                self.alignment_strip.scrollToItem(self.alignment_strip.item(0, column))
                return
