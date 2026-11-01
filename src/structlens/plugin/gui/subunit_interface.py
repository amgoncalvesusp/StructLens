"""Reference-only subunit and peptide interface analysis in the desktop GUI."""

from __future__ import annotations

import csv
import io
from concurrent.futures import Future
from pathlib import Path
from threading import Event
from typing import Any

from structlens.application.report_exports import _safe
from structlens.application.subunit_interface import (
    InterfaceContact,
    SubunitInterface,
    analyze_subunit_interface,
    eligible_interface_chains,
)
from structlens.core.models import ResidueId
from structlens.core.reports.safe_io import atomic_write_bytes

from .loaded_source import LoadedSource
from .qt_compat import user_role
from .qt_context import _ANALYSIS_EXECUTOR, _configure_table, _residue_label

_PAGE_SIZE = 100
_MAX_CSV_BYTES = 100 * 1024 * 1024


def _residue_fields(residue: ResidueId, side: str) -> dict[str, str | None]:
    return {
        f"structure_{side}": residue.structure_id,
        f"model_{side}": residue.model_id,
        f"chain_{side}": residue.chain_id,
        f"number_{side}": residue.auth_seq_id,
        f"insertion_{side}": residue.insertion_code,
        f"name_{side}": residue.residue_name,
    }


def export_interface_csv(path: Path, result: SubunitInterface) -> None:
    fields = ["record", "source_sha256", "model", "cutoff_angstrom", "observed_a", "observed_b",
              "unavailable_a", "unavailable_b", "excluded_a", "excluded_b", "structure_a", "model_a", "chain_a", "number_a", "insertion_a",
              "name_a", "structure_b", "model_b", "chain_b", "number_b", "insertion_b", "name_b", "atom_a",
              "atom_b", "distance_angstrom", "atom_pair_count", "partner_count", "interaction_type"]
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fields)
    writer.writeheader()

    def write(**row: object) -> None:
        if output.tell() > _MAX_CSV_BYTES:
            raise ValueError("Interface CSV exceeds the supported export size.")
        writer.writerow({key: _safe(value) for key, value in row.items()})

    write(record="analysis", source_sha256=result.source_sha256, model=result.model_id,
          chain_a=result.chain_a, chain_b=result.chain_b, cutoff_angstrom=result.cutoff_angstrom,
          observed_a=result.observed_a, observed_b=result.observed_b,
          unavailable_a=result.unavailable_a, unavailable_b=result.unavailable_b,
          excluded_a=result.excluded_a, excluded_b=result.excluded_b)
    for contact in result.contacts:
        pair = {**_residue_fields(contact.residue_a, "a"), **_residue_fields(contact.residue_b, "b")}
        write(record="pair", **pair, atom_a=contact.atom_a, atom_b=contact.atom_b,
              distance_angstrom=contact.distance_angstrom, atom_pair_count=contact.atom_pair_count)
        for interaction in contact.interactions:
            write(record="interaction", **pair, atom_a=interaction.atom_a, atom_b=interaction.atom_b,
                  distance_angstrom=interaction.distance_angstrom, interaction_type=interaction.interaction_type.value)
    for row in result.residues:
        side = "a" if row.residue_id.chain_id == result.chain_a else "b"
        write(record="residue", **_residue_fields(row.residue_id, side), partner_count=row.partner_count,
              atom_pair_count=row.atom_pair_count, distance_angstrom=row.min_distance_angstrom)
    atomic_write_bytes(path, output.getvalue().encode("utf-8-sig"), max_bytes=_MAX_CSV_BYTES,
                       label="interface CSV", replace_existing=True)


def open_subunit_interface(panel: Any) -> None:
    source = panel.reference_loaded_source
    if source is None or panel.reference_edit.text().strip() != str(source.path):
        panel._set_status("Load the reference complex before inspecting its subunit interfaces.")
        return
    previous = getattr(panel, "interface_window", None)
    if previous is not None:
        previous.dialog.close()
        previous.dialog.deleteLater()
    window = InterfaceWindow(panel, source)
    panel.interface_window = window
    window.dialog.show()
    window.dialog.raise_()


class InterfaceWindow:
    def __init__(self, panel: Any, source: LoadedSource) -> None:
        self.panel = panel
        self.source = source
        self.future: Future[SubunitInterface] | None = None
        self.cancel_event = Event()
        self.result: SubunitInterface | None = None
        self._displayed_pairs: tuple[InterfaceContact, ...] = ()
        w = panel.w
        self.dialog = w.QDialog(panel.widget)
        self.dialog.setWindowTitle("Subunit interfaces")
        self.dialog.setStyleSheet("QDialog, QTabWidget::pane, QTabWidget QWidget[pageCanvas=true] "
                                  "{ background: #111827; color: #dbe7f3; } QLabel { color: #dbe7f3; }"
                                  "QTabBar::tab { background: #19293e; color: #bfd2e8; padding: 6px 10px; }"
                                  "QTabBar::tab:selected { background: #24578a; color: white; }")
        self.dialog.setWindowModality(getattr(panel.c.Qt, "WindowModality", panel.c.Qt).WindowModal)
        self.dialog.resize(880, 550)
        outer = w.QVBoxLayout(self.dialog)
        plain = getattr(panel.c.Qt, "TextFormat", panel.c.Qt).PlainText
        header = w.QLabel(f"{source.path.name} · deposited coordinates · SHA-256 {source.snapshot.raw_sha256[:12]}")
        header.setTextFormat(plain)
        outer.addWidget(header)
        note = w.QLabel("Choose chains from a single model, including peptide chains. Contacts use heavy atoms in "
                        "the supplied coordinates. To study an assembly, load its coordinates. Partner counts "
                        "describe geometry and do not measure binding energy.")
        note.setWordWrap(True)
        outer.addWidget(note)
        settings = w.QHBoxLayout()
        self.model, self.chain_a, self.chain_b = w.QComboBox(), w.QComboBox(), w.QComboBox()
        self.cutoff = w.QDoubleSpinBox()
        self.cutoff.setRange(0.1, 20.0)
        self.cutoff.setDecimals(1)
        self.cutoff.setValue(5.0)
        self.cutoff.setSuffix(" Å")
        for title, control in (("Model", self.model), ("Chain A", self.chain_a), ("Chain B", self.chain_b),
                               ("Contact cutoff", self.cutoff)):
            label = w.QLabel(title)
            label.setBuddy(control)
            control.setAccessibleName(title)
            settings.addWidget(label)
            settings.addWidget(control)
        outer.addLayout(settings)
        self.tabs = w.QTabWidget()
        outer.addWidget(self.tabs, 1)
        self.pairs, self.pair_page = self._make_table(self.tabs, "Residue pairs", (
            "Residue A", "Residue B", "Closest atoms", "Min (Å)", "Atom pairs",
            "Geometric candidates · atoms A / B · distance (Å)"), self._fill_pairs)
        self.residues, self.residue_page = self._make_table(self.tabs, "Residues by partner count", (
            "Residue", "Distinct partners", "Atom pairs", "Minimum (Å)"), self._fill_residues)
        self.clear_filter_button = w.QPushButton("Show all residue pairs")
        self.clear_filter_button.setObjectName("secondaryButton")
        self.clear_filter_button.clicked.connect(lambda: self._display_pairs())
        self.tabs.widget(0).layout().addWidget(self.clear_filter_button)
        self.partners_button = w.QPushButton("Show selected residue’s partners")
        self.partners_button.setObjectName("secondaryButton")
        self.partners_button.clicked.connect(self._show_partners)
        self.tabs.widget(1).layout().addWidget(self.partners_button)
        self.residues.itemSelectionChanged.connect(lambda: self.partners_button.setEnabled(bool(self.residues.selectedItems())))
        self.residues.cellActivated.connect(lambda *_: self._show_partners())
        self.status = w.QLabel()
        self.status.setTextFormat(plain)
        self.status.setWordWrap(True)
        outer.addWidget(self.status)
        actions = w.QHBoxLayout()
        self.run_button, self.export_button = w.QPushButton("Calculate interface"), w.QPushButton("Export CSV…")
        self.run_button.setObjectName("primaryButton")
        self.export_button.setObjectName("secondaryButton")
        self.run_button.clicked.connect(self._run)
        self.export_button.clicked.connect(self._export)
        close = w.QPushButton("Close")
        close.setObjectName("secondaryButton")
        close.clicked.connect(self.dialog.close)
        for button in (self.run_button, self.export_button, close):
            actions.addWidget(button)
        outer.addLayout(actions)
        self.timer = panel.c.QTimer(self.dialog)
        self.timer.setInterval(40)
        self.timer.timeout.connect(self._poll)
        self.dialog.finished.connect(self._close)
        self.dialog.destroyed.connect(lambda *_args, event=self.cancel_event: event.set())
        for model in dict.fromkeys(str(chain.model_id) for chain in source.structure.chains):
            self.model.addItem(model, model)
        current = self.model.findData(panel.reference_model_combo.currentData())
        if current >= 0:
            self.model.setCurrentIndex(current)
        self.model.currentIndexChanged.connect(self._model_changed)
        self.chain_a.currentIndexChanged.connect(self._chain_changed)
        self.chain_b.currentIndexChanged.connect(self._reset)
        self.cutoff.valueChanged.connect(self._reset)
        self._model_changed()

    def _make_table(self, tabs: Any, title: str, headers: tuple[str, ...], on_change: Any) -> tuple[Any, Any]:
        w = self.panel.w
        page = w.QWidget()
        page.setProperty("pageCanvas", True)
        layout = w.QVBoxLayout(page)
        table = w.QTableWidget(0, len(headers))
        table.setAccessibleName(title)
        table.setHorizontalHeaderLabels(headers)
        _configure_table(table, w)
        table.resizeColumnsToContents()
        layout.addWidget(table)
        spin = w.QSpinBox()
        spin.setPrefix("Page ")
        spin.setRange(1, 1)
        spin.setAccessibleName(title + " page")
        spin.valueChanged.connect(on_change)
        layout.addWidget(spin)
        tabs.addTab(page, title)
        return table, spin

    def _chain_ids(self) -> tuple[str, ...]:
        return eligible_interface_chains(self.source.structure.chains, self.model.currentData())

    def _model_changed(self, *_: object) -> None:
        self.chain_a.blockSignals(True)
        self.chain_a.clear()
        for chain in self._chain_ids():
            self.chain_a.addItem(chain or "(blank)", chain)
        self.chain_a.blockSignals(False)
        self._chain_changed()

    def _chain_changed(self, *_: object) -> None:
        before = self.chain_b.currentData()
        self.chain_b.blockSignals(True)
        self.chain_b.clear()
        for chain in self._chain_ids():
            if chain != self.chain_a.currentData():
                self.chain_b.addItem(chain or "(blank)", chain)
        index = self.chain_b.findData(before)
        if index >= 0:
            self.chain_b.setCurrentIndex(index)
        self.chain_b.blockSignals(False)
        self._reset()

    def _reset(self, *_: object) -> None:
        self.result = None
        self._displayed_pairs = ()
        self.clear_filter_button.setEnabled(False)
        self.partners_button.setEnabled(False)
        for table, spin in ((self.pairs, self.pair_page), (self.residues, self.residue_page)):
            table.setRowCount(0)
            spin.setRange(1, 1)
        ready = self.chain_b.count() > 0
        self.export_button.setEnabled(False)
        self.run_button.setEnabled(ready)
        self.status.setText("Ready to calculate." if ready else
                            "Choose a model with two distinct, unambiguous protein or peptide chain IDs.")

    def _run(self) -> None:
        if self.future is not None or self.chain_b.currentData() is None:
            return
        self._reset()
        for control in (self.model, self.chain_a, self.chain_b, self.cutoff, self.run_button):
            control.setEnabled(False)
        self.status.setText("Calculating contacts…")
        self.future = _ANALYSIS_EXECUTOR.submit(analyze_subunit_interface, self.source.snapshot, self.model.currentData(),
                                                self.chain_a.currentData(), self.chain_b.currentData(), self.cutoff.value(),
                                                cancel_event=self.cancel_event)
        self.timer.start()

    def _poll(self) -> None:
        if self.future is None or not self.future.done():
            return
        self.timer.stop()
        finished, self.future = self.future, None
        for control in (self.model, self.chain_a, self.chain_b, self.cutoff, self.run_button):
            control.setEnabled(True)
        try:
            self.result = finished.result()
        except Exception as exc:
            self.status.setText(f"Interface analysis failed: {exc}")
            return
        r = self.result
        self.residue_page.setMaximum(max(1, (len(r.residues) + _PAGE_SIZE - 1) // _PAGE_SIZE))
        self._display_pairs()
        self._fill_residues()
        self.status.setText(f"{len(r.contacts)} contacting residue pairs at ≤ {r.cutoff_angstrom:g} Å. "
                            f"Residues with usable atoms: {r.chain_a} {r.observed_a}, {r.chain_b} {r.observed_b}. "
                            f"Without usable atoms: {r.chain_a} {r.unavailable_a}, {r.chain_b} {r.unavailable_b}. "
                            f"Excluded non-amino-acid residues: {r.chain_a} {r.excluded_a}, {r.chain_b} {r.excluded_b}.")
        self.export_button.setEnabled(True)

    def _show_partners(self) -> None:
        selected = self.residues.item(self.residues.currentRow(), 0)
        if selected is not None:
            residue = selected.data(user_role(self.panel.c))
            if isinstance(residue, ResidueId):
                self._display_pairs(residue)
                self.tabs.setCurrentIndex(0)

    def _display_pairs(self, residue: ResidueId | None = None) -> None:
        if self.result is None:
            return
        contacts = self.result.contacts
        self._displayed_pairs = contacts if residue is None else tuple(
            pair for pair in contacts if residue in (pair.residue_a, pair.residue_b)
        )
        self.pair_page.setMaximum(max(1, (len(self._displayed_pairs) + _PAGE_SIZE - 1) // _PAGE_SIZE))
        self.pair_page.setValue(1)
        self.clear_filter_button.setEnabled(residue is not None)
        self.clear_filter_button.setText(f"Showing partners of {_residue_label(residue)} · Show all pairs"
                                        if residue is not None else "Show all residue pairs")
        self._fill_pairs()

    def _fill_pairs(self, *_: object) -> None:
        if self.result is None:
            return
        start = (self.pair_page.value() - 1) * _PAGE_SIZE
        rows = self._displayed_pairs[start:start + _PAGE_SIZE]
        self.pairs.setRowCount(len(rows))
        for i, row in enumerate(rows):
            evidence = "\n".join(f"{item.interaction_type.value} · {item.atom_a} / {item.atom_b} · {item.distance_angstrom:.3f}"
                                 for item in row.interactions)
            self._write_row(self.pairs, i, (_residue_label(row.residue_a), _residue_label(row.residue_b),
                           f"{row.atom_a} / {row.atom_b}", f"{row.distance_angstrom:.3f}", row.atom_pair_count,
                           evidence or "Proximity only"))
        self.pairs.resizeRowsToContents()

    def _fill_residues(self, *_: object) -> None:
        if self.result is None:
            return
        start = (self.residue_page.value() - 1) * _PAGE_SIZE
        rows = self.result.residues[start:start + _PAGE_SIZE]
        self.residues.setRowCount(len(rows))
        for i, row in enumerate(rows):
            self._write_row(self.residues, i, (_residue_label(row.residue_id), row.partner_count, row.atom_pair_count,
                                               f"{row.min_distance_angstrom:.3f}"))
            self.residues.item(i, 0).setData(user_role(self.panel.c), row.residue_id)

    def _write_row(self, table: Any, index: int, values: tuple[Any, ...]) -> None:
        for column, value in enumerate(values):
            table.setItem(index, column, self.panel.w.QTableWidgetItem(str(value)))

    def _export(self) -> None:
        if self.result is None:
            return
        filename, _ = self.panel.w.QFileDialog.getSaveFileName(self.dialog, "Export subunit interface", "interface.csv",
                                                               "CSV (*.csv)")
        if filename:
            try:
                path = Path(filename)
                if path.resolve() == self.source.path.resolve():
                    raise ValueError("Choose a file other than the input structure.")
                export_interface_csv(path, self.result)
                self.status.setText(f"Exported {path.name}.")
            except (ValueError, OSError) as exc:
                self.status.setText(f"Export failed: {exc}")

    def _close(self, *_: object) -> None:
        self.timer.stop()
        self.cancel_event.set()
        if self.future is not None:
            self.future.cancel()
            self.future = None
