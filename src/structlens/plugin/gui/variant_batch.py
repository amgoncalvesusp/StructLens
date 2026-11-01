"""Window-modal reference-versus-variants queue backed by canonical pair reports."""

from __future__ import annotations

from collections import deque
from concurrent.futures import Future
from dataclasses import dataclass, replace
from html import escape
from pathlib import Path
from typing import Any, Literal

from structlens.application.dto import AnalysisReportRequest
from structlens.application.project_state import ProjectState
from structlens.application.report_service import ReportService
from structlens.core.models import AlignmentMode, ProteinChain

from .loaded_source import LoadedSource
from .model import CanonicalReportBinding
from .project_transaction import CanonicalProjectCandidate
from .qt_context import _ANALYSIS_EXECUTOR, _configure_table, _number
from .report_controller import AnalysisReportController


@dataclass(frozen=True, slots=True)
class Variant:
    path: Path
    source: LoadedSource | None = None
    binding: CanonicalReportBinding | None = None


def compare_variant(request: AnalysisReportRequest) -> CanonicalReportBinding:
    for path, snapshot in (
        (request.reference_selection.path, request.reference_snapshot),
        (request.target_selection.path, request.target_snapshot),
    ):
        if path is None:
            raise ValueError("Batch inputs require a source path.")
        AnalysisReportController.verify_current_source(path, snapshot)
    return CanonicalReportBinding.create(ReportService().analyze(request), request)


def open_variant_batch(panel: Any) -> None:
    if panel._future is not None or panel.model.busy:
        panel._set_status("Finish the current comparison before opening the variant batch.")
        return
    reference = panel.reference_loaded_source
    chain = panel._selected_chain(panel.reference_structure, panel.reference_chain_combo, panel.reference_model_combo)
    if reference is None or chain is None or panel.reference_edit.text().strip() != str(reference.path):
        panel._set_status("Load a reference file and select its model and chain before comparing variants.")
        return
    previous = panel.variant_batch
    key = (reference.path, reference.snapshot.raw_sha256, str(chain.model_id), chain.chain_id)
    if previous is None or previous.key != key:
        if previous is not None:
            previous.shutdown()
            previous.dialog.deleteLater()
        panel.variant_batch = VariantBatch(panel, reference, chain)
    panel.variant_batch.dialog.show()
    panel.variant_batch.dialog.raise_()


class VariantBatch:
    """Capture each file once, schedule one worker at a time, and retain completed bindings."""

    def __init__(self, panel: Any, reference: LoadedSource, chain: ProteinChain) -> None:
        self.panel = panel
        self.reference = reference
        self.chain = chain
        self.key = (reference.path, reference.snapshot.raw_sha256, str(chain.model_id), chain.chain_id)
        self.rows: list[Variant] = []
        self.jobs: deque[tuple[int, Literal["read", "compare"]]] = deque()
        self.active: tuple[int, Literal["read", "compare"]] | None = None
        self.future: Future[Any] | None = None
        self.cancelled = False
        w = panel.w
        self.dialog = w.QDialog(panel.widget)
        self.dialog.setWindowTitle("Compare variants against one reference")
        self.dialog.setWindowModality(getattr(panel.c.Qt, "WindowModality", panel.c.Qt).WindowModal)
        self.dialog.resize(820, 445)
        layout = w.QVBoxLayout(self.dialog)
        heading = w.QLabel(f"Pinned reference · {reference.path.name} · model {chain.model_id}, chain {chain.chain_id}")
        heading.setTextFormat(getattr(panel.c.Qt, "TextFormat", panel.c.Qt).PlainText)
        heading.setWordWrap(True)
        layout.addWidget(heading)
        note = w.QLabel("Add target files. Choose a model and chain when a file has several. "
                        "Click a completed pair to inspect or export its verified report. "
                        "Alignment settings and site definitions are pinned when the first target is added.")
        note.setWordWrap(True)
        layout.addWidget(note)
        self.table = w.QTableWidget(0, 7)
        self.table.setHorizontalHeaderLabels(
            ["Target file", "Model / chain", "Status", "Coverage", "Mapped", "Strict RMSD (Å)", "Changes"]
        )
        _configure_table(self.table, w)
        self.table.setColumnWidth(0, 155)
        self.table.setColumnWidth(1, 140)
        self.table.setColumnWidth(2, 160)
        self.table.cellClicked.connect(self.open_row)
        layout.addWidget(self.table, 1)
        actions = w.QHBoxLayout()
        self.add_button = w.QPushButton("Add target files…")
        self.run_button = w.QPushButton("Compare pending variants")
        self.cancel_button = w.QPushButton("Cancel queue")
        self.inspect_button = w.QPushButton("Open selected report")
        self.remove_button = w.QPushButton("Remove selected")
        for button in (self.add_button, self.remove_button, self.run_button, self.cancel_button, self.inspect_button):
            actions.addWidget(button)
        layout.addLayout(actions)
        self.notice = w.QLabel("No target files queued.")
        self.notice.setWordWrap(True)
        self.notice.setTextFormat(getattr(panel.c.Qt, "TextFormat", panel.c.Qt).PlainText)
        layout.addWidget(self.notice)
        self.settings = panel._settings()
        self.sites = panel._site_definitions
        self.manual_text = panel.manual_edit.toPlainText()
        self.add_button.clicked.connect(self.browse)
        self.remove_button.clicked.connect(self.remove_selected)
        self.run_button.clicked.connect(self.run)
        self.cancel_button.clicked.connect(self.cancel)
        self.inspect_button.clicked.connect(lambda: self.open_row(self.table.currentRow()))
        self.table.currentCellChanged.connect(lambda *_: self.refresh())
        self.dialog.finished.connect(lambda *_: self.cancel() if self.busy else None)
        self.timer = panel.c.QTimer(self.dialog)
        self.timer.setInterval(40)
        self.timer.timeout.connect(self.poll)
        self.refresh()

    @property
    def busy(self) -> bool:
        return self.active is not None or bool(self.jobs)

    def browse(self) -> None:
        paths, _ = self.panel.w.QFileDialog.getOpenFileNames(
            self.dialog, "Choose target structures", "", "Structures (*.pdb *.ent *.cif *.mmcif *.gz);;All files (*)"
        )
        self.add_paths([Path(path) for path in paths])

    def add_paths(self, paths: list[Path]) -> None:
        if self.busy:
            return
        if not self.rows:
            self.settings = self.panel._settings()
            self.sites = self.panel._site_definitions
            self.manual_text = self.panel.manual_edit.toPlainText()
        seen = {entry.path for entry in self.rows}
        for path in paths:
            path = path.absolute()
            if path in seen:
                continue
            seen.add(path)
            row = len(self.rows)
            self.rows.append(Variant(path))
            self.table.insertRow(row)
            self.cell(row, 0, path.name, str(path))
            self.cell(row, 2, "Queued")
            self.jobs.append((row, "read"))
        self.advance()

    def remove_selected(self) -> None:
        row = self.table.currentRow()
        if self.busy or row < 0:
            return
        self.rows.pop(row)
        self.table.removeRow(row)
        self.refresh()

    def cell(self, row: int, column: int, text: str, tip: str = "") -> None:
        item = self.panel.w.QTableWidgetItem(text)
        item.setToolTip("<qt>" + escape(tip or text).replace("\n", "<br>") + "</qt>")
        self.table.setItem(row, column, item)

    def scope(self, row: int) -> ProteinChain | None:
        chooser = self.table.cellWidget(row, 1)
        if chooser is None or chooser.currentData() is None:
            return None
        source = self.rows[row].source
        return source.structure.chains[int(chooser.currentData())] if source else None

    def loaded(self, row: int, source: LoadedSource) -> None:
        chains = source.structure.chains
        if not chains:
            raise ValueError("No selectable chains with coordinates.")
        self.rows[row] = replace(self.rows[row], source=source)
        chooser = self.panel.w.QComboBox()
        chooser.setAccessibleName(f"Model and chain for {source.path.name}")
        if len(chains) > 1:
            chooser.addItem("Choose model / chain", None)
        for index, chain in enumerate(chains):
            chooser.addItem(f"Model {chain.model_id} · {chain.chain_id}", index)
        self.table.setCellWidget(row, 1, chooser)
        chooser.currentIndexChanged.connect(lambda _index: self.scope_changed(chooser))
        self.scope_changed(chooser)

    def scope_changed(self, chooser: Any) -> None:
        row = next(i for i in range(len(self.rows)) if self.table.cellWidget(i, 1) is chooser)
        self.rows[row] = replace(self.rows[row], binding=None)
        self.cell(row, 2, "Ready" if self.scope(row) is not None else "Choose model / chain")
        for col in range(3, 7):
            self.cell(row, col, "—")
        self.refresh()

    def run(self) -> None:
        if self.busy:
            return
        self.jobs.extend((i, "compare") for i, entry in enumerate(self.rows)
                         if entry.binding is None and self.scope(i) is not None)
        self.advance()

    def request_for(self, row: int) -> AnalysisReportRequest:
        entry = self.rows[row]
        target_chain = self.scope(row)
        if entry.source is None or target_chain is None:
            raise ValueError("Select a model and chain for this target.")
        manual_pairs = (
            tuple(self.panel._manual_pairs(self.chain, target_chain, text=self.manual_text))
            if self.settings.alignment_mode is AlignmentMode.MANUAL else ()
        )
        return self.panel._build_snapshot_request(  # type: ignore[no-any-return]
            self.reference, entry.source, reference=self.chain, target=target_chain,
            analysis_settings=self.settings, site_definitions=self.sites, manual_pairs=manual_pairs,
        )

    def advance(self) -> None:
        while self.jobs:
            row, action = self.jobs.popleft()
            self.active = (row, action)
            try:
                if action == "read":
                    self.cell(row, 2, "Reading…")
                    self.future = _ANALYSIS_EXECUTOR.submit(LoadedSource.capture, self.rows[row].path)
                else:
                    self.cell(row, 2, "Comparing…")
                    self.future = _ANALYSIS_EXECUTOR.submit(compare_variant, self.request_for(row))
                self.timer.start()
                self.refresh()
                return
            except Exception as exc:
                self.cell(row, 2, f"Failed: {exc}")
        self.active = None
        self.future = None
        self.cancelled = False
        self.timer.stop()
        self.refresh()

    def poll(self) -> None:
        if self.active is None or self.future is None or not self.future.done():
            return
        row, action = self.active
        try:
            if self.cancelled:
                self.cell(row, 2, "Cancelled")
            elif action == "read":
                self.loaded(row, self.future.result())
            else:
                self.completed(row, self.future.result())
        except Exception as exc:
            self.cell(row, 2, f"Failed: {exc}")
        self.advance()

    def completed(self, row: int, binding: CanonicalReportBinding) -> None:
        self.rows[row] = replace(self.rows[row], binding=binding)
        analysis = binding.report.analysis
        if analysis is None:
            self.cell(row, 2, "No analysis · open diagnostics")
            return
        self.cell(row, 2, "Complete")
        values = (f"{analysis.sequence_coverage:.1%}", str(analysis.mapped_residue_count),
                  _number(analysis.strict_rmsd_angstrom), str(analysis.mutation_count))
        for col, value in enumerate(values, 3):
            self.cell(row, col, value)

    def cancel(self) -> None:
        self.cancelled = True
        while self.jobs:
            row, _ = self.jobs.popleft()
            self.cell(row, 2, "Cancelled")
        if self.active is not None and self.future is not None:
            self.cell(self.active[0], 2, "Cancelling…")
            self.future.cancel()
        self.refresh()

    def shutdown(self) -> None:
        self.cancelled = True
        self.jobs.clear()
        if self.future is not None:
            self.future.cancel()
        self.active = None
        self.future = None

    def refresh(self) -> None:
        busy = self.busy
        current = self.table.currentRow()
        self.add_button.setEnabled(not busy)
        self.remove_button.setEnabled(not busy and 0 <= current < len(self.rows))
        self.cancel_button.setEnabled(busy and not self.cancelled)
        ready = any(entry.binding is None and self.scope(i) is not None for i, entry in enumerate(self.rows))
        self.run_button.setEnabled(not busy and ready)
        self.inspect_button.setEnabled(not busy and 0 <= current < len(self.rows) and self.rows[current].binding is not None)
        for row in range(len(self.rows)):
            chooser = self.table.cellWidget(row, 1)
            if chooser is not None:
                chooser.setEnabled(not busy)
        complete = sum(entry.binding is not None for entry in self.rows)
        self.notice.setText(f"{complete} of {len(self.rows)} reports available · {self.settings.alignment_mode.value} alignment."
                            + (" Waiting for the current calculation to finish." if self.cancelled and busy else ""))

    def open_row(self, row: int, column: int = 0) -> None:
        if self.busy or row < 0 or column == 1:
            return
        entry = self.rows[row]
        if entry.binding is None or entry.source is None:
            return
        request = entry.binding.request
        project = ProjectState(settings=request.analysis_settings, comparison_mode=entry.binding.report.comparison_mode)
        candidate = CanonicalProjectCandidate(
            project, entry.binding, self.reference, entry.source, self.panel._state_from_controls()
        )
        try:
            self.panel._commit_canonical_project(candidate, None)
        except Exception as exc:
            self.notice.setText(f"Could not open verified report: {exc}")
            return
        self.dialog.accept()
