"""Ligand-centred 2D contact diagram and its linked evidence table."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

from structlens.application.ligand_contacts import (
    PROXIMITY_ANGSTROM,
    Contact,
    DockLensContacts,
    import_docklens_contacts,
    ligand_proximity,
    residue_tag,
)
from structlens.core.parsing import ParsedStructure

from .qt_context import QtMixinContext, _button, _configure_table, _label, _residue_label

_PER_PAGE = 8
_PROFILE_LABELS = {
    "plip": "PLIP-like",
    "luna": "LUNA-like",
    "dsv": "DSV-like",
    "luna_dsv": "LUNA-like × DSV-like (conservative)",
}


@dataclass(frozen=True, slots=True)
class ContactView:
    report_id: str
    receptor_hash: str
    ligand_id: str | None
    imported: DockLensContacts | None
    pose_id: str | None
    row: int
    page: int


class LigandContactsMixin(QtMixinContext):
    def _build_ligand_contacts_panel(self, content: Any) -> None:
        group = self.w.QGroupBox("Ligand contact map", self.widget)
        layout = self.w.QVBoxLayout(group)
        self._contact_rows: tuple[Contact, ...] = ()
        self._contact_indices: tuple[int | None, ...] = ()
        self._contact_nodes: dict[int, Any] = {}
        self._contact_page = 0
        self._contact_report_id: str | None = None
        self._docklens_contacts: DockLensContacts | None = None
        self._contact_center = "Ligand"

        controls = self.w.QHBoxLayout()
        self.contact_source = self.w.QComboBox()
        self.contact_source.setAccessibleName("Contact source or external DockLens pose")
        self.contact_source.addItem("Deposited ligand · reference coordinates", None)
        self.contact_source.currentIndexChanged.connect(lambda _: self._refresh_ligand_contacts())
        controls.addWidget(self.contact_source, 1)
        self.contact_import_button = _button(self.w, "Import DockLens result…", "secondaryButton")
        self.contact_import_button.setToolTip(
            "Import docklens-result-v1 for the completed single-model PDB receptor and verify the original poses file."
        )
        self.contact_import_button.clicked.connect(self._import_ligand_contacts)
        controls.addWidget(self.contact_import_button)
        layout.addLayout(controls)

        self.ligand_contact_status = _label(self.w, "Compare structures to inspect ligand contacts.", "helpText")
        plain = getattr(self.c.Qt, "TextFormat", self.c.Qt).PlainText
        self.ligand_contact_status.setTextFormat(plain)
        layout.addWidget(self.ligand_contact_status)
        self.ligand_contact_scene = self.w.QGraphicsScene(group)
        self.ligand_contact_scene.setSceneRect(0, 0, 510, 185)
        self.ligand_contact_scene.setBackgroundBrush(self.g.QColor("#142133"))
        self.ligand_contact_view = self.w.QGraphicsView(self.ligand_contact_scene)
        self.ligand_contact_view.setAccessibleName("Contact map; equivalent contacts can be selected in the table")
        self.ligand_contact_view.setMinimumHeight(212)
        layout.addWidget(self.ligand_contact_view)
        self.ligand_contact_scene.selectionChanged.connect(self._contact_node_clicked)

        paging = self.w.QHBoxLayout()
        self.contact_prev = _button(self.w, "Previous", "secondaryButton")
        self.contact_next = _button(self.w, "Next", "secondaryButton")
        self.contact_prev.clicked.connect(lambda: self._draw_contact_map(self._contact_page - 1))
        self.contact_next.clicked.connect(lambda: self._draw_contact_map(self._contact_page + 1))
        self.contact_page_label = _label(self.w, "", "inlineNote")
        for widget in (self.contact_prev, self.contact_page_label, self.contact_next):
            paging.addWidget(widget)
        paging.addStretch()
        layout.addLayout(paging)
        self.ligand_contact_table = self.w.QTableWidget(0, 6)
        self.ligand_contact_table.setHorizontalHeaderLabels(
            ["Residue", "Interaction", "Receptor atom", "Ligand atom", "Distance (Å)", "Evidence"]
        )
        _configure_table(self.ligand_contact_table, self.w)
        self.ligand_contact_table.currentCellChanged.connect(lambda row, *_: self._contact_row_changed(row))
        layout.addWidget(self.ligand_contact_table)
        self.contact_evidence_button = _button(self.w, "Open residue evidence", "secondaryButton")
        self.contact_evidence_button.clicked.connect(self._open_contact_evidence)
        layout.addWidget(self.contact_evidence_button)
        content.addWidget(group)
        self.site_ligand_combo.currentIndexChanged.connect(lambda _: self._refresh_ligand_contacts())
        self._refresh_ligand_contacts()

    def _refresh_ligand_contacts(self) -> None:
        if not hasattr(self, "contact_source"):
            return
        binding = self.model.binding
        report_id = binding.report.report_id if binding else None
        if self._contact_report_id != report_id:
            self._contact_report_id = report_id
            self._docklens_contacts = None
            blocked = self.contact_source.blockSignals(True)
            while self.contact_source.count() > 1:
                self.contact_source.removeItem(1)
            self.contact_source.blockSignals(blocked)
        reference: ParsedStructure | None = getattr(self, "_site_reference", None)
        ready = bool(
            binding
            and reference
            and binding.report.analysis
            and binding.request.reference_selection.selection_id == reference.selection.selection_id
        )
        self.contact_import_button.setEnabled(
            ready
            and reference is not None
            and reference.selection.format.value == "pdb"
            and len(reference.metadata.available_models) == 1
        )
        if not ready or reference is None or binding is None or binding.report.analysis is None:
            self._set_contact_rows(
                (), (), "Compare the selected reference model and chain to inspect contacts.", "Ligand"
            )
            return
        analysis = binding.report.analysis
        source = self.contact_source.currentData()
        if source is None:
            ligand_id = self.site_ligand_combo.currentData()
            try:
                neighborhood = ligand_proximity(reference, str(ligand_id or ""))
            except ValueError as error:
                self._set_contact_rows((), (), str(error), "Ligand")
                return
            rows = neighborhood.contacts
            eligible = {c.reference: c.alignment_index for c in analysis.correspondences if c.reference is not None}
            indices = tuple(eligible.get(row.residue_id) if row.residue_id is not None else None for row in rows)
            status = (
                f"StructLens · nearest observed heavy atoms ≤ {PROXIMITY_ANGSTROM:g} Å · {len(rows)} residues. "
                f"{neighborhood.outside} residues outside this cutoff among measured atoms; "
                f"{neighborhood.unavailable} without observed heavy atoms. Proximity alone does not identify a chemical interaction."
            )
            center = str(ligand_id)
        else:
            imported = self._docklens_contacts
            if imported is None:
                return
            pose = next(p for p in imported.poses if p.pose_id == source)
            rows = pose.rows
            candidates: dict[str, list[int | None]] = defaultdict(list)
            for c in analysis.correspondences:
                if c.reference:
                    # DockLens PDB parser drops insertion codes; even a unique inserted position is not proven.
                    candidates[residue_tag(c.reference, insertion=False)].append(
                        None if c.reference.insertion_code else c.alignment_index
                    )
            indices = tuple(candidates[r.residue][0] if len(candidates[r.residue]) == 1 else None for r in rows)
            center = pose.ligand_id
            status = (
                f"DockLens {imported.version} · {_PROFILE_LABELS[imported.profile]} profile · counting unit: {imported.counting_unit}. "
                f"External docked pose {pose.pose_id} · {len(rows)} reported contacts · run {imported.run_id} "
                f"({imported.errors} input error(s)) · {indices.count(None)} unmapped in the selected chain. "
                "This pose is separate from the deposited ligand."
            )
        self._set_contact_rows(rows, indices, status, center)

    def _set_contact_rows(
        self, rows: tuple[Contact, ...], indices: tuple[int | None, ...], status: str, center: str
    ) -> None:
        self._contact_rows, self._contact_indices, self._contact_center = rows, indices, center
        self.ligand_contact_status.setText(status)
        self._draw_contact_map(0)
        if self.model.selection is not None and self.model.report and self.model.report.analysis:
            selected = self.model.selection.reference_residue
            for c in self.model.report.analysis.correspondences:
                if selected is not None and selected == c.reference:
                    self._sync_ligand_contact_selection(c.alignment_index)
                    break

    def _fill_contact_table(self, rows: tuple[Contact, ...], selected: int) -> None:
        table = self.ligand_contact_table
        blocked = table.blockSignals(True)
        table.setRowCount(0)
        table.setRowCount(len(rows))
        for i, row in enumerate(rows):
            label = _residue_label(row.residue_id) if row.residue_id else row.residue
            for col, value in enumerate(
                (
                    label,
                    row.family,
                    row.receptor_atom,
                    row.ligand_atom,
                    "Unavailable" if row.distance is None else f"{row.distance:.3f}",
                    row.evidence,
                )
            ):
                cell = self.w.QTableWidgetItem(value)
                table.setItem(i, col, cell)
        table.setCurrentCell(selected, 0 if selected >= 0 else -1)
        table.blockSignals(blocked)

    def _draw_contact_map(self, page: int, selected: int = -1) -> None:
        scene = self.ligand_contact_scene
        table = self.ligand_contact_table
        total = max(1, (len(self._contact_rows) + _PER_PAGE - 1) // _PER_PAGE)
        self._contact_page = max(0, min(page, total - 1))
        start = self._contact_page * _PER_PAGE
        page_rows = self._contact_rows[start : start + _PER_PAGE]
        self._fill_contact_table(page_rows, selected)
        self.contact_evidence_button.setEnabled(
            0 <= selected < len(page_rows) and self._contact_indices[start + selected] is not None
        )
        self.contact_page_label.setText(
            f"Page {self._contact_page + 1} of {total} · {len(self._contact_rows)} contacts"
        )
        self.contact_prev.setEnabled(self._contact_page > 0)
        self.contact_next.setEnabled(self._contact_page + 1 < total)
        blocked = scene.blockSignals(True)
        try:
            scene.clear()
            self._contact_nodes = {}
            if not self._contact_rows:
                text = scene.addSimpleText("No contact records to display")
                text.setBrush(self.g.QColor("#dce6f3"))
                text.setPos(12, 12)
                return
            self._contact_map_box(193, 74, 124, self._contact_center, "Ligand", None)
            for row, contact in enumerate(page_rows):
                left = row % 2 == 0
                x, y = (4 if left else 362), (row // 2) * 45 + 3
                line = scene.addLine(
                    193 if left else 317, 91, 148 if left else 362, y + 17, self.g.QPen(self.g.QColor("#82a8bb"))
                )
                line.setZValue(-1)
                node = self._contact_map_box(x, y, 144, table.item(row, 0).text(), contact.family, row)
                self._contact_nodes[row] = node
                node.setSelected(row == table.currentRow())
        finally:
            scene.blockSignals(blocked)

    def _contact_map_box(self, x: int, y: int, width: int, title: str, subtitle: str, row: int | None) -> Any:
        scene = self.ligand_contact_scene
        rect = scene.addRect(
            0,
            0,
            width,
            36,
            self.g.QPen(self.g.QColor("#7fa4bf")),
            self.g.QBrush(self.g.QColor("#184e62" if row is None else "#253a50")),
        )
        rect.setPos(x, y)
        for i, label in enumerate((title, subtitle)):
            text = self.w.QGraphicsSimpleTextItem(label, rect)
            font = text.font()
            font.setPointSize(8)
            text.setFont(font)
            metrics = self.g.QFontMetrics(font)
            mode = getattr(self.c.Qt, "TextElideMode", self.c.Qt).ElideRight
            text.setText(metrics.elidedText(label, mode, width - 10))
            text.setBrush(self.g.QColor("#e7f1f5"))
            text.setPos(5, 1 + 16 * i)
        if row is not None:
            flags = getattr(self.w.QGraphicsItem, "GraphicsItemFlag", self.w.QGraphicsItem)
            rect.setFlag(flags.ItemIsSelectable, True)
            rect.setData(0, row)
        return rect

    def _contact_node_clicked(self) -> None:
        chosen = self.ligand_contact_scene.selectedItems()
        if chosen:
            self.ligand_contact_table.setCurrentCell(int(chosen[0].data(0)), 0)

    def _contact_row_changed(self, row: int) -> None:
        if not 0 <= row < self.ligand_contact_table.rowCount():
            self.contact_evidence_button.setEnabled(False)
            return
        scene = self.ligand_contact_scene
        blocked = scene.blockSignals(True)
        scene.clearSelection()
        self._contact_nodes[row].setSelected(True)
        scene.blockSignals(blocked)
        index = self._contact_indices[self._contact_page * _PER_PAGE + row]
        self.contact_evidence_button.setEnabled(index is not None)
        if index is not None:
            self._select_linked_residue(index)
        else:
            self.residue_table.setCurrentCell(-1, -1)
            self._selection_controller = self._selection_controller.clear_selection()
            self.model = replace(self.model, selection=None)
            self._sync_alignment_selection(None)
            self._highlight_chart_residue(None)
            self.evidence_card_label.setText("This DockLens endpoint does not map uniquely to the selected chain.")
            self._refresh_workflow_state()

    def _sync_ligand_contact_selection(self, alignment_index: int | None) -> None:
        if not hasattr(self, "ligand_contact_table"):
            return
        table = self.ligand_contact_table
        current = table.currentRow()
        if (
            alignment_index is not None
            and current >= 0
            and self._contact_indices[self._contact_page * _PER_PAGE + current] == alignment_index
        ):
            return
        blocked = table.blockSignals(True)
        if alignment_index is not None and alignment_index in self._contact_indices:
            row = self._contact_indices.index(alignment_index)
            page, offset = divmod(row, _PER_PAGE)
            self._draw_contact_map(page, offset)
        else:
            table.setCurrentCell(-1, -1)
            self.ligand_contact_scene.clearSelection()
            self.contact_evidence_button.setEnabled(False)
        table.blockSignals(blocked)

    def _open_contact_evidence(self) -> None:
        row = self.ligand_contact_table.currentRow()
        if row >= 0 and (index := self._contact_indices[self._contact_page * _PER_PAGE + row]) is not None:
            self._open_residue_filter("all", index)

    def _import_ligand_contacts(self) -> None:
        binding = self.model.binding
        if binding is None or not self.contact_import_button.isEnabled():
            return
        result, _ = self.w.QFileDialog.getOpenFileName(self.widget, "DockLens paired result", "", "JSON (*.json)")
        if not result:
            return
        poses, _ = self.w.QFileDialog.getOpenFileName(self.widget, "Original poses file used by DockLens")
        if not poses:
            return
        try:
            imported = import_docklens_contacts(
                Path(result), binding.request.reference_snapshot.raw_sha256, Path(poses)
            )
        except (ValueError, OSError) as exc:
            self.ligand_contact_status.setText(f"DockLens import failed: {exc}")
            return
        self._set_contact_source(imported, imported.poses[0].pose_id if imported.poses else None)
        self._set_status(f"Verified DockLens receptor and poses SHA-256 · {len(imported.poses)} pose(s).")

    def _set_contact_source(self, imported: DockLensContacts | None, pose_id: str | None) -> None:
        self._docklens_contacts = imported
        combo = self.contact_source
        blocked = combo.blockSignals(True)
        try:
            while combo.count() > 1:
                combo.removeItem(1)
            for pose in imported.poses if imported else ():
                combo.addItem(f"DockLens · {pose.ligand_id} / {pose.pose_id}", pose.pose_id)
            combo.setCurrentIndex(max(0, combo.findData(pose_id)))
        finally:
            combo.blockSignals(blocked)
        self._refresh_ligand_contacts()

    def _capture_contact_view(self) -> ContactView | None:
        binding = self.model.binding
        if not binding:
            return None
        return ContactView(
            binding.report.report_id,
            binding.request.reference_snapshot.raw_sha256,
            self.site_ligand_combo.currentData(),
            self._docklens_contacts,
            self.contact_source.currentData(),
            self.ligand_contact_table.currentRow(),
            self._contact_page,
        )

    def _restore_contact_view(self, view: ContactView | None) -> None:
        binding = self.model.binding
        if (
            not binding
            or not view
            or view.report_id != binding.report.report_id
            or (view.receptor_hash != binding.request.reference_snapshot.raw_sha256)
        ):
            self._set_contact_source(None, None)
            return
        self._contact_report_id = view.report_id
        combo = self.site_ligand_combo
        if combo.findData(view.ligand_id) >= 0:
            combo.setCurrentIndex(combo.findData(view.ligand_id))
        self._set_contact_source(view.imported, view.pose_id)
        self._draw_contact_map(view.page, view.row)
