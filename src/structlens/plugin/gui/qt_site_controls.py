"""Editable site requests for the selected reference snapshot."""

from __future__ import annotations

import re
from typing import Any

from structlens.application.report_evidence import reference_ligand_atoms
from structlens.application.report_input import canonical_selection, canonical_site_definition, canonical_snapshot
from structlens.application.site_service import define_site
from structlens.core.parsing import ParsedStructure, load_structure_evidence
from structlens.core.sites import SiteDefinition, SiteDefinitionMode

from .qt_context import SCIENTIFIC_SECTIONS, QtMixinContext, _button, _find_residue, _label, _residue_label


class SiteControlsMixin(QtMixinContext):
    _site_reference: ParsedStructure | None = None
    _site_reference_key: tuple[str, str, str] | None = None

    def _build_site_controls(self, content: Any) -> None:
        group = self.w.QGroupBox("Site definition", self.widget)
        grid = self.w.QGridLayout(group)
        self.site_name_edit = self.w.QLineEdit(group)
        self.site_name_edit.setMaxLength(80)
        self.site_name_edit.setPlaceholderText("e.g. ATP pocket")
        self.site_mode_combo = self.w.QComboBox(group)
        for title, mode in (("Key residues", SiteDefinitionMode.KEY_RESIDUES),
                            ("Ligand radius", SiteDefinitionMode.LIGAND_RADIUS),
                            ("Residue radius", SiteDefinitionMode.RESIDUE_RADIUS)):
            self.site_mode_combo.addItem(title, mode.value)
        self.site_ligand_combo = self.w.QComboBox(group)
        self.site_ligand_combo.addItem("Load a reference structure", None)
        self.site_residues_edit = self.w.QLineEdit(group)
        self.site_residues_edit.setPlaceholderText("A:70, A:73, A:166")
        self.site_radius_spin = self.w.QDoubleSpinBox(group)
        self.site_radius_spin.setRange(0.1, 20.0)
        self.site_radius_spin.setValue(5.0)
        self.site_radius_spin.setSuffix(" Å")
        fields = (("Site name", self.site_name_edit), ("Definition", self.site_mode_combo),
                  ("Reference positions", self.site_residues_edit), ("Radius", self.site_radius_spin))
        for index, (label, field) in enumerate(fields):
            row, col = divmod(index, 2)
            grid.addWidget(_label(self.w, label, "fieldLabel"), row * 2, col)
            grid.addWidget(field, row * 2 + 1, col)
            field.setAccessibleName(label)
        grid.addWidget(_label(self.w, "Reference ligand", "fieldLabel"), 4, 0, 1, 2)
        grid.addWidget(self.site_ligand_combo, 5, 0, 1, 2)
        self.site_ligand_combo.setAccessibleName("Reference ligand")
        self.site_preview_label = _label(self.w, "", "inlineNote")
        self.site_preview_label.setTextFormat(getattr(self.c.Qt, "TextFormat", self.c.Qt).PlainText)
        grid.addWidget(self.site_preview_label, 6, 0, 1, 2)
        self.site_define_button = _button(self.w, "Add site to analysis", "secondaryButton")
        self.site_define_button.clicked.connect(self._define_site_from_controls)
        grid.addWidget(self.site_define_button, 7, 0, 1, 2)
        self.configured_site_combo = self.w.QComboBox(group)
        self.configured_site_combo.setAccessibleName("Configured sites")
        self.site_remove_button = _button(self.w, "Remove selected site", "secondaryButton")
        self.site_remove_button.clicked.connect(self._remove_selected_site)
        grid.addWidget(self.configured_site_combo, 8, 0)
        grid.addWidget(self.site_remove_button, 8, 1)
        self.site_status_label = _label(self.w, "Add a site to calculate its metrics in the next comparison.", "inlineNote")
        grid.addWidget(self.site_status_label, 9, 0, 1, 2)
        self.site_ligand_combo.currentIndexChanged.connect(self._refresh_site_preview)
        self.site_mode_combo.currentIndexChanged.connect(self._refresh_site_preview)
        self.site_residues_edit.textChanged.connect(self._refresh_site_preview)
        self.site_radius_spin.valueChanged.connect(self._refresh_site_preview)
        content.addWidget(group)
        self._refresh_configured_sites()
        self._refresh_site_preview()

    def _refresh_site_reference(self) -> None:
        source = self.reference_loaded_source
        model_id = self._combo_data(self.reference_model_combo)
        chain_id = self._combo_data(self.reference_chain_combo)
        key = (source.snapshot.content_id, model_id, chain_id) if source and model_id and chain_id else None
        if key == self._site_reference_key:
            return
        self._site_reference_key = key
        self._site_reference = None
        self._site_definitions = ()
        self._refresh_configured_sites()
        message = "No ligand in this reference selection"
        if key and source:
            selection = source.selection(model_id=key[1], chain_id=key[2])
            try:
                self._site_reference = load_structure_evidence(
                    canonical_snapshot(source.snapshot, "reference"), selection=canonical_selection(selection, "reference")
                )
            except (OSError, ValueError, RuntimeError) as error:
                message = f"Ligand information unavailable: {error}"
        combo = self.site_ligand_combo
        previous = combo.blockSignals(True)
        try:
            combo.clear()
            ligands = [item for item in self._site_reference.components if item.is_ligand] if self._site_reference else []
            if not ligands:
                combo.addItem(message if key else "Load a reference structure", None)
            for item in ligands:
                position = f"{item.author_chain_id}:{item.auth_seq_id}{item.insertion_code or ''}"
                label = f"{item.residue_name} · {position} · model {item.model_id}"
                if item.label_chain_id and item.label_chain_id != item.author_chain_id:
                    label += f" · label chain {item.label_chain_id}"
                if item.entity_id:
                    label += f" · entity {item.entity_id}"
                combo.addItem(label, item.component_id)
                combo.setItemData(combo.count() - 1, item.component_id, getattr(self.c.Qt, "ItemDataRole", self.c.Qt).ToolTipRole)
        finally:
            combo.blockSignals(previous)
        if not previous:
            combo.currentIndexChanged.emit(combo.currentIndex())
        else:
            self._refresh_site_preview()

    def _draft_site(self) -> SiteDefinition:
        mode = SiteDefinitionMode(self.site_mode_combo.currentData())
        reference = self._selected_chain(self.reference_structure, self.reference_chain_combo, self.reference_model_combo)
        if reference is None:
            raise ValueError("Load a reference chain before defining a site.")
        name = self.site_name_edit.text().strip()
        if not name:
            number = 1
            while any(item.name == f"Site {number}" for item in self._site_definitions):
                number += 1
            name = f"Site {number}"
        slug = re.sub(r"[^\w-]+", "-", name.casefold()).strip("-")
        if not slug:
            raise ValueError("Use at least one letter or number in the site name.")
        site_id = f"site-{slug}"
        radius = float(self.site_radius_spin.value())
        if mode is SiteDefinitionMode.LIGAND_RADIUS:
            ligand_id = self.site_ligand_combo.currentData()
            if not ligand_id:
                raise ValueError("Choose a deposited ligand identifier in the reference selection.")
            return SiteDefinition(site_id, name, mode, ligand_id=ligand_id, radius_angstrom=radius)
        positions = tuple(part.strip() for part in self.site_residues_edit.text().split(",") if part.strip())
        if mode is SiteDefinitionMode.RESIDUE_RADIUS and len(positions) != 1:
            raise ValueError("Residue-radius sites require one reference center position, for example A:166.")
        if not positions:
            raise ValueError("Key-residue sites need at least one reference position.")
        residues = tuple(_find_residue(reference, position) for position in positions)
        if any(item is None for item in residues):
            raise ValueError("One or more key positions are not present in the selected reference chain.")
        valid = tuple(dict.fromkeys(item for item in residues if item is not None))
        if mode is SiteDefinitionMode.KEY_RESIDUES:
            return SiteDefinition(site_id, name, mode, reference_residues=valid)
        return SiteDefinition(site_id, name, mode, center_residue=valid[0], radius_angstrom=radius)

    def _refresh_site_preview(self, *_: object) -> None:
        mode = self.site_mode_combo.currentData()
        self.site_residues_edit.setEnabled(mode != SiteDefinitionMode.LIGAND_RADIUS.value)
        self.site_radius_spin.setEnabled(mode != SiteDefinitionMode.KEY_RESIDUES.value)
        try:
            draft = self._draft_site()
            chain = self._selected_chain(self.reference_structure, self.reference_chain_combo, self.reference_model_combo)
            atoms = reference_ligand_atoms(self._site_reference) if self._site_reference else {}
            if self._site_reference:
                if len(self._site_reference.protein_structure.chains) != 1:
                    raise ValueError("Select one reference chain to preview its site.")
                draft = canonical_site_definition(draft)
                chain = self._site_reference.protein_structure.chains[0]
            members = define_site(draft, chain.residue_records, ligand_atoms=atoms)
            labels = ", ".join(_residue_label(item.residue_id) for item in members[:25])
            remaining = f" (+{len(members) - 25} more)" if len(members) > 25 else ""
            radius = f" within {draft.radius_angstrom:g} Å" if draft.radius_angstrom is not None else ""
            text = f"{len(members)} reference residue(s){radius}: {labels or 'none'}{remaining}"
        except (ValueError, TypeError) as error:
            text = str(error)
        self.site_preview_label.setText(text)

    def _define_site_from_controls(self) -> None:
        try:
            draft = self._draft_site()
            if any(item.site_id == draft.site_id or item.name.casefold() == draft.name.casefold()
                   for item in self._site_definitions):
                raise ValueError("This site name is already in use. Choose another name or remove the existing site.")
        except (TypeError, ValueError) as error:
            self._site_unavailable(str(error))
            return
        self._site_definitions = (*self._site_definitions, draft)
        self._refresh_configured_sites()
        self.configured_site_combo.setCurrentIndex(self.configured_site_combo.findData(draft.site_id))
        self.nav.setCurrentRow(SCIENTIFIC_SECTIONS.index("Sites"))
        self.site_status_label.setText("Site added. Run Compare structures to update the results.")
        self._refresh_workflow_state()

    def _remove_selected_site(self) -> None:
        site_id = self.configured_site_combo.currentData()
        if site_id is None:
            return
        self._site_definitions = tuple(site for site in self._site_definitions if site.site_id != site_id)
        self._refresh_configured_sites()
        self.site_status_label.setText("Site removed. Run Compare structures to update the results.")
        self._refresh_workflow_state()

    def _refresh_configured_sites(self) -> None:
        self.configured_site_combo.clear()
        for item in self._site_definitions:
            self.configured_site_combo.addItem(item.name, item.site_id)
        self.configured_site_combo.setEnabled(bool(self._site_definitions))
        self.site_remove_button.setEnabled(bool(self._site_definitions))
