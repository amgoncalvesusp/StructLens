"""Focused reports presentation/routing mixin for the Qt panel."""

# ruff: noqa: F403,F405

from __future__ import annotations

from .presentation import present_report
from .qt_context import *  # noqa: F401,F403


class ReportMixin(QtMixinContext):
    """Cohesive GUI-only reports behavior composed into PanelController."""

    def _invalidate_analysis_views(self) -> None:
        """Drop derived v0.3 views whenever their authoritative analysis changes."""

        self.model = replace(self.model, report=None, binding=None, analysis=None, selection=None)
        self._report_presentation = None
        self._selection_controller = AnalysisSelectionController()
        self._v03_export_records = {}
        self._chart_datasets = {}
        self.result_decision.setText("No comparison yet.")
        for value in self.result_labels.values():
            value.setText("—")
        self._render_msa_result(None)
        self._v03_bundle_payloads = {
            "msa_summary": None,
            "conservation": None,
            "interactions": None,
            "sites": None,
            "evidence": None,
            "vectors": None,
        }
        self.mutation_table.setRowCount(0)
        self.residue_table.setRowCount(0)
        self.structure_result_table.setRowCount(0)
        self.structure_result_summary.setText(
            "No structure comparison result yet. Run comparison to populate this tab."
        )
        self._render_site_metrics(())
        self.chart_preview_status.setText("Chart unavailable. Run the corresponding scientific service first.")
        self._clear_chart_layout(self.chart_preview_layout)
        self._clear_chart_layout(self.sequence_chart_layout)
        self._update_chart_export_state(self.chart_combo.currentText())
        if hasattr(self, "quality_table"):
            self.quality_table.setRowCount(0)
            self.quality_summary_label.setText("No QC report loaded.")
        if hasattr(self, "pocket_table"):
            self.pocket_table.setRowCount(0)
            self.pocket_status_label.setText("Pocket evidence appears after Compare.")
            self.pocket_capability_label.setText("No canonical pocket report loaded.")
            self.pocket_detail_label.setText("Select a pocket candidate to inspect volume evidence.")
            self.pocket_detect_button.setEnabled(False)
            self.pocket_measure_button.setEnabled(False)
            self._pocket_presentation = None

    def _populate_result(self, result: AnalysisResult) -> None:
        self.result_decision.setText(
            f"<b>{result.alignment_decision}</b><br>Reference: {result.reference_id} · Target: {result.target_id}"
        )
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
        self.residue_summary.setText(
            f"{len(result.correspondences)} aligned positions · {result.mapped_residue_count} mapped Cα pairs"
        )
        self._fill_structure_result(result)
        self._fill_results_history()
        self._render_sequence_result(result)
        self._render_selected_chart()
        self._fill_mutations(result)
        self._fill_residues(result)
        self._update_legend()

    def _populate_report(
        self,
        report: AnalysisReport,
        *,
        request: AnalysisReportRequest | None = None,
    ) -> None:
        """Validate, then render every available view from one report."""

        # Candidate construction is deliberately before every model, history,
        # selection, presentation and widget mutation.
        presentation = present_report(report)
        binding = None if request is None else CanonicalReportBinding.create(report, request, presentation)
        self._commit_report(report, presentation, binding)

    def _commit_report(
        self,
        report: AnalysisReport,
        presentation: ReportPresentation,
        binding: CanonicalReportBinding | None,
    ) -> None:
        """Commit a prevalidated report projection to the panel."""

        self._invalidate_analysis_views()
        self._report_history = (*self._report_history, report)
        self.model = self.model.with_report(report) if binding is None else self.model.with_binding(binding)
        self._selection_controller = AnalysisSelectionController().bind_report(report)
        self._report_presentation = presentation
        self._v03_bundle_payloads = {
            name: value if value is not None else None for name, value in presentation.bundle_payloads.items()
        }
        self._render_report_presentation(presentation)
        self._update_legend()
        if binding is None:
            self._report_controller.accept_report(report)
        else:
            self._report_controller.accept_binding(binding)
            self._report_request = binding.request
        self._pending_report_request = None
        self._set_busy(False)
        self._set_status(f"Analysis report complete · {report.report_id}")
        self.nav.setCurrentRow(_STRUCTURES_PAGE_INDEX)

    def _fill_report_structure_result(self, report: AnalysisReport) -> None:
        analysis = report.analysis
        if analysis is None:
            return
        values = (
            analysis.reference_id,
            analysis.target_id,
            analysis.alignment_decision,
            _number(analysis.strict_rmsd_angstrom),
            _number(analysis.refined_rmsd_angstrom),
            _number(analysis.tm_score),
            str(analysis.mapped_residue_count),
            str(len(analysis.excluded_alignment_indices)),
            str(analysis.method_provenance.method_id if analysis.method_provenance else "report"),
        )
        self.structure_result_table.setRowCount(1)
        for column, value in enumerate(values):
            self.structure_result_table.setItem(0, column, self.w.QTableWidgetItem(value))
        self.structure_result_summary.setText(
            f"{analysis.reference_id} → {analysis.target_id} · decision={analysis.alignment_decision} · "
            f"{analysis.mapped_residue_count} mapped residue pair(s); excluded positions remain explicit."
        )

    def _fill_report_mutations(self, report: AnalysisReport) -> None:
        values = report.analysis.mutations if report.analysis is not None else ()
        self.mutation_table.setRowCount(0)
        for event in values:
            row = self.mutation_table.rowCount()
            self.mutation_table.insertRow(row)
            cells = (
                str(event.alignment_index),
                event.kind.value,
                _residue_label(event.reference),
                _residue_label(event.target),
                event.canonical_notation,
                _number(event.blosum62_score),
                _number(event.grantham_distance),
                event.physicochemical_class or "—",
            )
            for column, value in enumerate(cells):
                self.mutation_table.setItem(row, column, self.w.QTableWidgetItem(value))

    def _fill_report_residues(self, report: AnalysisReport) -> None:
        values = report.analysis.correspondences if report.analysis is not None else ()
        self.residue_table.setRowCount(0)
        for item in values:
            row = self.residue_table.rowCount()
            self.residue_table.insertRow(row)
            cells = (
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
            for column, value in enumerate(cells):
                self.residue_table.setItem(row, column, self.w.QTableWidgetItem(value))

    def _fill_report_history(self) -> None:
        self.results_table.setRowCount(len(self._report_history))
        for row, report in enumerate(self._report_history):
            analysis = report.analysis
            if analysis is None:
                continue
            cells = (
                analysis.reference_id,
                analysis.target_id,
                analysis.alignment_decision,
                _number(analysis.sequence_identity),
                _number(analysis.sequence_coverage),
                _number(analysis.strict_rmsd_angstrom),
                _number(analysis.refined_rmsd_angstrom),
                _number(analysis.tm_score),
                str(analysis.mapped_residue_count),
                str(len(analysis.mutations)),
            )
            for column, value in enumerate(cells):
                self.results_table.setItem(row, column, self.w.QTableWidgetItem(value))
        self.results_history_status.setText(
            f"{len(self._report_history)} completed canonical report(s) · values are authoritative and descriptive."
        )

    def _fill_report_evidence_label(self, report: AnalysisReport) -> None:
        if report.evidence_cards:
            card = report.evidence_cards[0]
            self.evidence_card_label.setText(
                f"Reference {_residue_label(card.residue_id)} → Target {card.target_id or 'unavailable'} · "
                f"quality={card.quality.overall_status}; sequence, structure, interaction, site, and QC evidence are linked."
            )
        else:
            self.evidence_card_label.setText(
                "Evidence Cards unavailable — no authoritative card was reported for this comparison."
            )

    def _update_report_evidence_label(self, report: AnalysisReport) -> None:
        self._fill_report_evidence_label(report)

    def _fill_structure_result(self, result: AnalysisResult) -> None:
        self.structure_result_table.setRowCount(1)
        values = (
            result.reference_id,
            result.target_id,
            result.alignment_decision,
            _number(result.strict_rmsd_angstrom),
            _number(result.refined_rmsd_angstrom),
            _number(result.tm_score),
            str(result.mapped_residue_count),
            str(len(result.excluded_alignment_indices)),
            _backend_label(result),
        )
        for column, value in enumerate(values):
            self.structure_result_table.setItem(0, column, self.w.QTableWidgetItem(value))
        self.structure_result_summary.setText(
            f"{result.reference_id} → {result.target_id} · decision={result.alignment_decision} · "
            f"{result.mapped_residue_count} mapped residue pair(s); excluded positions remain explicit."
        )

    def _fill_results_history(self) -> None:
        self.results_table.setRowCount(len(self._analysis_history))
        for row, result in enumerate(self._analysis_history):
            values = (
                result.reference_id,
                result.target_id,
                result.alignment_decision,
                f"{result.sequence_identity:.3f}",
                f"{result.sequence_coverage:.3f}",
                _number(result.strict_rmsd_angstrom),
                _number(result.refined_rmsd_angstrom),
                _number(result.tm_score),
                str(result.mapped_residue_count),
                str(result.mutation_count),
            )
            for column, value in enumerate(values):
                self.results_table.setItem(row, column, self.w.QTableWidgetItem(value))
        count = len(self._analysis_history)
        self.results_history_status.setText(
            f"{count} completed analysis result(s) · values are authoritative and descriptive; unavailable metrics remain —."
            if count
            else "No completed analyses yet. Results from each comparison will be retained here."
        )

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
        if self.model.report is not None and self.model.report.analysis is not None:
            mutations = self.model.report.analysis.mutations
            if row < len(mutations):
                self._focus_report_alignment_index(mutations[row].alignment_index)
            return
        result = self.model.analysis
        if result is None or row >= len(result.mutations):
            return
        alignment_index = result.mutations[row].alignment_index
        self._focus_alignment_index(result, alignment_index)

    def _focus_residue(self, row: int, _: int) -> None:
        if self.model.report is not None and self.model.report.analysis is not None:
            correspondences = self.model.report.analysis.correspondences
            if row < len(correspondences):
                self._focus_report_alignment_index(correspondences[row].alignment_index)
            return
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

    def _focus_report_alignment_index(self, alignment_index: int) -> None:
        if self.model.report is None or self.model.report.analysis is None:
            return
        try:
            self._selection_controller = self._selection_controller.select_alignment_index(alignment_index)
            item = self._selection_controller.lookup_correspondence(alignment_index)
            card = None
            try:
                card = self._selection_controller.lookup_evidence_card(alignment_index)
            except (IndexError, KeyError, ValueError):
                pass
        except (IndexError, KeyError, ValueError) as exc:
            self._show_error(f"Could not select alignment {alignment_index}: {exc}")
            return
        quality = card.quality.overall_status if card is not None else "unavailable"
        self.evidence_card_label.setText(
            f"Reference {_residue_label(item.reference)} → Target {_residue_label(item.target)}\n"
            f"Status: {item.status.value} · Cα displacement: {_number(item.ca_displacement_angstrom)} Å · "
            f"Evidence Card quality: {quality}. Missing values remain unavailable."
        )
        self._set_status(f"Residue {alignment_index} selected in the canonical report")

    # -------------------------------------------------------------- persistence


__all__ = ["ReportMixin"]
