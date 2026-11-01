"""Qt widgets and rendering for the immutable report presentation."""

# ruff: noqa: F403,F405

from __future__ import annotations

import math

from structlens.application.chart_data import DistanceChartCells

from .qt_context import *  # noqa: F401,F403


class PresentationMixin(QtMixinContext):
    """Render canonical report widgets from ``ReportPresentation`` only."""

    def _build_presentation_sections(self, content: Any) -> None:
        group = self.w.QGroupBox("Scientific details", self.widget)
        layout = self.w.QVBoxLayout(group)
        layout.setContentsMargins(18, 16, 18, 16)
        self.report_section_labels: dict[str, Any] = {}
        for name, title in _SECTION_TITLES:
            section = self.w.QGroupBox(title, group)
            section_layout = self.w.QVBoxLayout(section)
            value = _label(self.w, "Run Compare structures to view available results.", "inlineNote")
            value.setWordWrap(True)
            section_layout.addWidget(value)
            layout.addWidget(section)
            self.report_section_labels[name] = value
        content.addWidget(group)

    def _render_report_presentation(self, presentation: ReportPresentation) -> None:
        if not isinstance(presentation, ReportPresentation):
            raise TypeError("presentation must be a ReportPresentation")
        for name, title in _SECTION_TITLES:
            section = getattr(presentation.sections, name)
            self.report_section_labels[name].setText(_section_text(title, section))
        self._render_presentation_summary(presentation)
        self._render_presentation_msa(presentation)
        self._render_presentation_mutations(presentation)
        self._render_presentation_correspondences(presentation)
        self._render_presentation_sites(presentation)
        self._render_quality()
        self._render_pocket_presentation()
        self.evidence_card_label.setText("Select a residue to inspect its evidence.")
        self._presentation_history = (*self._presentation_history, presentation)
        self._render_presentation_history()
        self._install_presentation_charts(presentation)
        self._refresh_ligand_contacts()

    def _render_presentation_summary(self, presentation: ReportPresentation) -> None:
        rows = presentation.sections.summary.rows
        row = rows[0] if rows else {}
        reference = str(row.get("reference", "Unavailable"))
        target = str(row.get("target", "Unavailable"))
        decision = str(row.get("decision", "Unavailable"))
        self._set_report_decision(decision, reference, target)
        for key in (
            "sequence_identity",
            "sequence_similarity",
            "sequence_coverage",
            "strict_rmsd_angstrom",
            "refined_rmsd_angstrom",
            "tm_score",
            "mapped_residue_count",
        ):
            self.result_labels[key].setText(_display_metric(key, row.get(key)))
        analysis = self.model.report.analysis if self.model.report is not None else None
        self._set_mutation_counts(analysis.mutations if analysis is not None else None)
        self.structure_result_table.setRowCount(1)
        structure_values = (
            reference,
            target,
            decision,
            _display_metric("strict_rmsd_angstrom", row.get("strict_rmsd_angstrom")),
            _display_metric("refined_rmsd_angstrom", row.get("refined_rmsd_angstrom")),
            _display_metric("tm_score", row.get("tm_score")),
            str(row.get("mapped_residue_count", "Unavailable")),
            "See correspondence status",
            "comparison report",
        )
        for column, value in enumerate(structure_values):
            self.structure_result_table.setItem(0, column, self.w.QTableWidgetItem(value))
        self.structure_result_summary.setText(presentation.sections.summary.summary)

    def _render_presentation_msa(self, presentation: ReportPresentation) -> None:
        section = presentation.sections.msa
        self.msa_table.setRowCount(len(section.rows))
        for row_index, row in enumerate(section.rows):
            for column, key in enumerate(("structure", "aligned_sequence", "source")):
                self.msa_table.setItem(row_index, column, self.w.QTableWidgetItem(str(row.get(key, ""))))
        self.msa_summary_label.setText(section.summary)
        self.sequence_chart_status.setText(section.summary)

    def _render_presentation_mutations(self, presentation: ReportPresentation) -> None:
        rows = presentation.sections.mutations.rows
        self.mutation_table.setRowCount(len(rows))
        keys = ("index", "kind", "reference", "target", "notation", "blosum62", "grantham", "class")
        for row_index, row in enumerate(rows):
            for column, key in enumerate(keys):
                self.mutation_table.setItem(row_index, column, self.w.QTableWidgetItem(str(row.get(key, ""))))
        self.mutation_summary.setText(presentation.sections.mutations.summary)

    def _render_presentation_correspondences(self, presentation: ReportPresentation) -> None:
        analysis = self.model.report.analysis if self.model.report is not None else None
        self._set_alignment_rows(analysis.correspondences if analysis else ())
        rows = presentation.sections.correspondences.rows
        self.residue_table.setRowCount(len(rows))
        keys = (
            "index",
            "reference",
            "target",
            "status",
            "ca_displacement_angstrom",
            "backbone_rmsd_angstrom",
            "sidechain_rmsd_angstrom",
            "all_heavy_atom_rmsd_angstrom",
            "outlier",
            "key",
        )
        for row_index, row in enumerate(rows):
            for column, key in enumerate(keys):
                value = row.get(key, "")
                display = _display_metric(key, value) if key.endswith("_angstrom") else str(value)
                self.residue_table.setItem(row_index, column, self.w.QTableWidgetItem(display))
        self.residue_summary.setText(presentation.sections.correspondences.summary)

    def _render_presentation_sites(self, presentation: ReportPresentation) -> None:
        rows = presentation.sections.sites.rows
        self.site_metrics_table.setRowCount(len(rows))
        keys = (
            "site_id",
            "structure_id",
            "mapped_residue_count",
            "coverage_fraction",
            "global_frame_backbone_rmsd_angstrom",
            "site_fitted_backbone_rmsd_angstrom",
            "centroid_displacement_angstrom",
            "radius_of_gyration_angstrom",
            "atomic_envelope_volume_angstrom3",
            "sasa_angstrom2",
            "polar_residue_fraction",
            "charged_residue_fraction",
        )
        for row_index, row in enumerate(rows):
            for column, key in enumerate(keys):
                value = str(row.get(key, "Unavailable — value not reported"))
                self.site_metrics_table.setItem(row_index, column, self.w.QTableWidgetItem(value))
        self.site_status_label.setText(presentation.sections.sites.summary)

    def _render_presentation_history(self) -> None:
        self.results_table.setRowCount(len(self._presentation_history))
        for row_index, presentation in enumerate(self._presentation_history):
            rows = presentation.sections.summary.rows
            row = rows[0] if rows else {}
            values = (
                row.get("reference", "Unavailable"),
                row.get("target", "Unavailable"),
                row.get("decision", "Unavailable"),
                _display_metric("sequence_identity", row.get("sequence_identity")),
                _display_metric("sequence_coverage", row.get("sequence_coverage")),
                _display_metric("strict_rmsd_angstrom", row.get("strict_rmsd_angstrom")),
                _display_metric("refined_rmsd_angstrom", row.get("refined_rmsd_angstrom")),
                _display_metric("tm_score", row.get("tm_score")),
                row.get("mapped_residue_count", "Unavailable"),
                row.get("mutation_count", "Unavailable"),
            )
            for column, value in enumerate(values):
                self.results_table.setItem(row_index, column, self.w.QTableWidgetItem(str(value)))
        self.results_history_status.setText(
            f"{len(self._presentation_history)} completed comparison(s). Full precision is retained in exports."
        )

    def _install_presentation_charts(self, presentation: ReportPresentation) -> None:
        datasets: dict[str, ChartDataset | MatrixDataset] = {}
        chart = presentation.chart_datasets["structural_deviation_profile"]
        usable = chart["points"]
        use_alignment_axis = any(item["position"] is None and item["metadata"]["reference_residue"] is not None for item in usable)
        points = tuple(
            (float(item["metadata"]["alignment_index"] + 1) if use_alignment_axis else
             float(item["position"]) if item["position"] is not None else float("nan"),
             None if item["value"] is None else float(item["value"]))
            for item in usable
        )
        labels = tuple(str(item["label"]) for item in usable)
        metadata = tuple(item["metadata"] for item in usable)
        unit = str(chart["unit"])
        datasets["Structural deviation profile"] = ChartDataset(
            "structural_deviation_profile",
            "Structural deviation profile",
            "Alignment position" if use_alignment_axis else "Reference position",
            f"Cα displacement ({unit})",
            unit,
            (ChartSeries("Cα displacement", points, labels, metadata),),
            "Authoritative report values; unavailable positions remain unavailable.",
        )
        msa = presentation.chart_datasets["msa_conservation_profile"]
        report = self.model.report
        correspondences = report.analysis.correspondences if report and report.analysis else ()
        by_residue = {rid: c.alignment_index for c in correspondences for rid in (c.reference, c.target) if rid}
        msa_columns = report.msa.columns if report and report.msa else ()
        msa_metadata = []
        for column in msa_columns:
            # MSA column numbers can differ from the structural alignment.
            indices = {by_residue[cell.residue.residue_id] for cell in column.cells
                       if cell.residue and cell.residue.residue_id in by_residue}
            index = by_residue.get(column.reference_residue) if column.reference_residue else (
                next(iter(indices)) if len(indices) == 1 else None
            )
            msa_metadata.append({"alignment_index": index})
        msa_points = tuple(
            (float(item["index"] + 1), None if item["value"] is None else float(item["value"]))
            for item in msa["columns"]
        )
        datasets["MSA conservation profile"] = ChartDataset(
            "msa_conservation_profile",
            "MSA conservation profile",
            "Alignment column",
            f"Conservation ({msa['unit']})",
            str(msa["unit"]),
            (ChartSeries("Conservation", msa_points, tuple(c.reference_label for c in msa_columns), tuple(msa_metadata)),),
            presentation.sections.msa.summary,
        )
        if report is not None and report.distance_map is not None:
            datasets["Distance-difference heatmap"] = MatrixDataset(
                "distance_difference_heatmap", "Distance-difference heatmap", "Reference residue", "Reference residue",
                DistanceChartCells(report.distance_map),
                "Target − reference (Å). Blue indicates shorter distances, red longer distances, and grey missing values. Select a cell to inspect both residues.",
            )
        if report is not None and report.analysis is not None:
            legacy = _legacy_analysis(report)
            if legacy is not None:
                datasets["Mutation / conservation matrix"] = mutation_conservation_matrix(legacy)
        self._chart_datasets = datasets
        self._render_selected_chart()
        self._update_chart_export_state(self.chart_combo.currentText())

    def _render_quality(self) -> None:
        report = self.model.report
        if report is None or not hasattr(self, "quality_table"):
            return
        try:
            filter_value = QualityFilter(str(self.quality_filter.currentData() or QualityFilter.ALL.value))
        except ValueError:
            filter_value = QualityFilter.ALL
        reference = QualityPresenter.present(report.input_quality.reference, filter_value=filter_value)
        target = QualityPresenter.present(report.input_quality.target, filter_value=filter_value)
        rows = tuple(("reference", row) for row in reference.rows) + tuple(
            ("target", row) for row in target.rows
        )
        self.quality_table.setRowCount(len(rows))
        for index, (role, row) in enumerate(rows):
            location = row.residue_id or row.atom_id or "—"
            values = (role, row.severity.value, row.code, location, row.message, row.remediation or "—")
            for column, value in enumerate(values):
                self.quality_table.setItem(index, column, self.w.QTableWidgetItem(str(value)))
        summary = reference.summary
        target_summary = target.summary
        self.quality_summary_label.setText(
            f"Reference: {summary.error_count} error(s), {summary.warning_count} warning(s), "
            f"{summary.info_count} info · Target: {target_summary.error_count} error(s), "
            f"{target_summary.warning_count} warning(s), {target_summary.info_count} info"
        )

    def _render_pocket_presentation(self) -> None:
        report = self.model.report
        if report is None or not hasattr(self, "pocket_table"):
            return
        snapshot = report.pockets
        if snapshot is None:
            self.pocket_table.setRowCount(0)
            self.pocket_status_label.setText("Pocket evidence is unavailable in this canonical report.")
            self.pocket_capability_label.setText("Rerun Compare after selecting both structures.")
            self.pocket_detect_button.setEnabled(False)
            self.pocket_measure_button.setEnabled(False)
            return
        presentation = PocketPresenter.present(snapshot)
        self._pocket_presentation = presentation
        self.pocket_table.setRowCount(len(presentation.candidate_rows))
        for index, row in enumerate(presentation.candidate_rows):
            detail = presentation.detail_for(row.candidate_id)
            values = (
                row.candidate_id[:12],
                row.role,
                row.state.value if row.state is not None else "unmatched",
                str(row.sphere_count),
                str(row.lining_residue_count),
                _number(detail.fine_volume_angstrom3 if detail else None),
                _number(detail.absolute_sensitivity_angstrom3 if detail else None),
            )
            for column, value in enumerate(values):
                item = self.w.QTableWidgetItem(str(value))
                item.setData(self.c.Qt.UserRole, row.candidate_id)
                self.pocket_table.setItem(index, column, item)
        self.pocket_status_label.setText(
            f"{len(presentation.candidate_rows)} candidate(s) · {len(presentation.match_rows)} explicit match state(s) · "
            f"report {report.report_id[:12]}"
        )
        capability = "; ".join(presentation.disabled_reasons)
        if not capability:
            capability = "Detection and volume measurement are read-only here; rerun Compare to compute a new canonical report."
        self.pocket_capability_label.setText(
            f"Status: {presentation.availability.value}. {capability}"
        )
        # These controls intentionally do not start a second GUI-local
        # calculation.  The current report remains immutable; a new Compare
        # is required to recompute detection or volume evidence.
        self.pocket_detect_button.setEnabled(False)
        self.pocket_measure_button.setEnabled(False)

    def _select_pocket_candidate(self, row: int, _column: int) -> None:
        presentation = getattr(self, "_pocket_presentation", None)
        if presentation is None:
            return
        item = self.pocket_table.item(row, 0)
        if item is None:
            return
        candidate_id = str(item.data(self.c.Qt.UserRole) or "")
        detail = presentation.detail_for(candidate_id)
        if detail is None:
            return
        self.pocket_detail_label.setText(
            f"{detail.role.title()} candidate {detail.candidate_id[:12]} · "
            f"coarse { _number(detail.coarse_volume_angstrom3)} Å³ · "
            f"fine { _number(detail.fine_volume_angstrom3)} Å³ · "
            f"sensitivity { _number(detail.absolute_sensitivity_angstrom3)} Å³\n"
            f"{detail.method_explanation}"
        )

    def _request_pocket_detection(self) -> None:
        if getattr(self, "_pocket_presentation", None) is None:
            self._show_error("Pocket detection is unavailable until a canonical comparison report exists.")
            return
        self._set_status("Pocket detection is part of Compare; rerun Compare to change pocket settings.")

    def _request_pocket_measurement(self) -> None:
        presentation = getattr(self, "_pocket_presentation", None)
        if presentation is None or not presentation.can_measure:
            reasons = presentation
            message = "Volume measurement is unavailable for the current pocket evidence."
            if reasons is not None and reasons.disabled_reasons:
                message = reasons.disabled_reasons[0]
            self._show_error(message)
            return
        self._set_status("Volume measurements are included in the canonical pocket report.")


def _display_metric(key: str, value: Any) -> str:
    """Format the view only; reports and exports keep their full precision."""
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return str(value) if value is not None else "Unavailable — value not reported"
    if not math.isfinite(numeric):
        return "Unavailable — value not reported"
    if key.startswith("sequence_"):
        return f"{numeric:.1%}"
    if key in {"mapped_residue_count", "mutation_count"}:
        return str(int(numeric))
    return f"{numeric:.3f}"


def _section_text(title: str, section: Any) -> str:
    rows = "\n".join(" · ".join(f"{key}={value}" for key, value in row.items()) for row in section.rows)
    base = f"{title} · status={section.status.value} · {section.summary}"
    return f"{base}\n{rows}" if rows else base


_SECTION_TITLES = (
    ("summary", "Report summary"),
    ("quality", "Input QC"),
    ("msa", "MSA and mutations"),
    ("correspondences", "Correspondences"),
    ("evidence_cards", "Evidence Card"),
    ("interactions", "Interaction evidence"),
    ("sites", "Site evidence"),
    ("distance_map", "Distance map"),
    ("vectors", "Vector evidence"),
    ("diagnostics", "Diagnostics"),
    ("export", "Export status"),
    ("pymol", "PyMOL data"),
)


__all__ = ["PresentationMixin"]
