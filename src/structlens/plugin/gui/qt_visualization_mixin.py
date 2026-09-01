"""Focused visualization presentation/routing mixin for the Qt panel."""

# ruff: noqa: F403,F405

from __future__ import annotations

from .qt_context import *  # noqa: F401,F403


class VisualizationMixin(QtMixinContext):
    """Cohesive GUI-only visualization behavior composed into PanelController."""

    def _state_from_controls(self) -> VisualizationState:
        return VisualizationState(
            highlight_filter=HighlightFilter(str(self.filter_combo.currentData())),
            color_mode=ColorMode(str(self.color_combo.currentData())),
            representation=Representation(str(self.representation_combo.currentData())),
            show_labels=self.labels_check.isChecked(),
            show_reference=self.reference_check.isChecked(),
            show_target=self.target_check.isChecked(),
            local_radius_angstrom=float(self.radius_spin.value()),
            preset=self.preset_combo.currentText(),
        )

    def _visualization_changed(self, *_: object) -> None:
        state = self._state_from_controls()
        self.model = self.model.with_visualization(state)
        self._update_legend()

    def _apply_preset(self, name: str) -> None:
        if not name:
            return
        state = self._renderer.apply_preset(name)
        self._set_combo_value(self.filter_combo, state.highlight_filter.value)
        self._set_combo_value(self.color_combo, state.color_mode.value)
        self._set_combo_value(self.representation_combo, state.representation.value)
        self.labels_check.setChecked(state.show_labels)
        self.reference_check.setChecked(state.show_reference)
        self.target_check.setChecked(state.show_target)
        self.model = self.model.with_visualization(state)
        self._update_legend()

    def _apply_visualization(self) -> None:
        binding = self.model.binding
        result = binding.report if binding is not None else self.model.analysis
        if result is None:
            self._show_error("Run a comparison before applying a PyMOL view.")
            return
        if self.command is None:
            self._show_error("Apply a visualization from inside a PyMOL session.")
            return
        state = self._state_from_controls()
        if binding is None and (
            self.target_object_name is None or (state.show_reference and self.reference_object_name is None)
        ):
            self._show_error(
                "The selected sources are not loaded as PyMOL objects. Use Browse in PyMOL or choose Use object… first."
            )
            return
        verified = None
        if binding is not None:
            try:
                verified = self._canonical_binding()
            except (TypeError, ValueError) as exc:
                self._show_error(f"Could not apply canonical PyMOL view: {exc}")
                return
        self.model = self.model.with_visualization(state)
        selected = len(self._visualization_service.select(result, state))
        if verified is not None:
            apply_report = getattr(self._pymol, "apply_report", None)
            if apply_report is not None:
                apply_report(
                    verified.report,
                    snapshots=verified.snapshots,
                    state=state,
                    reference_object=self.reference_object_name,
                    target_object=self.target_object_name,
                )
            else:
                # Compatibility for pre-report adapter doubles/integrations.
                # The production PyMOLAdapter always takes the report path.
                self._pymol.apply(
                    verified.report,
                    state=state,
                    reference_object=self.reference_object_name,
                    target_object=self.target_object_name,
                )
        else:
            self._pymol.apply(
                result,
                state=state,
                reference_object=self.reference_object_name,
                target_object=self.target_object_name,
            )
        self._set_status(f"PyMOL view applied · {selected} correspondence rows selected")

    def _reset_visualization(self) -> None:
        self._pymol.reset()
        self._set_status("StructLens-owned PyMOL selections removed")

    def _update_legend(self) -> None:
        binding = self.model.binding
        result = binding.report if binding is not None else self.model.analysis
        if result is None:
            self.legend_label.setText("Legend appears after a comparison.")
            self.visualization_count.setText("0 rows selected")
            return
        state = self._state_from_controls()
        selected = self._visualization_service.select(result, state)
        correspondences = (
            (() if result.analysis is None else result.analysis.correspondences)
            if isinstance(result, AnalysisReport)
            else result.correspondences
        )
        self.visualization_count.setText(f"{len(selected)} of {len(correspondences)} rows selected")
        if state.color_mode is ColorMode.CA_DISPLACEMENT:
            values = [
                item.ca_displacement_angstrom for item in correspondences if item.ca_displacement_angstrom is not None
            ]
            if values:
                legend = displacement_legend(min(values), max(values))
                self.legend_label.setText(
                    f"Legend · {legend.title} · {legend.unit} · {legend.minimum:.2f}–{legend.maximum:.2f}"
                )
                return
        if state.color_mode is ColorMode.BACKBONE_RMSD:
            values = [
                item.backbone_rmsd_angstrom for item in correspondences if item.backbone_rmsd_angstrom is not None
            ]
            if values:
                legend = backbone_rmsd_legend(min(values), max(values))
                self.legend_label.setText(
                    f"Legend · {legend.title} · {legend.unit} · {legend.minimum:.2f}–{legend.maximum:.2f}"
                )
                return
        self.legend_label.setText(
            "Legend · status colors are paired with text labels; no scientific meaning is inferred from color alone."
        )

    # -------------------------------------------------------------- result UI

    def _render_msa_chart(self, alignment: MultipleSequenceAlignment) -> None:
        series = ChartSeries(
            "Alignment conservation",
            tuple((float(column.index + 1), column.conservation_score) for column in alignment.columns),
            tuple(column.reference_label for column in alignment.columns),
        )
        self._msa_chart_dataset = ChartDataset(
            "msa_conservation_profile",
            "MSA conservation profile",
            "Alignment column",
            "Alignment conservation",
            "fraction",
            (series,),
            "Authoritative alignment conservation; gaps and ambiguous residues are not amino-acid observations.",
        )
        self._render_dataset(
            self._msa_chart_dataset,
            self.sequence_chart_layout,
            self.sequence_chart_status,
            canvas_attribute="_sequence_chart_canvas",
            unavailable_message="MSA chart data is available, but no renderer is installed.",
        )

    def _render_selected_chart(self) -> None:
        label = self.chart_combo.currentText()
        result = self.model.analysis
        dataset: ChartDataset | MatrixDataset | None = None
        if label == "Structural deviation profile" and result is not None:
            dataset = structural_deviation_profile(result)
        elif label == "MSA conservation profile":
            dataset = self._chart_datasets.get(label) or self._msa_chart_dataset
        else:
            dataset = self._chart_datasets.get(label)
        if dataset is None:
            self.chart_preview_status.setText("Chart unavailable. Run the corresponding scientific service first.")
            self._clear_chart_layout(self.chart_preview_layout)
            return
        self._render_dataset(
            dataset,
            self.chart_preview_layout,
            self.chart_preview_status,
            canvas_attribute="_chart_canvas",
            unavailable_message=(
                "Authoritative chart data is ready; install the optional charts dependency (matplotlib) "
                "to render the interactive preview."
            ),
        )

    def _render_dataset(
        self,
        dataset: ChartDataset | MatrixDataset,
        layout: Any,
        status: Any,
        *,
        canvas_attribute: str,
        unavailable_message: str,
    ) -> None:
        self._clear_chart_layout(layout)
        chart_classes = load_chart_classes()
        if chart_classes is None:
            setattr(self, canvas_attribute, None)
            status.setText(f"{unavailable_message} Values remain available for XLSX/CSV export.")
            return
        FigureCanvasQTAgg, Figure = chart_classes
        figure = Figure(figsize=(7.0, 2.8), dpi=100, tight_layout=True)
        axes = figure.add_subplot(111)
        if isinstance(dataset, ChartDataset):
            y_label = _unit_label(dataset.y_label, dataset.unit)
            for series in dataset.series:
                points = [(x, y) for x, y in series.points if y is not None]
                if points:
                    axes.plot(
                        [point[0] for point in points],
                        [point[1] for point in points],
                        marker="o",
                        linewidth=1.4,
                        label=series.name,
                    )
            axes.set_xlabel(dataset.x_label)
            axes.set_ylabel(y_label)
            axes.set_title(dataset.title)
            if len(dataset.series) > 1:
                axes.legend()
        else:
            rows = list(dict.fromkeys(cell.row for cell in dataset.cells))
            columns = list(dict.fromkeys(cell.column for cell in dataset.cells))
            values = {(cell.row, cell.column): cell.value for cell in dataset.cells}
            image: list[list[float]] = []
            for row in rows:
                image_row: list[float] = []
                for column in columns:
                    value = values.get((row, column))
                    image_row.append(float("nan") if value is None else float(value))
                image.append(image_row)
            if image and columns:
                axes.imshow(image, **_matrix_image_kwargs(dataset, values.values()))
                axes.set_xticks(range(len(columns)), columns, rotation=45, ha="right")
                axes.set_yticks(range(len(rows)), rows)
            axes.set_xlabel(dataset.column_label)
            axes.set_ylabel(dataset.row_label)
            axes.set_title(dataset.title)
        canvas = FigureCanvasQTAgg(figure)
        layout.addWidget(canvas)
        setattr(self, canvas_attribute, canvas)
        unit = dataset.unit if isinstance(dataset, ChartDataset) else None
        unit_text = f" · unit={unit}" if unit else ""
        status.setText(f"{dataset.title}{unit_text} · {dataset.interpretation}")

    def _clear_chart_layout(self, layout: Any) -> None:
        while layout.count():
            item = layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.setParent(None)
                widget.deleteLater()

    def _render_sequence_result(self, result: AnalysisResult) -> None:
        self._render_dataset(
            mutation_conservation_matrix(result),
            self.sequence_chart_layout,
            self.sequence_chart_status,
            canvas_attribute="_sequence_chart_canvas",
            unavailable_message="Mutation/conservation chart unavailable for this result.",
        )

    def _update_chart_explanation(self, label: str) -> None:
        explanations = {
            "Structural deviation profile": "Reference position on X; selectable Cα displacement, backbone RMSD, side-chain RMSD, or local RMSD on Y (Å).",
            "Mutation / conservation matrix": "Rows are structures and columns are reference-aligned positions; cell text preserves mutation identity.",
            "Pairwise similarity heatmap": "Mirrors one stored value per pair for sequence identity, TM-score, RMSD, or key-site RMSD.",
            "Sequence–structure relationship": "Sequence identity (%) versus TM-score by default; points link back to a selected pair.",
            "Structural conservation profile": "Cα positional variability (Å) with a separate position-coverage track.",
            "Key-residue comparison": "Explicit key reference residues compared across targets with units and missing mappings visible.",
            "MSA conservation profile": "Alignment-column conservation with gap and ambiguous-residue tracks; gaps are not amino acids.",
            "Sequence logo": "Letter heights are pᵢ × information content, with canonical amino acids only.",
            "Interaction difference matrix": "Reference-normalized conserved, gained, lost, and target-only-unmapped interactions.",
            "Site comparison": "Global-frame and site-fitted RMSD remain separate; envelope volume is not a cavity volume.",
            "Distance-difference heatmap": "Target internal distance minus reference internal distance (Å), with missing pairs masked.",
            "Evidence Card figure": "Publication-ready descriptive evidence; no impact or functional score is inferred.",
        }
        self.chart_explanation.setText(explanations.get(label, "Charts consume the authoritative analysis state."))
        self._update_chart_export_state(label)
        if hasattr(self, "chart_preview_status"):
            self._render_selected_chart()

    def set_chart_datasets(self, datasets: Mapping[str, ChartDataset | MatrixDataset]) -> None:
        """Stage authoritative chart datasets for all v0.3 publication exports."""

        self._chart_datasets = dict(datasets)
        self._update_chart_export_state(self.chart_combo.currentText())
        self._render_selected_chart()
        if self._chart_datasets:
            self.nav.setCurrentRow(_CHARTS_PAGE_INDEX)

    def _selected_chart_dataset(self, result: AnalysisResult | None) -> ChartDataset | MatrixDataset | None:
        label = self.chart_combo.currentText()
        if label == "Structural deviation profile" and result is not None:
            return structural_deviation_profile(result)
        dataset = self._chart_datasets.get(label)
        if dataset is None and label == "MSA conservation profile":
            dataset = self._msa_chart_dataset
        if dataset is None:
            self._show_error(f"The {label} dataset is unavailable. Run its scientific service before exporting.")
        return dataset


__all__ = ["VisualizationMixin"]


def _unit_label(label: str, unit: str | None) -> str:
    if unit is None or unit in label:
        return label
    return f"{label} ({unit})"
