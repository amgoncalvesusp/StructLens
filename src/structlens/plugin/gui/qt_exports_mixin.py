"""Focused exports presentation/routing mixin for the Qt panel."""

# ruff: noqa: F403,F405

from __future__ import annotations

from collections.abc import Mapping

from .project_transaction import parse_visualization_state
from .qt_context import *  # noqa: F401,F403


class ExportMixin(QtMixinContext):
    """Cohesive GUI-only exports behavior composed into PanelController."""

    def _save_project(self) -> None:
        path, _ = self.w.QFileDialog.getSaveFileName(
            self.widget, "Save StructLens project", "structlens_project.json", "StructLens project (*.json)"
        )
        if not path:
            return
        if self.model.report is not None:
            try:
                binding = self._canonical_binding()
                save_canonical_project(
                    binding,
                    path,
                    visualization_state=_state_dict(self._state_from_controls()),
                )
                self._set_status(f"Verified canonical project saved · {Path(path).name}")
            except (OSError, TypeError, ValueError, ProjectSchemaError) as exc:
                self._show_error(f"Could not save project: {exc}")
            return
        result = self.model.analysis
        reference_source = self.reference_edit.text().strip() or None
        target_source = self.target_edit.text().strip() or None
        project = ProjectState(
            reference_source=(None if self.reference_object_name else reference_source),
            target_sources=(() if self.target_object_name is not None else ((target_source,) if target_source else ())),
            settings=self._settings(),
            analysis_results=(self._analysis_history or ((result,) if result is not None else ())),
            visualization_state=_state_dict(self._state_from_controls()),
            source_objects={
                key: value
                for key, value in (
                    ("reference", self.reference_object_name),
                    ("target", self.target_object_name),
                )
                if value is not None
            },
            comparison_mode=ComparisonMode(str(self.comparison_combo.currentData() or ComparisonMode.PAIRWISE.value)),
        ).with_source_hashes()
        try:
            project.save(path)
            self._set_status(f"Project saved · {Path(path).name}")
        except OSError as exc:
            self._show_error(f"Could not save project: {exc}")

    def _open_project(self) -> None:
        path, _ = self.w.QFileDialog.getOpenFileName(
            self.widget, "Open StructLens project", "", "StructLens project (*.json)"
        )
        if not path:
            return
        try:
            candidate = stage_project(path)
            if isinstance(candidate, CanonicalProjectCandidate):
                self._commit_canonical_project(candidate, Path(path))
            else:
                self._commit_legacy_project(candidate.project, Path(path))
        except (OSError, ValueError, BundleValidationError, ProjectSchemaError) as exc:
            self._show_error(f"Could not open project: {exc}")

    def _commit_canonical_project(self, candidate: CanonicalProjectCandidate, path: Path) -> None:
        project = candidate.project
        binding = candidate.binding
        reference_chain, target_chain = self._validate_project_selections(candidate)
        previous = self._capture_panel_transaction_state()
        try:
            self._invalidate_analysis_views()
            self._clear_source("reference")
            self._clear_source("target")
            self._install_loaded_source("reference", candidate.reference_source)
            self._install_loaded_source("target", candidate.target_source)
            self._set_combo_value(self.reference_model_combo, binding.request.reference_selection.model_id)
            self._set_combo_value(self.target_model_combo, binding.request.target_selection.model_id)
            self._model_changed("reference")
            self._model_changed("target")
            self._set_combo_value(self.reference_chain_combo, reference_chain)
            self._set_combo_value(self.target_chain_combo, target_chain)
            self._sync_source_model()
            self._apply_project_settings(project.settings)
            self._set_combo_value(self.comparison_combo, project.comparison_mode.value)
            self._apply_visualization_state(candidate.visualization_state)
            self._site_definitions = binding.request.site_definitions
            self._analysis_history = ()
            self._report_history = ()
            self._presentation_history = ()
            self._pending_report_request = None
            self._commit_report(binding.report, binding.presentation, binding)
            self._set_status(f"Verified canonical project opened · {path.name}")
        except Exception:
            self._restore_panel_transaction_state(previous)
            raise

    def _capture_panel_transaction_state(self) -> dict[str, Any]:
        """Capture the exact accepted/UI state for canonical-open rollback."""

        return {
            "model": self.model,
            "reference_source": self.reference_loaded_source,
            "target_source": self.target_loaded_source,
            "reference_model": self._combo_data(self.reference_model_combo),
            "target_model": self._combo_data(self.target_model_combo),
            "reference_chain": self._combo_data(self.reference_chain_combo),
            "target_chain": self._combo_data(self.target_chain_combo),
            "settings": self._settings(),
            "comparison": str(self.comparison_combo.currentData() or ComparisonMode.PAIRWISE.value),
            "visualization": self._state_from_controls(),
            "report_request": self._report_request,
            "pending_request": self._pending_report_request,
            "report_presentation": self._report_presentation,
            "report_history": self._report_history,
            "presentation_history": self._presentation_history,
            "selection_controller": self._selection_controller,
            "analysis_history": self._analysis_history,
            "site_metrics": self._site_metrics,
            "site_definitions": self._site_definitions,
            "bundle_payloads": self._v03_bundle_payloads,
            "export_records": self._v03_export_records,
            "chart_datasets": self._chart_datasets,
            "section_texts": {name: label.text() for name, label in self.report_section_labels.items()},
            "footer_status": self.footer_status.text(),
            "header_status": self.header_status.text(),
            "navigation": self.nav.currentRow(),
        }

    def _restore_panel_transaction_state(self, previous: Mapping[str, Any]) -> None:
        """Restore an accepted panel snapshot after a failed canonical commit."""

        self._restore_loaded_source("reference", previous["reference_source"])
        self._restore_loaded_source("target", previous["target_source"])
        self._set_combo_value(self.reference_model_combo, str(previous["reference_model"] or ""))
        self._set_combo_value(self.target_model_combo, str(previous["target_model"] or ""))
        self._model_changed("reference")
        self._model_changed("target")
        self._set_combo_value(self.reference_chain_combo, str(previous["reference_chain"] or ""))
        self._set_combo_value(self.target_chain_combo, str(previous["target_chain"] or ""))
        self._apply_project_settings(previous["settings"])
        self._set_combo_value(self.comparison_combo, str(previous["comparison"]))
        self._apply_visualization_state(previous["visualization"])

        model = previous["model"]
        presentation = previous["report_presentation"]
        if presentation is not None:
            history = previous["presentation_history"]
            self._presentation_history = history[:-1]
            self._render_report_presentation(presentation)
        elif model.analysis is not None:
            self._populate_result(model.analysis)
        else:
            self._invalidate_analysis_views()

        self.model = model
        self._report_request = previous["report_request"]
        self._pending_report_request = previous["pending_request"]
        self._report_presentation = presentation
        self._report_history = previous["report_history"]
        self._presentation_history = previous["presentation_history"]
        self._selection_controller = previous["selection_controller"]
        self._analysis_history = previous["analysis_history"]
        self._site_metrics = previous["site_metrics"]
        self._site_definitions = previous["site_definitions"]
        self._v03_bundle_payloads = previous["bundle_payloads"]
        self._v03_export_records = previous["export_records"]
        self._chart_datasets = previous["chart_datasets"]
        for name, text in previous["section_texts"].items():
            self.report_section_labels[name].setText(str(text))
        self._render_site_metrics(previous["site_metrics"])
        self.footer_status.setText(str(previous["footer_status"]))
        self.header_status.setText(str(previous["header_status"]))
        self.nav.setCurrentRow(int(previous["navigation"]))

    def _restore_loaded_source(self, role: str, loaded: LoadedSource | None) -> None:
        if loaded is None:
            self._clear_source(role)
        else:
            self._install_loaded_source(role, loaded)

    @staticmethod
    def _validate_project_selections(candidate: CanonicalProjectCandidate) -> tuple[str, str]:
        """Validate every source selection before the panel is mutated."""

        request = candidate.binding.request
        selected: list[str] = []
        for role, source, selection in (
            ("reference", candidate.reference_source, request.reference_selection),
            ("target", candidate.target_source, request.target_selection),
        ):
            if len(selection.author_chain_ids) != 1:
                raise ProjectSchemaError(f"canonical {role} selection must contain exactly one author chain")
            chain_id = selection.author_chain_ids[0]
            if not any(
                str(chain.model_id) == selection.model_id and chain.chain_id == chain_id
                for chain in source.structure.chains
            ):
                raise ProjectSchemaError(
                    f"canonical {role} model {selection.model_id} / chain {chain_id} is absent from its snapshot"
                )
            selected.append(chain_id)
        return selected[0], selected[1]

    def _commit_legacy_project(self, project: ProjectState, path: Path) -> None:
        self._invalidate_analysis_views()
        self._analysis_history = tuple(project.analysis_results)
        self._apply_project_settings(project.settings)
        self._set_combo_value(self.comparison_combo, project.comparison_mode.value)
        self._apply_visualization_payload(project.visualization_state)
        if project.reference_source:
            if not self._load_source("reference", Path(project.reference_source)):
                raise ValueError("legacy reference source could not be loaded")
        elif project.source_objects.get("reference") and not self._load_named_pymol_object(
            "reference", project.source_objects["reference"]
        ):
            raise ValueError("saved legacy reference PyMOL object is unavailable")
        if project.target_sources:
            if not self._load_source("target", Path(project.target_sources[0])):
                raise ValueError("legacy target source could not be loaded")
        elif project.source_objects.get("target") and not self._load_named_pymol_object(
            "target", project.source_objects["target"]
        ):
            raise ValueError("saved legacy target PyMOL object is unavailable")
        if project.analysis_results:
            self.model = self.model.with_analysis(project.analysis_results[-1])
            self._populate_result(project.analysis_results[-1])
            self.nav.setCurrentRow(_STRUCTURES_PAGE_INDEX)
        else:
            self._fill_results_history()
        self._report_request = None
        self._set_status(f"Legacy project opened · unverified · {path.name}")

    def _apply_project_settings(self, settings: AnalysisSettings) -> None:
        self._set_combo_value(self.mode_combo, settings.alignment_mode.value)
        self.identity_spin.setValue(settings.minimum_sequence_identity)
        self.coverage_spin.setValue(settings.minimum_sequence_coverage)
        self.refined_check.setChecked(settings.refined_rmsd)
        self.cutoff_spin.setValue(settings.refinement_cutoff_angstrom)
        self.usalign_edit.setText(settings.usalign_executable or "")

    def _apply_visualization_payload(self, payload: Any) -> None:
        if not payload:
            return
        state = parse_visualization_state(payload)
        self._apply_visualization_state(state)

    def _apply_visualization_state(self, state: VisualizationState) -> None:
        self._set_combo_value(self.preset_combo, state.preset)
        self._set_combo_value(self.filter_combo, state.highlight_filter.value)
        self._set_combo_value(self.color_combo, state.color_mode.value)
        self._set_combo_value(self.representation_combo, state.representation.value)
        self.radius_spin.setValue(state.local_radius_angstrom)
        self.labels_check.setChecked(state.show_labels)
        self.reference_check.setChecked(state.show_reference)
        self.target_check.setChecked(state.show_target)
        self.model = self.model.with_visualization(state)
        self._update_legend()

    # ---------------------------------------------------------------- exports

    def set_v03_export_records(self, **records: Any) -> None:
        """Stage authoritative v0.3 records for the XLSX exporter.

        Scientific services calculate these records; the GUI only stores a
        shallow copy and routes them to ``export_v03_xlsx`` when requested.
        """

        frozen = freeze_json(records)
        self._v03_export_records = dict(frozen)

    def set_v03_bundle_payloads(
        self,
        *,
        msa_summary: Mapping[str, Any] | None = None,
        conservation: Mapping[str, Any] | None = None,
        interactions: Mapping[str, Any] | None = None,
        sites: Mapping[str, Any] | None = None,
        evidence: Mapping[str, Any] | None = None,
        vectors: Mapping[str, Any] | None = None,
    ) -> None:
        """Stage authoritative v0.3 payloads for the next PyMOL export.

        The GUI deliberately accepts already calculated, JSON-ready mappings
        rather than invoking a scientific service.  ``None`` means that the
        corresponding analysis is unavailable and the bundle writer will
        omit that optional entry instead of fabricating an empty result.
        Shallow copies prevent later top-level mutations by a producer from
        changing the staged export.
        """

        payloads = {
            "msa_summary": msa_summary,
            "conservation": conservation,
            "interactions": interactions,
            "sites": sites,
            "evidence": evidence,
            "vectors": vectors,
        }
        self._v03_bundle_payloads = {
            name: None if payload is None else freeze_json(payload) for name, payload in payloads.items()
        }

    def _v03_bundle_kwargs(self) -> dict[str, Mapping[str, Any] | None]:
        """Return optional v0.3 payloads without calculating or normalizing them."""

        return {
            name: None if payload is None else thaw_json(payload) for name, payload in self._v03_bundle_payloads.items()
        }

    def _export_xlsx(self) -> None:
        # A canonical report is authoritative even when an older integration
        # staged v0.3 records on the controller.  Route it through the same
        # snapshot-aware exporter as every other canonical format.
        if self.model.report is not None:
            self._export("xlsx", export_analysis_xlsx, "XLSX")
            return
        if self._v03_export_records:
            result = self.model.analysis
            if result is None:
                self._show_error("Run a comparison before exporting results.")
                return
            path, _ = self.w.QFileDialog.getSaveFileName(
                self.widget,
                "Export v0.3 XLSX",
                "structlens_v03_result.xlsx",
                "XLSX (*.xlsx)",
            )
            if not path:
                return
            try:
                _compat_symbol("export_v03_xlsx", export_v03_xlsx)(path, **self._v03_export_records)
                self._set_status(f"v0.3 XLSX export written · {Path(path).name}")
            except (OSError, ValueError) as exc:
                self._show_error(f"Could not export v0.3 XLSX: {exc}")
            return
        self._export("xlsx", export_analysis_xlsx, "XLSX")

    def _export_csv(self) -> None:
        self._export("csv", export_analysis_csv, "CSV")

    def _export_tsv(self) -> None:
        self._export("tsv", _compat_symbol("export_analysis_tsv", export_analysis_tsv), "TSV")

    def _export_json(self) -> None:
        self._export("json", export_analysis_json, "JSON")

    def _export_pymol_bundle(self) -> None:
        if self.model.report is not None:
            try:
                self._canonical_binding()
            except (TypeError, ValueError) as exc:
                if self.model.binding is None and self._report_request is None:
                    self._show_error("Run a comparison before exporting a PyMOL bundle.")
                else:
                    self._show_error(f"Could not export canonical PyMOL bundle: {exc}")
                return
            self._show_error("Canonical snapshot-native PyMOL bundle export is unavailable until Task 15.")
            return
        result = self.model.analysis
        if result is None:
            self._show_error("Run a comparison before exporting a PyMOL bundle.")
            return
        if self.reference_structure is None or self.target_structure is None:
            self._show_error("Load reference and target coordinate files before exporting a PyMOL bundle.")
            return
        default_name = f"StructLens_{result.reference_id}_pairwise.structlens-pymol"
        path, _ = self.w.QFileDialog.getSaveFileName(
            self.widget,
            "Export StructLens-PyMOL bundle",
            default_name,
            "StructLens-PyMOL bundle (*.structlens-pymol)",
        )
        if not path:
            return
        try:
            _compat_symbol("write_pymol_bundle", write_pymol_bundle)(
                path,
                reference=self.reference_structure,
                targets={result.target_id: self.target_structure},
                analysis=result,
                provenance=dict(result.provenance),
                method_provenance=result.method_provenance,
                **self._v03_bundle_kwargs(),
            )
            self._set_status(f"Validated PyMOL bundle written · {Path(path).name}")
        except (OSError, ValueError, BundleValidationError) as exc:
            self._show_error(f"Could not export PyMOL bundle: {exc}")

    def _open_in_pymol(self) -> None:
        if self.model.report is not None:
            try:
                self._canonical_binding()
            except (TypeError, ValueError) as exc:
                self._show_error(f"Could not open canonical PyMOL bundle: {exc}")
                return
            self._show_error("Canonical snapshot-native external PyMOL launch is unavailable until Task 15.")
            return
        result = self.model.analysis
        if result is None:
            self._show_error("Run a comparison before opening PyMOL.")
            return
        if self.reference_structure is None or self.target_structure is None:
            self._show_error("Load reference and target coordinate files before opening PyMOL.")
            return
        with NamedTemporaryFile(prefix="structlens_", suffix=".structlens-pymol", delete=False) as handle:
            bundle_path = Path(handle.name)
        self._temporary_paths.append(bundle_path)
        try:
            _compat_symbol("write_pymol_bundle", write_pymol_bundle)(
                bundle_path,
                reference=self.reference_structure,
                targets={result.target_id: self.target_structure},
                analysis=result,
                provenance=dict(result.provenance),
                method_provenance=result.method_provenance,
                **self._v03_bundle_kwargs(),
            )
            launcher = PyMOLLauncher(self.pymol_edit.text().strip() or None)
            launch = launcher.launch_bundle(bundle_path)
            self._set_status(
                f"PyMOL launched with validated bundle · {bundle_path.name}"
                if launch.process_id is not None
                else f"Validated bundle prepared · {bundle_path.name}"
            )
            self._refresh_pymol_status()
        except Exception as exc:
            bundle_path.unlink(missing_ok=True)
            self._temporary_paths = [path for path in self._temporary_paths if path != bundle_path]
            self._show_error(f"Could not open in PyMOL: {exc}")

    def _export_chart_xlsx(self) -> None:
        result = self.model.analysis
        if result is None and self.model.report is None:
            self._show_error("Run a comparison before exporting chart data.")
            return
        dataset = self._selected_chart_dataset(result)
        if dataset is None:
            return
        if self.model.report is not None:
            try:
                self._canonical_binding()
            except (TypeError, ValueError) as exc:
                self._show_error(f"Could not export chart data: {exc}")
                return
        path, _ = self.w.QFileDialog.getSaveFileName(
            self.widget, "Export chart data", "structlens_chart.xlsx", "XLSX (*.xlsx)"
        )
        if not path:
            return
        try:
            _compat_symbol("export_chart_xlsx", export_chart_xlsx)(dataset, path)
            self._set_status(f"Chart data exported · {Path(path).name}")
        except (OSError, ValueError) as exc:
            self._show_error(f"Could not export chart data: {exc}")

    def _export_chart_image(self, suffix: str, dpi: int) -> None:
        result = self.model.analysis
        if result is None and self.model.report is None:
            self._show_error("Run a comparison before exporting a chart image.")
            return
        dataset = self._selected_chart_dataset(result)
        if dataset is None:
            return
        if self.model.report is not None:
            try:
                self._canonical_binding()
            except (TypeError, ValueError) as exc:
                self._show_error(f"Could not export chart image: {exc}")
                return
        path, _ = self.w.QFileDialog.getSaveFileName(
            self.widget,
            f"Export {suffix.upper()} chart ({dpi} dpi)",
            f"structlens_chart.{suffix}",
            f"{suffix.upper()} (*.{suffix})",
        )
        if not path:
            return
        try:
            _compat_symbol("export_chart_image", export_chart_image)(dataset, path, dpi=dpi)
            self._set_status(f"Chart image exported · {Path(path).name} · {dpi} dpi")
        except (OSError, RuntimeError, ValueError) as exc:
            self._show_error(f"Could not export chart image: {exc}")

    def _export(self, suffix: str, exporter: Any, label: str) -> None:
        report = self.model.report
        result = self.model.analysis
        if report is None and result is None:
            self._show_error("Run a comparison before exporting results.")
            return
        binding: CanonicalReportBinding | None = None
        if report is not None:
            # Validate the local canonical binding before opening a save dialog
            # or invoking any writer.  A stale render must have no external UI
            # side effect and must remain recoverable from the current report.
            try:
                binding = self._canonical_binding()
            except (OSError, ValueError) as exc:
                message = str(exc)
                if message == "canonical report has no verified source binding":
                    self._show_error(f"Could not export {label}: {message}")
                else:
                    self._show_error(f"Could not export {label}: canonical report binding is invalid ({message})")
                return
        path, _ = self.w.QFileDialog.getSaveFileName(
            self.widget, f"Export {label}", f"structlens_result.{suffix}", f"{label} (*.{suffix})"
        )
        if not path:
            return
        try:
            if binding is not None:
                exporter(binding.report, path, snapshots=binding.snapshots)
            else:
                exporter(result, path)
            self._set_status(f"{label} export written · {Path(path).name}")
        except (OSError, ValueError) as exc:
            self._show_error(f"Could not export {label}: {exc}")

    def _canonical_binding(self) -> CanonicalReportBinding:
        report = self.model.report
        request = self._report_request
        if report is None:
            raise ValueError("no canonical report is available")
        binding = self.model.binding
        if binding is None or request is None:
            raise ValueError("canonical report has no verified source binding")
        if binding.report is not report or request is not binding.request:
            raise ValueError("canonical report provenance is stale; a verified binding is required")
        if self._report_presentation is not binding.presentation:
            raise ValueError("canonical report presentation is stale; a verified binding is required")
        CanonicalReportBinding.create(binding.report, binding.request, binding.presentation)
        return binding

    # --------------------------------------------------------------- feedback

    def _set_status(self, message: str) -> None:
        self.header_status.setText("Busy" if self.model.busy else "Ready")
        self.footer_status.setText(message)

    def _show_error(self, message: str) -> None:
        self.model = self.model.with_error(message)
        self._set_busy(False)
        self._set_status(message)
        self.header_status.setText("Needs attention")
        # Hidden/offscreen panels still receive durable model/status feedback;
        # a modal dialog there would have no user able to dismiss it.
        if self.widget.isVisible() and hasattr(self.w.QMessageBox, "warning"):
            self.w.QMessageBox.warning(self.widget, "StructLens", message)

    def _update_mode_help(self, label: str) -> None:
        key = label.split(" ·", 1)[0]
        self.mode_help.setText(WORKFLOW_HELP.get(key, WORKFLOW_HELP["Auto"]))
        self.manual_group.setVisible(key == "Manual")

    def _update_comparison_help(self, label: str) -> None:
        key = label.split(" ·", 1)[0]
        self.comparison_help.setText(
            {
                "Pairwise": "One reference and one target. Use it for detailed residue inspection and focused WT/mutant analysis.",
            }.get(key, "Choose a comparison topology to see its outputs and computational implications.")
        )

    def _define_site_from_controls(self) -> None:
        """Store a validated site request; metrics belong to the report service."""

        positions = tuple(item.strip() for item in self.site_residues_edit.text().split(",") if item.strip())
        mode_value = self.site_mode_combo.currentData()
        if not isinstance(mode_value, str):
            self._site_unavailable("Choose a supported site definition mode.")
            return
        try:
            mode = SiteDefinitionMode(mode_value)
        except ValueError:
            self._site_unavailable("Choose a supported site definition mode.")
            return
        if mode is SiteDefinitionMode.KEY_RESIDUES and not positions:
            self._site_unavailable("Key-residue sites need at least one reference position.")
            return
        reference = self._selected_chain(
            self.reference_structure, self.reference_chain_combo, self.reference_model_combo
        )
        if reference is None:
            self._site_unavailable("Load a reference chain before defining a site.")
            return

        reference_residues = tuple(
            residue for position in positions if (residue := _find_residue(reference, position)) is not None
        )
        if mode is SiteDefinitionMode.KEY_RESIDUES and len(reference_residues) != len(positions):
            self._site_unavailable("One or more key positions are not present in the selected reference chain.")
            return

        ligand_id = self.site_ligand_edit.text().strip() or None
        if mode is SiteDefinitionMode.LIGAND_RADIUS and not ligand_id:
            self._site_unavailable("Ligand-radius sites require a ligand identifier.")
            return
        center_residue = (
            reference_residues[0] if mode is SiteDefinitionMode.RESIDUE_RADIUS and reference_residues else None
        )
        if mode is SiteDefinitionMode.RESIDUE_RADIUS and center_residue is None:
            self._site_unavailable("Residue-radius sites require one reference center position, for example A:166.")
            return
        try:
            definition = SiteDefinition(
                "panel-site",
                name="Panel site",
                mode=mode,
                reference_residues=reference_residues,
                center_residue=center_residue,
                ligand_id=ligand_id,
                radius_angstrom=(
                    None if mode is SiteDefinitionMode.KEY_RESIDUES else float(self.site_radius_spin.value())
                ),
            )
        except (TypeError, ValueError) as exc:
            self._site_unavailable(f"Invalid site definition: {exc}")
            return

        self._site_definitions = tuple(
            item for item in self._site_definitions if item.site_id != definition.site_id
        ) + (definition,)
        self.nav.setCurrentRow(SCIENTIFIC_SECTIONS.index("Sites"))
        self.site_status_label.setText("Site definition stored for the next Compare; metrics are report-generated.")

    def _site_unavailable(self, message: str) -> None:
        """Keep site absence explicit instead of presenting fabricated zeros."""

        self.set_site_metrics(())
        self.site_status_label.setText(f"Unavailable · {message}")
        self.nav.setCurrentRow(SCIENTIFIC_SECTIONS.index("Sites"))
        self._show_error(message)

    def _render_site_metrics(self, metrics: tuple[SiteMetrics, ...] | list[SiteMetrics]) -> None:
        """Render metrics supplied by the report/service without local recomputation."""

        self._site_metrics = tuple(metrics)
        self.site_metrics_table.setRowCount(len(self._site_metrics))
        for row, item in enumerate(self._site_metrics):
            values = (
                item.site_id,
                item.structure_id,
                str(item.mapped_residue_count),
                f"{item.coverage_fraction:.3f}",
                _number(item.global_frame_backbone_rmsd_angstrom),
                _number(item.site_fitted_backbone_rmsd_angstrom),
                _number(item.centroid_displacement_angstrom),
                _number(item.radius_of_gyration_angstrom),
                _number(item.atomic_envelope_volume_angstrom3),
                _number(item.sasa_angstrom2),
                _fraction_number(item.polar_residue_fraction),
                _fraction_number(item.charged_residue_fraction),
            )
            for column, value in enumerate(values):
                self.site_metrics_table.setItem(row, column, self.w.QTableWidgetItem(value))
        if self._site_metrics:
            self.site_status_label.setText(
                f"{len(self._site_metrics)} authoritative site metric record(s) loaded; units are explicit."
            )
        else:
            self.site_status_label.setText(
                "Site metrics appear after an analysis service run; this panel never estimates them locally."
            )

    def set_site_metrics(self, metrics: tuple[SiteMetrics, ...] | list[SiteMetrics]) -> None:
        """Compatibility adapter for callers that provide site metrics directly."""

        self._render_site_metrics(metrics)

    def _update_chart_export_state(self, label: str) -> None:
        """Keep export controls aligned with the chart profile they produce."""

        supported = (
            label == "Structural deviation profile"
            or label in self._chart_datasets
            or (label == "MSA conservation profile" and self._msa_chart_dataset is not None)
        )
        for button in self.chart_export_buttons:
            button.setEnabled(supported)
            if supported:
                button.setToolTip("Export the selected authoritative chart dataset.")
            else:
                button.setToolTip(
                    "This chart profile is available as structured data in the application API; "
                    "GUI image/XLSX export is currently implemented for Structural deviation profile."
                )

    def _set_combo_value(self, combo: Any, value: str) -> None:
        for index in range(combo.count()):
            if str(combo.itemData(index)) == value or combo.itemText(index) == value:
                combo.setCurrentIndex(index)
                return

    def _refresh_pymol_status(self) -> None:
        try:
            executable = PyMOLLauncher(self.pymol_edit.text().strip() or None).locate()
        except Exception:
            self.pymol_status.setText(
                "Ready to export · PyMOL executable not configured or not found. Export remains available even without PyMOL."
            )
            return
        self.pymol_status.setText(
            f"Ready to open in PyMOL · executable resolved to {executable.name}. Scientific calculations stay in StructLens."
        )

    def close(self) -> None:
        if self._future is not None:
            self._cancel_event.set()
            self._future.cancel()
            self._future = None
        if self._poll_timer is not None:
            self._poll_timer.stop()
            self._poll_timer = None
        self._pymol.reset()
        for path in self._temporary_paths:
            path.unlink(missing_ok=True)

    l__ = ["PanelController", "build_panel"]


__all__ = ["ExportMixin"]
