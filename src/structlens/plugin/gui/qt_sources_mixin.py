"""Focused sources presentation/routing mixin for the Qt panel."""

# ruff: noqa: F403,F405

from __future__ import annotations

from .qt_context import *  # noqa: F401,F403


class SourceMixin(QtMixinContext):
    """Cohesive GUI-only sources behavior composed into PanelController."""

    def _browse_source(self, role: str) -> None:
        path, _ = self.w.QFileDialog.getOpenFileName(
            self.widget,
            f"Choose {role} structure",
            "",
            "Structure files (*.pdb *.ent *.cif *.mmcif *.gz);;All files (*)",
        )
        if path:
            edit = self.reference_edit if role == "reference" else self.target_edit
            edit.setText(path)
            self._load_source(role, Path(path))

    def _load_sources_from_edits(self) -> bool:
        loaded = True
        if self.reference_edit.text().strip():
            loaded = self._load_source("reference", Path(self.reference_edit.text().strip())) and loaded
        else:
            self._clear_source("reference")
            loaded = False
        if self.target_edit.text().strip():
            loaded = self._load_source("target", Path(self.target_edit.text().strip())) and loaded
        else:
            self._clear_source("target")
            loaded = False
        return loaded

    def _load_source(
        self,
        role: str,
        path: Path,
        *,
        pymol_object_name: str | None = None,
    ) -> bool:
        try:
            loaded = LoadedSource.capture(path, pymol_object_name=pymol_object_name)
            structure = loaded.structure
            if not structure.chains:
                raise ValueError("source contains no coordinate records or residue chains")
            object_name = pymol_object_name
            loaded = LoadedSource(loaded.path, loaded.snapshot, structure, object_name)
            self._invalidate_analysis_views()
            self._clear_source(role)
            self._install_loaded_source(role, loaded)
            self._set_status(f"{role.title()} loaded · choose chains or run comparison")
            return True
        except Exception as exc:  # parser errors are user-facing recovery states
            self._invalidate_analysis_views()
            self._clear_source(role)
            self._show_error(f"Could not load {role}: {exc}")
            return False

    def _install_loaded_source(self, role: str, loaded: LoadedSource) -> None:
        """Commit a fully parsed source to widgets after validation succeeds."""

        if loaded.pymol_object_name is None:
            loaded = LoadedSource(
                loaded.path,
                loaded.snapshot,
                loaded.structure,
                self._load_snapshot_into_pymol(role, loaded),
            )
        structure = loaded.structure
        if role == "reference":
            self.reference_loaded_source = loaded
            self.reference_structure = structure
            self.reference_object_name = loaded.pymol_object_name
            self.reference_edit.setText(str(loaded.path))
            self._populate_models(self.reference_model_combo, structure)
            self._populate_chains(
                self.reference_chain_combo,
                structure,
                model_id=self._combo_data(self.reference_model_combo),
            )
            self.reference_meta.setText(_structure_meta(structure))
        else:
            self.target_loaded_source = loaded
            self.target_structure = structure
            self.target_object_name = loaded.pymol_object_name
            self.target_edit.setText(str(loaded.path))
            self._populate_models(self.target_model_combo, structure)
            self._populate_chains(
                self.target_chain_combo,
                structure,
                model_id=self._combo_data(self.target_model_combo),
            )
            self.target_meta.setText(_structure_meta(structure))
        self._sync_source_model()

    def _load_snapshot_into_pymol(self, role: str, loaded: LoadedSource) -> str | None:
        """Materialize the captured logical bytes, never the mutable source path."""

        if self.command is None:
            return None
        load = getattr(self.command, "load", None)
        if load is None:
            return None
        suffix = ".pdb" if loaded.snapshot.logical_format == "pdb" else ".cif"
        with NamedTemporaryFile(prefix="structlens_snapshot_", suffix=suffix, delete=False) as handle:
            handle.write(loaded.snapshot.decompressed_bytes)
            temporary = Path(handle.name)
        self._temporary_paths.append(temporary)
        object_name = selection_name(
            "panel",
            f"{Path(loaded.snapshot.display_name).stem}_{loaded.snapshot.content_id[:12]}",
            role,
        )
        try:
            get_names = getattr(self.command, "get_names", None)
            names = {str(name) for name in get_names("objects")} if get_names is not None else set()
            object_name = _unique_pymol_name(object_name, names)
            load(str(temporary), object_name)
            return str(object_name)
        except Exception:
            return None
        finally:
            # The logical snapshot is already retained in LoadedSource; the
            # materialized coordinate file must not survive a load failure or
            # wait for panel destruction.
            temporary.unlink(missing_ok=True)

    def _clear_source(self, role: str) -> None:
        if role == "reference":
            self.reference_loaded_source = None
            self.reference_structure = None
            self.reference_object_name = None
            self.reference_model_combo.clear()
            self.reference_chain_combo.clear()
            self.reference_meta.setText("Not loaded")
        else:
            self.target_loaded_source = None
            self.target_structure = None
            self.target_object_name = None
            self.target_model_combo.clear()
            self.target_chain_combo.clear()
            self.target_meta.setText("Not loaded")
        self._sync_source_model()

    def _use_pymol_object(self, role: str) -> None:
        if self.command is None:
            self._show_error("PyMOL object sources are available only inside a PyMOL session.")
            return
        get_names = getattr(self.command, "get_names", None)
        save = getattr(self.command, "save", None)
        if get_names is None or save is None:
            self._show_error("The PyMOL command proxy cannot list and save objects.")
            return
        names = [str(name) for name in get_names("objects")]
        if not names:
            self._show_error("No PyMOL objects are available. Load a structure first.")
            return
        selected, accepted = self.w.QInputDialog.getItem(
            self.widget,
            f"Choose {role} object",
            "PyMOL object",
            names,
            0,
            False,
        )
        if not accepted or not selected:
            return
        with NamedTemporaryFile(prefix="structlens_", suffix=".pdb", delete=False) as handle:
            temporary = Path(handle.name)
        try:
            save(str(temporary), selected)
            self._temporary_paths.append(temporary)
            self._load_source(role, temporary)
        except Exception as exc:
            self._show_error(f"Could not read PyMOL object {selected}: {exc}")
        finally:
            temporary.unlink(missing_ok=True)

    def _load_named_pymol_object(self, role: str, object_name: str) -> bool:
        save = getattr(self.command, "save", None)
        if save is None:
            return False
        with NamedTemporaryFile(prefix="structlens_", suffix=".pdb", delete=False) as handle:
            temporary = Path(handle.name)
        try:
            save(str(temporary), object_name)
            self._temporary_paths.append(temporary)
            return self._load_source(role, temporary, pymol_object_name=object_name)
        except Exception:
            return False
        finally:
            temporary.unlink(missing_ok=True)

    def _load_file_into_pymol(self, role: str, path: Path) -> str | None:
        if self.command is None:
            return None
        load = getattr(self.command, "load", None)
        if load is None:
            return None
        object_name = selection_name("panel", path.stem, role)
        try:
            get_names = getattr(self.command, "get_names", None)
            names = {str(name) for name in get_names("objects")} if get_names is not None else set()
            if object_name in names:
                delete = getattr(self.command, "delete", None)
                if delete is not None:
                    delete(object_name)
                else:
                    object_name = selection_name("panel", f"{path.stem}_{path.stat().st_mtime_ns:x}", role)
            load(str(path), object_name)
            return str(object_name)
        except Exception:
            return None

    def _populate_models(self, combo: Any, structure: ProteinStructure) -> None:
        combo.blockSignals(True)
        combo.clear()
        model_ids = tuple(dict.fromkeys(str(chain.model_id) for chain in structure.chains))
        for model_id in model_ids:
            combo.addItem(f"Model {model_id}", model_id)
        combo.blockSignals(False)
        if combo.count():
            combo.setCurrentIndex(0)

    def _populate_chains(
        self,
        combo: Any,
        structure: ProteinStructure,
        *,
        model_id: str | None = None,
    ) -> None:
        combo.blockSignals(True)
        combo.clear()
        for chain in structure.chains:
            if model_id is not None and str(chain.model_id) != model_id:
                continue
            count = len(chain.residue_records or chain.residues)
            combo.addItem(f"Chain {chain.chain_id} · {count} residues", chain.chain_id)
        combo.blockSignals(False)
        if combo.count():
            combo.setCurrentIndex(0)

    def _chain_changed(self, _: int) -> None:
        self._sync_source_model()

    def _model_changed(self, role: str) -> None:
        structure = self.reference_structure if role == "reference" else self.target_structure
        model_combo = self.reference_model_combo if role == "reference" else self.target_model_combo
        chain_combo = self.reference_chain_combo if role == "reference" else self.target_chain_combo
        if structure is not None:
            self._populate_chains(chain_combo, structure, model_id=self._combo_data(model_combo))
        self._sync_source_model()

    def _sync_source_model(self) -> None:
        self.model = self.model.with_sources(
            reference_path=self.reference_edit.text().strip() or None,
            target_path=self.target_edit.text().strip() or None,
            reference_chain_id=self._combo_data(self.reference_chain_combo),
            target_chain_id=self._combo_data(self.target_chain_combo),
        )

    def _combo_data(self, combo: Any) -> str | None:
        return _combo_data(combo)

    # --------------------------------------------------------------- analysis

    def _settings(self) -> AnalysisSettings:
        mode = str(self.mode_combo.currentData() or AlignmentMode.AUTO.value)
        return AnalysisSettings(
            alignment_mode=AlignmentMode(mode),
            minimum_sequence_identity=float(self.identity_spin.value()),
            minimum_sequence_coverage=float(self.coverage_spin.value()),
            refined_rmsd=bool(self.refined_check.isChecked()),
            refinement_cutoff_angstrom=float(self.cutoff_spin.value()),
            usalign_executable=self.usalign_edit.text().strip() or None,
        )

    def _selected_chain(
        self,
        structure: ProteinStructure | None,
        combo: Any,
        model_combo: Any | None = None,
    ) -> ProteinChain | None:
        return _selected_chain_for_gui(structure, combo, model_combo)

    def _start_analysis(self) -> None:
        reference_loaded = self.reference_loaded_source
        target_loaded = self.target_loaded_source
        if reference_loaded is None or target_loaded is None:
            if "build_request" in self._report_controller.__dict__:
                self._start_injected_legacy_request()
                return
            self._show_error("Load one reference and one target chain before comparing.")
            self.nav.setCurrentRow(0)
            return
        if self.reference_edit.text().strip() != str(reference_loaded.path) or self.target_edit.text().strip() != str(
            target_loaded.path
        ):
            self._show_error("A source path was edited after loading; load it before comparing.")
            self.nav.setCurrentRow(0)
            return
        reference = self._selected_chain(
            self.reference_structure, self.reference_chain_combo, self.reference_model_combo
        )
        target = self._selected_chain(self.target_structure, self.target_chain_combo, self.target_model_combo)
        if reference is None or target is None:
            self._show_error("Load one reference and one target chain before comparing.")
            self.nav.setCurrentRow(0)
            return
        try:
            self._report_controller.verify_current_source(reference_loaded.path, reference_loaded.snapshot)
            self._report_controller.verify_current_source(target_loaded.path, target_loaded.snapshot)
        except (OSError, SnapshotError) as exc:
            self._show_error(f"Could not capture verified source input: {exc}")
            self.nav.setCurrentRow(_STRUCTURES_PAGE_INDEX)
            return
        except (TypeError, ValueError) as exc:
            self._show_error(str(exc))
            self.nav.setCurrentRow(_STRUCTURES_PAGE_INDEX)
            return
        try:
            settings = self._settings()
            manual_pairs = (
                tuple(self._manual_pairs(reference, target)) if settings.alignment_mode is AlignmentMode.MANUAL else ()
            )
            request = self._build_snapshot_request(
                reference_loaded,
                target_loaded,
                reference=reference,
                target=target,
                analysis_settings=settings,
                site_definitions=self._site_definitions,
                manual_pairs=manual_pairs,
            )
        except (TypeError, ValueError) as exc:
            self._show_error(str(exc))
            self.nav.setCurrentRow(_STRUCTURES_PAGE_INDEX)
            return
        self._submit_report_request(request)

    def _submit_report_request(self, request: AnalysisReportRequest) -> None:
        self.model = self.model.with_busy("Building canonical analysis report…")
        self._pending_report_request = request
        self._set_busy(True)
        self._cancel_event = Event()

        self._future = _ANALYSIS_EXECUTOR.submit(
            self._report_controller.execute,
            request,
        )
        self._poll_timer = self.c.QTimer(self.widget)
        self._poll_timer.setInterval(50)
        self._poll_timer.timeout.connect(self._poll_analysis)
        self._poll_timer.start()

    def _start_injected_legacy_request(self) -> None:
        """Retain the narrow synthetic request-builder compatibility seam."""

        if not self._load_sources_from_edits():
            self.nav.setCurrentRow(0)
            return
        reference = self._selected_chain(
            self.reference_structure, self.reference_chain_combo, self.reference_model_combo
        )
        target = self._selected_chain(self.target_structure, self.target_chain_combo, self.target_model_combo)
        if reference is None or target is None:
            self._show_error("Load one reference and one target chain before comparing.")
            return
        try:
            settings = self._settings()
            manual_pairs = (
                tuple(self._manual_pairs(reference, target)) if settings.alignment_mode is AlignmentMode.MANUAL else ()
            )
            request = self._report_controller.build_request(
                self.reference_edit.text().strip(),
                self.target_edit.text().strip(),
                reference_chain_id=reference.chain_id,
                target_chain_id=target.chain_id,
                reference_model_id=reference.model_id,
                target_model_id=target.model_id,
                analysis_settings=settings,
                site_definitions=self._site_definitions,
                manual_pairs=manual_pairs,
            )
        except (TypeError, ValueError) as exc:
            self._show_error(str(exc))
            return
        self._submit_report_request(request)

    def _build_snapshot_request(
        self,
        reference_loaded: LoadedSource,
        target_loaded: LoadedSource,
        *,
        reference: ProteinChain,
        target: ProteinChain,
        analysis_settings: AnalysisSettings,
        site_definitions: tuple[SiteDefinition, ...],
        manual_pairs: tuple[tuple[ResidueId, ResidueId], ...],
    ) -> AnalysisReportRequest:
        """Build from loaded snapshots, retaining an old injected builder seam."""

        compatibility_builder = self._report_controller.__dict__.get("build_request")
        if compatibility_builder is not None:
            return compatibility_builder(  # type: ignore[no-any-return]
                reference_loaded.path,
                target_loaded.path,
                reference_chain_id=reference.chain_id,
                target_chain_id=target.chain_id,
                reference_model_id=reference.model_id,
                target_model_id=target.model_id,
                analysis_settings=analysis_settings,
                site_definitions=site_definitions,
                manual_pairs=manual_pairs,
            )
        return self._report_controller.build_request_from_snapshots(
            reference_loaded.snapshot,
            target_loaded.snapshot,
            reference_selection=reference_loaded.selection(
                model_id=str(reference.model_id), chain_id=reference.chain_id
            ),
            target_selection=target_loaded.selection(model_id=str(target.model_id), chain_id=target.chain_id),
            analysis_settings=analysis_settings,
            site_definitions=site_definitions,
            manual_pairs=manual_pairs,
        )

    def _poll_analysis(self) -> None:
        if self._future is None or not self._future.done():
            return
        future = self._future
        self._future = None
        if self._poll_timer is not None:
            self._poll_timer.stop()
        self._poll_timer = None
        try:
            self._analysis_finished(future.result())
        except AnalysisCancelledError:
            self._analysis_cancelled()
        except Exception as exc:
            self._analysis_failed(str(exc))

    def _manual_pairs(self, reference: ProteinChain, target: ProteinChain) -> list[tuple[ResidueId, ResidueId]]:
        pairs: list[tuple[ResidueId, ResidueId]] = []
        for line_number, raw_line in enumerate(self.manual_edit.toPlainText().splitlines(), 1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            pieces = re.split(r"\s*(?:->|=|,)\s*", line, maxsplit=1)
            if len(pieces) != 2:
                raise ValueError(f"Manual pair line {line_number} must use ref -> target syntax")
            reference_id = _find_residue(reference, pieces[0].strip())
            target_id = _find_residue(target, pieces[1].strip())
            if reference_id is None or target_id is None:
                raise ValueError(f"Manual pair line {line_number} names a residue not in the selected chains")
            pairs.append((reference_id, target_id))
        if not pairs:
            raise ValueError("Manual mode needs at least one residue pair")
        return pairs

    def show_report(self, report: AnalysisReport) -> None:
        """ReportView adapter used by the headless report controller."""

        self._analysis_finished(report)

    def show_presentation(self, presentation: ReportPresentation) -> None:
        self._report_presentation = presentation

    def show_status(self, status: str) -> None:
        self._set_status(status)

    def _analysis_finished(self, result: AnalysisResult | AnalysisReport) -> None:
        if isinstance(result, AnalysisReport):
            request = self._pending_report_request or self._report_request
            try:
                self._populate_report(result, request=request)
            except (TypeError, ValueError) as exc:
                self._pending_report_request = None
                self._report_controller.discard_pending()
                self._set_busy(False)
                self._show_error(f"Comparison result was rejected before acceptance: {exc}")
            return
        self._invalidate_analysis_views()
        self._analysis_history = (*self._analysis_history, result)
        self.model = self.model.with_analysis(result)
        self._set_busy(False)
        self._populate_result(result)
        self._set_status(f"Analysis complete · {len(result.correspondences)} aligned positions")
        # The comparison result is immediately visible on the Structures tab;
        # Results remains the compiled cross-analysis view.
        self.nav.setCurrentRow(_STRUCTURES_PAGE_INDEX)

    def _analysis_failed(self, message: str) -> None:
        self._pending_report_request = None
        self._report_controller.discard_pending()
        self.model = self.model.with_error(f"Comparison failed: {message}")
        self._set_busy(False)
        self._show_error(f"Comparison failed: {message}")

    def _analysis_cancelled(self) -> None:
        self._pending_report_request = None
        self._report_controller.discard_pending()
        self.model = self.model.with_status("Comparison cancelled.")
        self._set_busy(False)
        self._set_status("Comparison cancelled; no result was changed")

    def _set_busy(self, busy: bool) -> None:
        self.progress.setVisible(busy)
        self.cancel_button.setVisible(busy)
        self.compare_button.setEnabled(not busy)
        self.run_button.setEnabled(not busy)
        self.compare_button.setText("Comparing…" if busy else "Compare")
        self.run_button.setText("Comparing…" if busy else "Compare")

    def _cancel_analysis(self) -> None:
        self._cancel_event.set()
        cancelled = False
        if self._future is not None:
            cancelled = bool(self._future.cancel())
            if cancelled:
                self._future = None
        if cancelled and self._poll_timer is not None:
            self._poll_timer.stop()
            self._poll_timer = None
        if cancelled:
            self._analysis_cancelled()
        else:
            self.model = self.model.with_status(
                "Cancellation requested; the active report calculation has no cancellation hook."
            )
            self._set_status("Cancellation requested; waiting for the canonical report to finish")

    # -------------------------------------------------------------- rendering


__all__ = ["SourceMixin"]


def _unique_pymol_name(base: str, existing: set[str]) -> str:
    candidate = base
    suffix = 2
    while candidate in existing:
        candidate = f"{base}_{suffix}"
        suffix += 1
    return candidate
