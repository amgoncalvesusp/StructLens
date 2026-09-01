"""Command-proxy adapter for reversible, namespaced StructLens views."""

from __future__ import annotations

from collections.abc import Callable, Iterable
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any, cast

from structlens.application.report_snapshot_io import verify_report_inputs
from structlens.core.models import (
    AnalysisResult,
    CorrespondenceStatus,
    ResidueCorrespondence,
    ResidueId,
)
from structlens.core.parsing import InputSelection, SourceSnapshot, load_structure_legacy_snapshot
from structlens.core.reports import AnalysisReport, AnalysisSnapshot, CorrespondenceSnapshot
from structlens.plugin.visualization.renderer import VisualizationRenderer, VisualizationState

from .selections import quote_identifier, residue_selection, selection_name

_STATUS_COLORS = {
    CorrespondenceStatus.CONSERVED: "slate",
    CorrespondenceStatus.SUBSTITUTION: "orange",
    CorrespondenceStatus.INSERTION: "yellow",
    CorrespondenceStatus.DELETION: "red",
    CorrespondenceStatus.NONSTANDARD: "magenta",
    CorrespondenceStatus.UNMAPPED: "gray70",
}
_REPRESENTATIONS = {
    "sticks": ("sticks",),
    "spheres": ("spheres",),
    "sticks_spheres": ("sticks", "spheres"),
    "cartoon_sticks": ("cartoon", "sticks"),
    "surface": ("surface",),
    "labels": (),
}


class PyMOLAdapter:
    """Render a result through an injected PyMOL ``cmd``-compatible proxy."""

    def __init__(self, command: object | None = None, *, project_id: str = "default") -> None:
        self.command = command
        self.project_id = project_id
        self._owned_selections: set[str] = set()
        self._owned_objects: set[str] = set()
        self._previous_view: tuple[object, ...] | None = None

    def focus_residue(self, item: ResidueCorrespondence, target_id: str) -> str:
        """Select and zoom one correspondence without changing unrelated objects."""

        if self.command is None:
            return ""
        name = selection_name(self.project_id, target_id, "focus")
        get_view = getattr(self.command, "get_view", None)
        if get_view is not None:
            view = get_view()
            if isinstance(view, (tuple, list)):
                self._previous_view = tuple(view)
        expressions = [residue_selection(residue) for residue in (item.reference, item.target) if residue]
        expression = " or ".join(expressions)
        self._call("select", name, expression)
        self._owned_selections.add(name)
        self._call("zoom", name)
        return name

    def apply(
        self,
        result: AnalysisResult | AnalysisReport,
        *,
        target_id: str | None = None,
        state: VisualizationState | None = None,
        reference_object: str | None = None,
        target_object: str | None = None,
    ) -> None:
        """Apply a legacy/render-only view to already loaded objects."""

        if self.command is None:
            return
        view = state or VisualizationState()
        analysis = result.analysis if isinstance(result, AnalysisReport) else result
        if analysis is None:
            self.reset()
            return
        self.reset()
        self._apply_analysis(
            analysis,
            target_id=target_id,
            state=view,
            reference_object=reference_object,
            target_object=target_object,
        )

    def apply_report(
        self,
        report: AnalysisReport,
        *,
        snapshots: tuple[SourceSnapshot, SourceSnapshot],
        state: VisualizationState | None = None,
        reference_object: str | None = None,
        target_object: str | None = None,
    ) -> None:
        """Apply a canonical view from exact snapshots and selected chains.

        The object-name arguments are capture/UI hints only.  They are never
        used as canonical coordinate authority.
        """

        del reference_object, target_object
        if self.command is None:
            return
        if not isinstance(report, AnalysisReport):
            raise TypeError("report must be an AnalysisReport")
        if len(snapshots) != 2 or any(not isinstance(item, SourceSnapshot) for item in snapshots):
            raise TypeError("snapshots must contain exact reference and target SourceSnapshot values")
        verify_report_inputs(report, snapshots, (), None)
        if report.analysis is None:
            self.reset()
            return
        view = state or VisualizationState()
        self.reset()
        try:
            reference_scope = self._materialize_selected_snapshot(
                snapshots[0], report.reference_selection, report.report_id, "reference"
            )
            target_scope = self._materialize_selected_snapshot(
                snapshots[1], report.target_selection, report.report_id, "target"
            )
            self._apply_analysis(
                report.analysis,
                target_id=None,
                state=view,
                reference_object=reference_scope,
                target_object=target_scope,
            )
        except Exception:
            self.reset()
            raise

    def _apply_analysis(
        self,
        analysis: AnalysisResult | AnalysisSnapshot,
        *,
        target_id: str | None,
        state: VisualizationState,
        reference_object: str | None,
        target_object: str | None,
    ) -> None:
        """Render an analysis after its coordinate scopes are established."""

        target = target_id or analysis.target_id
        renderer = VisualizationRenderer()
        selected: tuple[ResidueCorrespondence | CorrespondenceSnapshot, ...]
        if isinstance(analysis, AnalysisResult):
            selected = renderer.filtered_correspondences(analysis.correspondences, state)
        else:
            selected = cast(
                tuple[CorrespondenceSnapshot, ...],
                renderer.filtered_correspondences(cast(Any, analysis.correspondences), state),
            )
        if not selected:
            return
        target_name = selection_name(self.project_id, target, f"{state.highlight_filter.value}_target")
        target_expression = _join_residues(item.target for item in selected)
        if state.show_target and target_expression:
            target_scope = self._create_view_object(
                target,
                _scoped_expression(target_expression, target_object),
                "target",
            )
            if target_scope is not None:
                self._select(target_name, f"model {target_scope}")
                self._show(target_name, state.representation.value)
                self._owned_selections.add(target_name)
                self._color_target(analysis, selected, target_name, state, target_scope)
                if state.show_labels:
                    self._call("label", target_name, 'chain + ":" + resi + " " + resn')

        if state.show_reference:
            reference_name = selection_name(
                self.project_id, analysis.reference_id, f"{state.highlight_filter.value}_reference"
            )
            reference_expression = _join_residues(item.reference for item in selected)
            if reference_expression:
                reference_scope = self._create_view_object(
                    analysis.reference_id,
                    _scoped_expression(reference_expression, reference_object),
                    "reference",
                )
                if reference_scope is not None:
                    self._select(reference_name, f"model {reference_scope}")
                    self._call("show", "cartoon", reference_name)
                    self._call("color", "gray75", reference_name)
                    self._owned_selections.add(reference_name)

    def _materialize_selected_snapshot(
        self,
        snapshot: SourceSnapshot,
        selection: InputSelection,
        report_id: str,
        role: str,
    ) -> str:
        """Load exact bytes, then isolate the selected model state and chain."""

        if snapshot.content_id != selection.content_id:
            raise ValueError(f"{role} snapshot does not match the report selection")
        if len(selection.author_chain_ids) != 1:
            raise ValueError(f"{role} PyMOL visualization requires exactly one author chain")
        chain_id = selection.author_chain_ids[0]
        parsed = load_structure_legacy_snapshot(snapshot)
        model_ids = tuple(dict.fromkeys(str(chain.model_id) for chain in parsed.chains))
        if selection.model_id not in model_ids:
            raise ValueError(f"{role} selected model {selection.model_id} is absent from the snapshot")
        if not any(str(chain.model_id) == selection.model_id and chain.chain_id == chain_id for chain in parsed.chains):
            raise ValueError(
                f"{role} selected model {selection.model_id} / chain {chain_id} is absent from the snapshot"
            )
        source_state = model_ids.index(selection.model_id) + 1
        source_name = self._unique_object_name(report_id, f"{role}_snapshot")
        suffix = ".pdb" if snapshot.logical_format == "pdb" else ".cif"
        temporary: Path | None = None
        try:
            with NamedTemporaryFile(prefix="structlens_pymol_", suffix=suffix, delete=False) as handle:
                handle.write(snapshot.decompressed_bytes)
                temporary = Path(handle.name)
            load = getattr(self.command, "load", None)
            if load is None:
                raise RuntimeError("PyMOL command proxy cannot load canonical snapshots")
            load(str(temporary), source_name)
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
        self._owned_objects.add(source_name)

        selected_name = self._unique_object_name(report_id, f"{role}_model_{selection.model_id}_chain_{chain_id}")
        create = getattr(self.command, "create", None)
        if create is None:
            raise RuntimeError("PyMOL command proxy cannot isolate canonical model states")
        # PyMOL MODEL records are object states, not an atom-selection
        # ``model <number>`` predicate.  ``source_state`` performs the real
        # model isolation; the expression constrains the authoritative chain.
        create(
            selected_name,
            f"model {quote_identifier(source_name)} and chain {quote_identifier(chain_id)}",
            source_state,
            1,
        )
        self._owned_objects.add(selected_name)
        return selected_name

    def _unique_object_name(self, report_id: str, purpose: str) -> str:
        base = selection_name(self.project_id, report_id[:12], purpose)
        get_names = getattr(self.command, "get_names", None)
        names = {str(name) for name in get_names("objects")} if get_names is not None else set()
        names.update(self._owned_objects)
        candidate = base
        suffix = 2
        while candidate in names:
            candidate = f"{base}_{suffix}"
            suffix += 1
        return candidate

    def _color_target(
        self,
        result: AnalysisResult | AnalysisSnapshot,
        selected: tuple[ResidueCorrespondence | CorrespondenceSnapshot, ...],
        fallback_name: str,
        state: VisualizationState,
        target_scope: str,
    ) -> None:
        if state.color_mode.value == "mutation_status":
            for status in CorrespondenceStatus:
                entries = tuple(item for item in selected if item.status is status)
                expression = _scoped_expression(_join_residues(item.target for item in entries), target_scope)
                if expression:
                    name = selection_name(self.project_id, result.target_id, f"{status.value}_target")
                    self._select(name, expression)
                    self._call("color", _STATUS_COLORS[status], name)
                    self._owned_selections.add(name)
            return
        values = [_metric_value(item, state.color_mode.value) for item in selected]
        finite = [value for value in values if value is not None]
        minimum = min(finite) if finite else 0.0
        maximum = max(finite) if finite else 0.0
        for item, value in zip(selected, values, strict=True):
            if item.target is None:
                continue
            name = selection_name(
                self.project_id,
                result.target_id,
                f"row_{item.alignment_index}",
            )
            self._select(
                name,
                _scoped_expression(residue_selection(item.target), target_scope),
            )
            self._call("color", _scale_color(value, minimum, maximum), name)
            self._owned_selections.add(name)
        if not finite:
            self._call("color", "slate", fallback_name)

    def _create_view_object(self, structure_id: str, expression: str, purpose: str) -> str | None:
        """Copy selected atoms into an owned object before changing visuals."""

        create = getattr(self.command, "create", None)
        if create is None:
            return None
        name = selection_name(self.project_id, structure_id, f"{purpose}_view")
        create(name, expression)
        self._owned_objects.add(name)
        return name

    def reset(self) -> None:
        """Delete only names created by this adapter instance."""

        if self._previous_view is not None:
            self._call("set_view", self._previous_view)
            self._previous_view = None
        for name in tuple(self._owned_selections):
            self._call("delete", name)
        self._owned_selections.clear()
        for name in tuple(self._owned_objects):
            self._call("delete", name)
        self._owned_objects.clear()

    def _select(self, name: str, expression: str) -> None:
        self._call("select", name, expression)

    def _show(self, name: str, representation: str) -> None:
        for mode in _REPRESENTATIONS.get(representation, ("sticks",)):
            self._call("show", mode, name)

    def _call(self, method: str, *args: object) -> None:
        function: Callable[..., object] | None = getattr(self.command, method, None)
        if function is not None:
            function(*args)


def _join_residues(residues: Iterable[ResidueId | None]) -> str:
    expressions = [residue_selection(residue) for residue in residues if residue is not None]
    return " or ".join(expressions)


def _scoped_expression(expression: str, object_name: str | None) -> str:
    if object_name is None:
        return expression
    return f"model {object_name} and ({expression})"


def _metric_value(item: ResidueCorrespondence | CorrespondenceSnapshot, color_mode: str) -> float | None:
    if color_mode == "ca_displacement":
        return item.ca_displacement_angstrom
    if color_mode == "backbone_rmsd":
        return item.backbone_rmsd_angstrom
    return None


def _scale_color(value: float | None, minimum: float, maximum: float) -> str:
    if value is None:
        return "slate"
    if maximum <= minimum:
        return "yellow"
    ratio = (value - minimum) / (maximum - minimum)
    if ratio < 0.34:
        return "blue"
    if ratio < 0.67:
        return "yellow"
    return "red"


__all__ = ["PyMOLAdapter"]
