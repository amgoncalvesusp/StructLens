"""Command-proxy adapter for reversible, namespaced StructLens views."""

from __future__ import annotations

from collections.abc import Callable, Iterable

from structlens.core.models import (
    AnalysisResult,
    CorrespondenceStatus,
    ResidueCorrespondence,
    ResidueId,
)
from structlens.plugin.visualization.renderer import VisualizationRenderer, VisualizationState

from .selections import residue_selection, selection_name

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

    def __init__(
        self, command: object | None = None, *, project_id: str = "default"
    ) -> None:
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
        expressions = [
            residue_selection(residue)
            for residue in (item.reference, item.target)
            if residue
        ]
        expression = " or ".join(expressions)
        self._call("select", name, expression)
        self._owned_selections.add(name)
        self._call("zoom", name)
        return name

    def apply(
        self,
        result: AnalysisResult,
        *,
        target_id: str | None = None,
        state: VisualizationState | None = None,
        reference_object: str | None = None,
        target_object: str | None = None,
    ) -> None:
        """Apply a visualization state and remove only previous StructLens state."""

        if self.command is None:
            return
        view = state or VisualizationState()
        target = target_id or result.target_id
        selected = VisualizationRenderer().filtered_correspondences(
            result.correspondences, view
        )
        self.reset()
        if not selected:
            return
        target_name = selection_name(
            self.project_id, target, f"{view.highlight_filter.value}_target"
        )
        target_expression = _join_residues(item.target for item in selected)
        if view.show_target and target_expression:
            target_scope = self._create_view_object(
                target,
                _scoped_expression(target_expression, target_object),
                "target",
            )
            if target_scope is not None:
                self._select(target_name, f"model {target_scope}")
                self._show(target_name, view.representation.value)
                self._owned_selections.add(target_name)
                self._color_target(result, selected, target_name, view, target_scope)
                if view.show_labels:
                    self._call("label", target_name, 'chain + ":" + resi + " " + resn')

        if view.show_reference:
            reference_name = selection_name(
                self.project_id, result.reference_id, f"{view.highlight_filter.value}_reference"
            )
            reference_expression = _join_residues(item.reference for item in selected)
            if reference_expression:
                reference_scope = self._create_view_object(
                    result.reference_id,
                    _scoped_expression(reference_expression, reference_object),
                    "reference",
                )
                if reference_scope is not None:
                    self._select(reference_name, f"model {reference_scope}")
                    self._call("show", "cartoon", reference_name)
                    self._call("color", "gray75", reference_name)
                    self._owned_selections.add(reference_name)

    def _color_target(
        self,
        result: AnalysisResult,
        selected: tuple[ResidueCorrespondence, ...],
        fallback_name: str,
        state: VisualizationState,
        target_scope: str,
    ) -> None:
        if state.color_mode.value == "mutation_status":
            for status in CorrespondenceStatus:
                entries = tuple(item for item in selected if item.status is status)
                expression = _scoped_expression(
                    _join_residues(item.target for item in entries), target_scope
                )
                if expression:
                    name = selection_name(self.project_id, result.target_id, f"{status.value}_target")
                    self._select(name, expression)
                    self._call("color", _STATUS_COLORS[status], name)
                    self._owned_selections.add(name)
            return
        values = [
            _metric_value(item, state.color_mode.value)
            for item in selected
        ]
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

    def _create_view_object(
        self, structure_id: str, expression: str, purpose: str
    ) -> str | None:
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


def _metric_value(item: ResidueCorrespondence, color_mode: str) -> float | None:
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
