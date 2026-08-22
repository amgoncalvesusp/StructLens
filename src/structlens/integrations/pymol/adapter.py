"""Small command-proxy adapter that never imports PyMOL at module import time."""

from __future__ import annotations

from collections.abc import Callable

from structlens.core.models import AnalysisResult, ResidueCorrespondence

from .selections import residue_selection, selection_name


class PyMOLAdapter:
    def __init__(
        self, command: object | None = None, *, project_id: str = "default"
    ) -> None:
        self.command = command
        self.project_id = project_id
        self._owned_selections: set[str] = set()

    def focus_residue(self, item: ResidueCorrespondence, target_id: str) -> str:
        if self.command is None:
            return ""
        parts = [selection_name(self.project_id, target_id, "focus")]
        expressions = [
            residue_selection(residue)
            for residue in (item.reference, item.target)
            if residue
        ]
        expression = " or ".join(expressions)
        name = parts[0]
        self._call("select", name, expression)
        self._owned_selections.add(name)
        self._call("zoom", name)
        return name

    def apply(self, result: AnalysisResult, *, target_id: str | None = None) -> None:
        if self.command is None:
            return
        target = target_id or result.target_id
        name = selection_name(self.project_id, target, "mutations")
        expressions = [
            residue_selection(item.target)
            for item in result.correspondences
            if item.target is not None
            and item.status.value not in {"conserved", "unmapped"}
        ]
        if expressions:
            self._call("select", name, " or ".join(expressions))
            self._owned_selections.add(name)
            self._call("show", "sticks", name)

    def reset(self) -> None:
        for name in tuple(self._owned_selections):
            self._call("delete", name)
        self._owned_selections.clear()

    def _call(self, method: str, *args: object) -> None:
        function: Callable[..., object] | None = getattr(self.command, method, None)
        if function is not None:
            function(*args)


__all__ = ["PyMOLAdapter"]
