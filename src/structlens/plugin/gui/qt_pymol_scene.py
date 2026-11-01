"""Verified PyMOL handoff for the current report and residue selection."""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path

from structlens.integrations.pymol.launcher import PyMOLLauncher
from structlens.integrations.pymol.scene import report_scene_data, write_scene_data
from structlens.integrations.pymol.scene_runtime import render_scene

from .qt_context import Any, QtMixinContext


class PyMOLSceneMixin(QtMixinContext):
    def _scene_data(self, *, focus: bool) -> dict[str, Any]:
        binding = self._canonical_binding()
        selected = self.model.selection
        index = self._selection_controller.lookup_correspondence().alignment_index if selected is not None else None
        if focus and index is None:
            raise ValueError("Select a residue in the evidence table first.")
        state = self._state_from_controls()
        return report_scene_data(
            binding.report, binding.snapshots, index, start_focused=focus,
            show_reference=state.show_reference, show_target=state.show_target,
        )

    def _write_canonical_scene(self, path: Path, *, focus: bool) -> None:
        write_scene_data(path, self._scene_data(focus=focus))

    def _open_canonical_scene(self, *, focus: bool) -> None:
        try:
            data = self._scene_data(focus=focus)
            key = (data["report"], tuple(s["visible"] for s in data["sources"]))
            if self.command is not None:
                previous = getattr(self, "_last_host_scene", None)
                if not focus and previous and previous[0] == key and previous[1] in self.command.get_scene_list():
                    self.command.scene(previous[1], "recall", animate=0)
                else:
                    scenes = render_scene(self.command, data)
                    self._last_host_scene = (key, scenes[0])
            else:
                launcher = PyMOLLauncher(self.pymol_edit.text().strip() or None)
                launcher.locate()
                standard = self.c.QStandardPaths
                cache = Path(standard.writableLocation(getattr(standard, "StandardLocation", standard).CacheLocation))
                # PyMOL may read the file after this panel closes. Keep one immutable cache entry per view.
                marker = repr((key, [s["residue"] for s in data["sources"]], focus)).encode("utf-8")
                path = cache / "StructLens" / "pymol" / (sha256(marker).hexdigest() + ".pml")
                write_scene_data(path, data)
                launcher.launch_scene(path)
            self._set_status("Opened verified PyMOL scene. Use Scenes to switch between Overview and Residue.")
        except Exception as exc:
            self._show_error(f"Could not open canonical PyMOL scene: {exc}")
