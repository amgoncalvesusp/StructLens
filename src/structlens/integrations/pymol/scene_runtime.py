"""Standalone PyMOL code included verbatim in exported PML scenes."""

import base64
import os
import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from typing import Any


@contextmanager
def author_numbering(cmd: Any, is_cif: bool) -> Iterator[None]:
    """Canonical correspondence uses author chain, residue and insertion identifiers."""
    if not is_cif or cmd.get_setting_int("cif_use_auth"):
        yield
    else:
        cmd.set("cif_use_auth", 1)
        try:
            yield
        finally:
            cmd.set("cif_use_auth", 0)


def selected_indices(cmd: Any, raw: str, state: int, chain: str, conformers: list[dict[str, Any]]) -> list[int]:
    """Select exact author-chain atoms and canonical polymer alternate locations."""
    policy = {(c["resi"], c["resn"]): {tuple(a) for a in c["atoms"]} for c in conformers}
    rows: list[tuple[int, str, str, str, str, str]] = []
    cmd.iterate_state(state, raw, "rows.append((index, chain, resi, resn, name, alt))", space={"rows": rows})
    return [
        index for index, chain_id, resi, resn, atom, alt in rows
        if chain_id == chain and ((resi, resn) not in policy or (atom, alt) in policy[resi, resn])
    ]


def render_scene(cmd: Any, data: dict[str, Any]) -> list[str]:
    """Load two isolated chains and record general and optional residue scenes."""
    names = []
    created_scenes = []
    focus_parts = []
    view_before = cmd.get_view()
    previous = getattr(cmd, "_structlens_verified_scene", None)
    old_objects: list[str] = []
    old_enabled = set(cmd.get_names("objects", enabled_only=1))
    old_selections: list[str] = []
    if previous and previous["marker"] in cmd.get_names("selections"):
        # A static atom selection no longer covers an unrelated replacement that reuses the old object's name.
        for name in previous["objects"] + previous["selections"]:
            if name in cmd.get_names("all") and cmd.count_atoms(name + " and " + previous["marker"]):
                (old_objects if name in previous["objects"] else old_selections).append(name)
    marker = cmd.get_unused_name("SL_scene_owner")
    try:
        for source in data["sources"]:
            name = cmd.get_unused_name("SL_" + data["report"][:8] + "_" + source["role"], 0)
            raw = cmd.get_unused_name("SL_scene_source")
            subset = cmd.get_unused_name("SL_scene_atoms")
            try:
                with tempfile.TemporaryDirectory(prefix="structlens_scene_") as directory:
                    filename = os.path.join(directory, "snapshot." + source["format"])
                    with open(filename, "wb") as handle:
                        handle.write(base64.b64decode(source["bytes"]))
                    with author_numbering(cmd, source["format"] == "cif"):
                        cmd.load(filename, raw, format=source["format"], zoom=0)
                indices = selected_indices(cmd, raw, source["state"], source["chain"], source["conformers"])
                if not indices:
                    raise ValueError("Selected chain has no atoms in the specified model")
                cmd.select_list(subset, raw, indices, state=source["state"], mode="index")
                names.append(name)
                cmd.create(name, subset, source["state"], 1, zoom=0)
            finally:
                cmd.delete(subset)
                cmd.delete(raw)
            if source["role"] == "target":
                coords = cmd.get_coords(name, state=1)
                cmd.load_coords(coords @ data["rotation"] + data["translation"], name, state=1)
            cmd.hide("everything", name)
            cmd.show("cartoon", name)
            cmd.show("sticks", name + " and not hydro and not polymer.protein")
            cmd.color("marine" if source["role"] == "reference" else "orange", name)
            if source["residue"] is not None:
                atoms: list[tuple[int, str, str, str]] = []
                cmd.iterate(name, "atoms.append((index, chain, resi, resn))", space={"atoms": atoms})
                indices = [i for i, chain, resi, resn in atoms if [chain, resi, resn] == source["residue"]]
                if not indices:
                    raise ValueError("Selected report residue is absent from its model and chain")
                part = cmd.get_unused_name("SL_residue")
                focus_parts.append(part)
                cmd.select_list(part, name, indices, state=1, mode="index")
                cmd.disable(part)
            if not source["visible"]:
                cmd.disable(name)
        scope = " or ".join(names)
        cmd.select(marker, scope, enable=0)
        for name in old_objects:
            cmd.disable(name)
        cmd.orient(scope)
        stem = "SL_" + data["report"][:8]
        used = set(cmd.get_scene_list())
        key = stem
        serial = 1
        while key + "_Overview" in used or key + "_Residue" in used:
            serial += 1
            key = stem + "_" + str(serial)
        overview = key + "_Overview"
        cmd.scene(overview, "store", message="StructLens | " + data["report"] +
                  " | target fitted to reference; polymer conformers follow the report; other components as deposited")
        created_scenes.append(overview)
        if focus_parts:
            focus_scope = " or ".join(focus_parts)
            cmd.show("sticks", focus_scope)
            cmd.zoom(focus_scope, buffer=4)
            detail = key + "_Residue"
            cmd.scene(detail, "store", message=data["message"])
            created_scenes.append(detail)
        cmd.scene(created_scenes[-1] if data["focused"] else overview, "recall", animate=0)
    except Exception:
        for scene in created_scenes:
            cmd.scene(scene, "clear")
        for name in focus_parts + names + [marker]:
            cmd.delete(name)
        for name in old_objects:
            if name in old_enabled:
                cmd.enable(name)
        cmd.set_view(view_before)
        raise
    if previous:
        for scene, message in previous["scenes"].items():
            if scene in cmd.get_scene_list() and cmd.get_scene_message(scene) == message:
                cmd.scene(scene, "clear")
        for name in old_selections + old_objects:
            cmd.delete(name)
        if old_objects:
            cmd.delete(previous["marker"])
    cmd._structlens_verified_scene = {
        "objects": names, "selections": focus_parts, "marker": marker,
        "scenes": {scene: cmd.get_scene_message(scene) for scene in created_scenes},
    }
    return created_scenes
