"""Portable PyMOL scenes for a verified canonical report and captured sources."""

from __future__ import annotations

import base64
import json
from importlib.resources import files
from pathlib import Path
from typing import Any

from structlens.application.report_snapshot_io import verify_report_inputs
from structlens.core.models import ResidueId
from structlens.core.parsing import InputSelection, SourceSnapshot
from structlens.core.reports import AnalysisReport
from structlens.core.reports.safe_io import atomic_write_bytes

from .snapshot_selection import selected_coordinates


def _source(
    snapshot: SourceSnapshot, selection: InputSelection, residue: ResidueId | None, role: str, visible: bool
) -> dict[str, Any]:
    state, conformers = selected_coordinates(snapshot, selection)
    chain_id = selection.author_chain_ids[0]
    address = [residue.chain_id, f"{residue.auth_seq_id}{residue.insertion_code or ''}", residue.residue_name] if residue else None
    return {
        "bytes": base64.b64encode(snapshot.decompressed_bytes).decode("ascii"),
        "format": "pdb" if snapshot.logical_format == "pdb" else "cif",
        "state": state, "conformers": conformers,
        "chain": chain_id, "residue": address, "role": role, "visible": visible,
    }


def report_scene_data(
    report: AnalysisReport, snapshots: tuple[SourceSnapshot, SourceSnapshot],
    alignment_index: int | None = None, *, start_focused: bool = True,
    show_reference: bool = True, show_target: bool = True,
) -> dict[str, Any]:
    """Collect immutable scene inputs without resolving mutable source paths."""
    verify_report_inputs(report, snapshots, (), None)
    if len(snapshots) != 2 or report.provenance is None:
        raise ValueError("Both verified coordinate snapshots are required.")
    for role, snapshot, selection in zip(("reference", "target"), snapshots, (report.reference_selection, report.target_selection), strict=True):
        if snapshot.content_id != selection.content_id or snapshot.raw_sha256 != report.provenance.input_hashes[f"{role}_raw"]:
            raise ValueError(f"{role} snapshot does not match its report role")
    analysis = report.analysis
    if analysis is None or analysis.transform is None:
        raise ValueError("This report has no verified target-to-reference fit.")
    item = next((row for row in analysis.correspondences if row.alignment_index == alignment_index), None)
    if alignment_index is not None and item is None:
        raise ValueError("Selected position is absent from this report.")
    sources = [
        _source(snapshots[0], report.reference_selection, item.reference if item else None, "reference", show_reference),
        _source(snapshots[1], report.target_selection, item.target if item else None, "target", show_target),
    ]
    positions = ["gap" if entry["residue"] is None else "{}:{} {}".format(*entry["residue"]) for entry in sources]
    return {
        "report": report.report_id, "rotation": analysis.transform.rotation, "translation": analysis.transform.translation,
        "sources": sources, "focused": bool(item and start_focused),
        "message": "Reference " + positions[0] + " -> Target " + positions[1],
    }


def write_report_scene(
    path: str | Path, report: AnalysisReport, snapshots: tuple[SourceSnapshot, SourceSnapshot],
    alignment_index: int | None = None, *, start_focused: bool = True,
    show_reference: bool = True, show_target: bool = True,
) -> Path:
    """Export a self-contained script with a fixed renderer and encoded user data."""
    data = report_scene_data(report, snapshots, alignment_index, start_focused=start_focused,
                             show_reference=show_reference, show_target=show_target)
    return write_scene_data(path, data)


def write_scene_data(path: str | Path, data: dict[str, Any]) -> Path:
    encoded = base64.b64encode(json.dumps(data, separators=(",", ":")).encode("utf-8")).decode("ascii")
    lines = '\n'.join('"' + encoded[start:start + 120] + '"' for start in range(0, len(encoded), 120))
    renderer = files("structlens.integrations.pymol").joinpath("scene_runtime.py").read_text(encoding="utf-8")
    script = "python\nfrom pymol import cmd\nimport base64, json\ndata = json.loads(base64.b64decode(\n" + lines + "\n))\n" + renderer + "\nrender_scene(cmd, data)\npython end\n"
    dest = Path(path)
    atomic_write_bytes(dest, script.encode("ascii"), max_bytes=512 * 1024 * 1024, label="PyMOL scene", replace_existing=True)
    return dest
