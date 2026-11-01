"""Verified, standalone PyMOL scene handoff."""

import base64
import json
from dataclasses import replace
from pathlib import Path

import pytest
from tests.unit.plugin._task13_round4_fixtures import rich_report


def scene_data(path):
    text = path.read_text(encoding="ascii")
    encoded = text.split("data = json.loads(base64.b64decode(\n", 1)[1].split("\n))", 1)[0]
    return json.loads(base64.b64decode("".join(line.strip('"') for line in encoded.splitlines())))


def test_scene_uses_captured_sources_and_report_transform(tmp_path):
    from structlens.integrations.pymol.scene import write_report_scene

    request, report = rich_report(tmp_path)
    output = tmp_path / "scene [review];.pml"
    for selection in (report.reference_selection, report.target_selection):
        if selection.path:
            Path(selection.path).write_text("replacement after comparison")
    item = next(x for x in report.analysis.correspondences if x.reference and x.target)
    write_report_scene(output, report, (request.reference_snapshot, request.target_snapshot), item.alignment_index)
    data = scene_data(output)
    assert data["report"] == report.report_id
    assert data["rotation"] == [list(row) for row in report.analysis.transform.rotation]
    assert data["translation"] == list(report.analysis.transform.translation)
    assert data["focused"]
    assert [base64.b64decode(source["bytes"]) for source in data["sources"]] == [
        request.reference_snapshot.decompressed_bytes, request.target_snapshot.decompressed_bytes
    ]
    assert data["sources"][0]["residue"] == [item.reference.chain_id, str(item.reference.auth_seq_id) + (item.reference.insertion_code or ""), item.reference.residue_name]
    assert data["sources"][1]["residue"] == [item.target.chain_id, str(item.target.auth_seq_id) + (item.target.insertion_code or ""), item.target.residue_name]


def test_missing_fit_stale_snapshot_or_unmapped_index_cannot_overwrite_scene(tmp_path):
    from structlens.integrations.pymol.scene import write_report_scene

    request, report = rich_report(tmp_path)
    dest = tmp_path / "existing.pml"
    dest.write_text("preserve")
    from structlens.core.parsing import capture_snapshot

    other = tmp_path / "other.pdb"
    other.write_bytes(request.reference_snapshot.decompressed_bytes + b"REMARK unrelated source\n")
    with pytest.raises(ValueError):
        write_report_scene(dest, report, (capture_snapshot(other), request.target_snapshot))
    with pytest.raises(ValueError, match="position"):
        write_report_scene(dest, report, (request.reference_snapshot, request.target_snapshot), 999)
    no_fit = replace(report, analysis=replace(report.analysis, transform=None))
    with pytest.raises(ValueError, match="fit"):
        write_report_scene(dest, no_fit, (request.reference_snapshot, request.target_snapshot))
    assert dest.read_text() == "preserve"


def test_launcher_passes_scene_as_single_path(tmp_path, monkeypatch):
    from structlens.integrations.pymol.launcher import PyMOLLauncher

    scene = tmp_path / 'a; space.pml'
    scene.write_text("scene")
    calls = []
    monkeypatch.setattr("structlens.integrations.pymol.launcher.subprocess.Popen",
                        lambda args, **kwargs: calls.append((args, kwargs)) or type("Process", (), {"pid": 41})())
    PyMOLLauncher(scene).launch_scene(scene)
    assert calls == [([str(scene), str(scene.resolve())], {"shell": False})]


def test_applied_host_view_uses_complete_homogeneous_matrix():
    import numpy as np

    from structlens.integrations.pymol.transforms import apply_transform

    class Command:
        def transform_object(self, name, values, state, log, selection, homogenous):
            self.call = name, values, state, homogenous

    cmd = Command()
    rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    apply_transform(cmd, "target", rotation, (2, 3, 4))
    name, flat, state, homogeneous = cmd.call
    assert name == "target" and state == homogeneous == 1
    pos = np.array([1, 2, 3])
    matrix = np.array(flat).reshape(4, 4)
    np.testing.assert_allclose((matrix @ [*pos, 1])[:3], pos @ rotation + [2, 3, 4])


def test_scene_distinguishes_state_chain_and_one_sided_residues(tmp_path):
    from tests.unit.plugin._task13_round5_fixtures import multi_selected_report

    from structlens.core.models import CorrespondenceStatus
    from structlens.integrations.pymol.scene import report_scene_data

    request, report = multi_selected_report(tmp_path)
    snapshots = (request.reference_snapshot, request.target_snapshot)
    mapped = report.analysis.correspondences[0]
    for side in ("reference", "target"):
        changed = replace(mapped, **{side: None}, status=CorrespondenceStatus.INSERTION if side == "reference" else CorrespondenceStatus.DELETION)
        example = replace(report, analysis=replace(report.analysis, correspondences=(changed,)))
        sources = report_scene_data(example, snapshots, changed.alignment_index)["sources"]
        assert [(s["chain"], s["state"]) for s in sources] == [("B", 2), ("B", 1)]
        assert (sources[0]["residue"] is None) == (side == "reference")
        assert (sources[1]["residue"] is None) == (side == "target")


@pytest.mark.parametrize("suffix", ["pdb", "mmcif"])
def test_state_ordinal_includes_hetero_only_models(tmp_path, suffix):
    from structlens.core.parsing import InputSelection, capture_snapshot
    from structlens.integrations.pymol.scene import _source

    if suffix == "pdb":
        source = (
            "MODEL        8\nHETATM    1  O   HOH A  10       0.000   0.000   0.000  1.00 15.00           O\nENDMDL\n"
            "MODEL       24\nATOM      2  CA  ALA A   1       9.000   1.000   2.000  1.00 15.00           C\nENDMDL\nEND\n"
        )
    else:
        header = Path("tests/fixtures/parsing/enriched.mmcif").read_text().split("ATOM 1 ")[0]
        source = header + "HETATM 1 O O . HOH X 2 . ? 0.0 0.0 0.0 1.0 15.0 ? 10 HOH A O 8\nATOM 2 C CA . ALA X 1 1 ? 9.0 1.0 2.0 1.0 15.0 ? 1 ALA A CA 24\n"
    path = tmp_path / f"water-first.{suffix}"
    path.write_text(source, encoding="ascii")
    snapshot = capture_snapshot(path)
    selection = InputSelection(snapshot.content_id, snapshot.display_name, snapshot.logical_format, "24", author_chain_ids=("A",))
    entry = _source(snapshot, selection, None, "target", True)
    assert entry["state"] == 2
