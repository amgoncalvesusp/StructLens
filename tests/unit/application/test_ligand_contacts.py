"""Coordinate contacts and DockLens result interoperability."""

import json
from hashlib import sha256
from pathlib import Path

import pytest

from structlens.application.ligand_contacts import import_docklens_contacts, ligand_proximity
from structlens.application.report_input import canonical_selection, canonical_snapshot
from structlens.core.parsing import InputSelection, StructureFormat, capture_snapshot, load_structure_evidence

FIXTURE = Path(__file__).parents[2] / "fixtures" / "parsing" / "enriched.pdb"


def reference():
    source = capture_snapshot(FIXTURE)
    selection = InputSelection(
        source.content_id, source.display_name, StructureFormat.PDB, "1", author_chain_ids=("A",)
    )
    return load_structure_evidence(
        canonical_snapshot(source, "reference"), selection=canonical_selection(selection, "reference")
    )


def test_proximity_uses_exact_ligand_heavy_atoms_and_preserves_insertion_codes():
    parsed = reference()
    ligand = next(c for c in parsed.components if c.residue_name == "ATP")
    neighborhood = ligand_proximity(parsed, ligand.component_id)
    rows = neighborhood.contacts
    assert neighborhood.outside == 1
    assert neighborhood.unavailable == 0
    assert [(r.residue, r.distance) for r in rows if r.residue == "GLY100AA"] == [("GLY100AA", 4.0)]
    assert all(r.family == "Heavy-atom proximity" and r.distance <= 4.5 for r in rows)
    assert not any(r.residue == "ALA100A" for r in rows)
    with pytest.raises(ValueError):
        ligand_proximity(parsed, "ATP")


def document(digest):
    return {
        "schema": "docklens-result-v1",
        "run_id": "run-1",
        "generated_by": {"application": "DockLens", "version": "1.4.0"},
        "inputs": {"receptor_sha256": digest, "poses_sha256": sha256(b"pose").hexdigest()},
        "parameters": {"hbond_preset": "luna_dsv", "interaction_types": ["hbond"], "counting_unit": "semantic"},
        "totals": {"poses": 2, "interactions": 1, "errors": 0},
        "poses": [
            {"pose_id": "p1", "ligand_id": "ATP", "total_interactions": 1},
            {"pose_id": "p2", "ligand_id": "LIG", "total_interactions": 0},
        ],
        "interactions": [
            {
                "pose_id": "p1",
                "type": "hbond",
                "receptor_residue": "GLY100AA",
                "receptor_atom": "N",
                "ligand_atom": "O1",
                "distance_A": 3.2,
                "chemistry_basis": "explicit H",
                "chemistry_confidence": "supported",
            }
        ],
        "input_qc": [],
    }


def test_import_matches_both_hashes_preserves_pose_and_evidence(tmp_path):
    digest = "b" * 64
    path = tmp_path / "result.json"
    data = document(digest)
    path.write_text(json.dumps(data))
    poses = tmp_path / "poses.pdbqt"
    poses.write_bytes(b"pose")
    imported = import_docklens_contacts(path, digest, poses)
    assert imported.profile == "luna_dsv"
    assert imported.poses[0].rows[0].family == "hbond"
    assert "explicit H" in imported.poses[0].rows[0].evidence
    assert imported.poses[1].rows == ()
    with pytest.raises(ValueError, match="receptor"):
        import_docklens_contacts(path, "a" * 64, poses)
    poses.write_bytes(b"other pose")
    with pytest.raises(ValueError, match="poses"):
        import_docklens_contacts(path, digest, poses)
    poses.write_bytes(b"pose")
    for invalid in (-3, 10**400):
        data["interactions"][0]["distance_A"] = invalid
        path.write_text(json.dumps(data))
        with pytest.raises(ValueError, match="distance"):
            import_docklens_contacts(path, digest, poses)


def test_docklens_import_rejects_nonregular_local_paths(tmp_path):
    with pytest.raises(ValueError, match="regular"):
        import_docklens_contacts(tmp_path, "b" * 64, tmp_path)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("inputs", []),
        ("generated_by", {}),
        ("parameters", {}),
        (
            "poses",
            [
                {"pose_id": "p1", "ligand_id": "ATP", "total_interactions": 1},
                {"pose_id": "p1", "ligand_id": "ATP", "total_interactions": 0},
            ],
        ),
        ("totals", {"poses": 0}),
    ],
)
def test_malformed_docklens_result_reports_actionable_error(tmp_path, field, value):
    data = document("b" * 64)
    data[field] = value
    path = tmp_path / "result.json"
    path.write_text(json.dumps(data))
    poses = tmp_path / "poses.pdbqt"
    poses.write_bytes(b"pose")
    with pytest.raises(ValueError):
        import_docklens_contacts(path, "b" * 64, poses)
