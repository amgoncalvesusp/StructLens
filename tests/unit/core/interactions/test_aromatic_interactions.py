from __future__ import annotations

import pytest

from structlens.core.interactions import InteractionType, detect_interactions
from structlens.core.interactions.chemistry import (
    aromatic_ring_atom_names,
    cationic_group_atom_names,
    residue_is_aromatic,
    residue_is_cationic,
)
from structlens.core.models import AtomRecord, ResidueId, ResidueNumbering, ResidueRecord


def _residue(name: str, number: int, atoms: tuple[tuple[str, str, tuple[float, float, float]], ...]) -> ResidueRecord:
    residue_id = ResidueId("ref", "1", "A", str(number), None, name)
    return ResidueRecord(
        residue_id,
        ResidueNumbering(str(number), str(number), None),
        name,
        name[0],
        tuple(AtomRecord(atom, element, xyz) for atom, element, xyz in atoms),
    )


def _phenyl(center_x: float, number: int) -> ResidueRecord:
    atoms = tuple(
        (name, "C", (center_x + x, y, 0.0))
        for name, x, y in (
            ("CG", 0.0, 1.4),
            ("CD1", 1.21, 0.7),
            ("CD2", -1.21, 0.7),
            ("CE1", 1.21, -0.7),
            ("CE2", -1.21, -0.7),
            ("CZ", 0.0, -1.4),
        )
    )
    return _residue("PHE", number, atoms)


def test_aromatic_tables_are_explicit_and_reject_unknown_residues() -> None:
    assert aromatic_ring_atom_names("PHE") == (("CG", "CD1", "CE1", "CZ", "CE2", "CD2"),)
    assert cationic_group_atom_names("ARG") == (("NE", "CZ", "NH1", "NH2"),)
    assert residue_is_aromatic(_phenyl(0.0, 1))
    assert residue_is_cationic(_residue("LYS", 2, (("NZ", "N", (0.0, 0.0, 0.0)),)))
    assert not residue_is_aromatic(_residue("UNK", 3, (("C1", "C", (0.0, 0.0, 0.0)),)))


def test_detects_parallel_pi_stacking_from_ring_centroids_and_normals() -> None:
    records = detect_interactions((_phenyl(0.0, 1), _phenyl(4.0, 2)))

    pi_records = [record for record in records if record.interaction_type is InteractionType.PI_STACKING]
    assert len(pi_records) == 1
    assert pi_records[0].distance_angstrom == pytest.approx(4.0, abs=0.1)
    assert pi_records[0].angle_degrees == pytest.approx(0.0, abs=0.1)
    assert pi_records[0].evidence_mode == "aromatic_ring_geometry"


def test_detects_cation_pi_using_explicit_cation_group_geometry() -> None:
    aromatic = _phenyl(0.0, 1)
    lys = _residue("LYS", 2, (("NZ", "N", (0.0, 0.0, 4.0)),))

    records = detect_interactions((aromatic, lys))

    cation_pi = [record for record in records if record.interaction_type is InteractionType.CATION_PI]
    assert len(cation_pi) == 1
    assert cation_pi[0].distance_angstrom == pytest.approx(4.0, abs=0.1)
    assert cation_pi[0].angle_degrees == pytest.approx(0.0, abs=0.1)
    assert cation_pi[0].evidence_mode == "cation_ring_geometry"


def test_aromatic_interactions_remain_unavailable_when_ring_atoms_are_missing() -> None:
    incomplete = _residue("PHE", 1, (("CG", "C", (0.0, 0.0, 0.0)),))
    lys = _residue("LYS", 2, (("NZ", "N", (0.0, 0.0, 4.0)),))

    records = detect_interactions((incomplete, lys))

    assert all(record.interaction_type is not InteractionType.CATION_PI for record in records)
