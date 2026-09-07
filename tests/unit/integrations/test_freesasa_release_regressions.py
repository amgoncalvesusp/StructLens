from pathlib import Path

import pytest

from structlens.integrations.freesasa.adapter import FreeSASAAdapter, calculate_sasa

PDB = "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \nEND\n"


def test_sasa_from_file_uses_real_freesasa(tmp_path: Path) -> None:
    path = tmp_path / "single_atom.pdb"
    path.write_text(PDB, encoding="ascii")
    area = calculate_sasa(path)
    assert area is not None and area > 0
    assert path.read_text(encoding="ascii") == PDB


def test_pdb_text_and_file_calculate_the_same_area(tmp_path: Path) -> None:
    path = tmp_path / "single_atom.pdb"
    path.write_text(PDB, encoding="ascii")
    adapter = FreeSASAAdapter()
    assert adapter.calculate_pdb(PDB) == pytest.approx(adapter.calculate_file(path))


def test_missing_structure_reports_unavailable(tmp_path: Path) -> None:
    assert calculate_sasa(tmp_path / "missing.pdb") is None
