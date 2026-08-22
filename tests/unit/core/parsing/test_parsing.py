"""Contract tests for normalized PDB, mmCIF, FASTA, and gzip inputs."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from structlens.core.parsing import load_fasta, load_structure

_FIXTURES = Path(__file__).parents[3] / "fixtures" / "parsing"


def test_pdb_keeps_author_numbering_insertion_codes_altlocs_and_missing_atoms() -> None:
    structure = load_structure(_FIXTURES / "numbering_altloc.pdb")
    chain = structure.chains[0]
    first, inserted, nonstandard = chain.residue_records

    assert chain.sequence == "AGX"
    assert [residue.auth_seq_id for residue in chain.residues] == ["100", "100", "101"]
    assert inserted.residue_id.insertion_code == "A"
    assert [atom.name for atom in inserted.atoms] == ["CA"]
    assert first.atoms[1].name == "CA"
    assert first.atoms[1].coordinate == pytest.approx((12.0, 10.0, 10.0))
    assert nonstandard.is_standard is False
    assert nonstandard.one_letter is None


def test_mmcif_retains_author_and_label_numbering_without_changing_order() -> None:
    structure = load_structure(_FIXTURES / "numbering.mmcif")
    chain = structure.chains[0]

    assert chain.sequence == "AG"
    assert [record.numbering.auth_seq_id for record in chain.residue_records] == [
        "100",
        "100",
    ]
    assert [record.numbering.label_seq_id for record in chain.residue_records] == [
        "42",
        "43",
    ]
    assert chain.residue_records[1].numbering.insertion_code == "A"


def test_load_structure_accepts_gzip_compressed_pdb(tmp_path: Path) -> None:
    source = _FIXTURES / "numbering_altloc.pdb"
    compressed = tmp_path / "numbering_altloc.pdb.gz"
    with source.open("rb") as input_file, gzip.open(compressed, "wb") as output_file:
        output_file.write(input_file.read())

    structure = load_structure(compressed)

    assert structure.chains[0].sequence == "AGX"
    assert structure.source_path == str(compressed)


def test_load_fasta_reads_each_record_without_concatenating_them() -> None:
    sequences = load_fasta(_FIXTURES / "sequences.fasta")

    assert [(record.identifier, record.sequence) for record in sequences] == [
        ("reference", "ACDE"),
        ("target", "ATDE"),
    ]
