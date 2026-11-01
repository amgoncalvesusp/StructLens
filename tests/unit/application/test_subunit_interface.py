"""Inter-chain measurements in the deposited coordinate frame."""

import pytest

from structlens.application.subunit_interface import analyze_subunit_interface
from structlens.core.parsing import capture_snapshot


def pdb_atom(serial, name, res, chain, pos, x, element, *, occ=1, insertion=" ", alt=" "):
    return (
        f"ATOM  {serial:5d} {name:>4s}{alt}{res:3s} {chain}{pos:4d}{insertion}   "
        f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}{occ:6.2f}{15.0:6.2f}          {element:>2s}\n"
    )


@pytest.fixture
def complex_snapshot(tmp_path):
    lines = [
        pdb_atom(1, "NZ", "LYS", "A", 10, 0, "N", insertion="B"),
        pdb_atom(2, "NZ", "LYS", "A", 11, 8, "N"),
        pdb_atom(3, "CA", "GLY", "A", 12, 3, "C", occ=0),
        pdb_atom(4, "OD1", "ASP", "B", 30, 3, "O"),
        pdb_atom(5, "OD2", "ASP", "B", 30, 4, "O"),
        pdb_atom(6, "H", "ALA", "B", 31, 0, "H"),
        pdb_atom(7, "CA", "ALA", "B", 31, 30, "C"),
        pdb_atom(8, "CA", "ALA", "B", 32, 2, "C", occ=0.1, alt="A"),
        pdb_atom(9, "CA", "ALA", "B", 32, 20, "C", occ=0.9, alt="B"),
    ]
    path = tmp_path / "complex.pdb"
    path.write_text(
        "MODEL        1\n"
        + "".join(lines)
        + "ENDMDL\nMODEL        2\n"
        + pdb_atom(1, "CA", "ALA", "A", 1, 0, "C")
        + pdb_atom(2, "CA", "ALA", "B", 1, 30, "C")
        + "ENDMDL\nEND\n"
    )
    return capture_snapshot(path)


def test_pairs_are_unique_model_specific_and_observed(complex_snapshot):
    result = analyze_subunit_interface(complex_snapshot, "1", "A", "B")
    assert result.source_sha256 == complex_snapshot.raw_sha256
    assert [
        (p.residue_a.auth_seq_id, p.residue_b.auth_seq_id, p.distance_angstrom, p.atom_pair_count)
        for p in result.contacts
    ] == [("10", "30", 3.0, 2), ("11", "30", 4.0, 2)]
    first = result.contacts[0]
    assert first.residue_a.insertion_code == "B" and first.atom_a == "NZ" and first.atom_b == "OD1"
    salt = next(i for i in first.interactions if i.interaction_type.value == "salt_bridge")
    assert (salt.atom_a, salt.atom_b, salt.distance_angstrom) == ("NZ", "OD1", 3.0)
    assert (result.observed_a, result.unavailable_a, result.observed_b, result.unavailable_b) == (2, 1, 3, 0)
    assert [(r.residue_id.chain_id, r.partner_count, r.atom_pair_count) for r in result.residues] == [
        ("B", 2, 4),
        ("A", 1, 2),
        ("A", 1, 2),
    ]
    assert not analyze_subunit_interface(complex_snapshot, "2", "A", "B").contacts
    assert len(analyze_subunit_interface(complex_snapshot, "1", "B", "A", 3).contacts) == 1


def test_atom_chemistry_excludes_backbone_salt_and_backbone_hydrophobic(tmp_path):
    path = tmp_path / "chem.pdb"
    path.write_text(
        pdb_atom(1, "O", "LYS", "A", 1, 0, "O")
        + pdb_atom(2, "OD1", "ASP", "B", 1, 2, "O")
        + pdb_atom(3, "CA", "LEU", "A", 2, 15, "C")
        + pdb_atom(4, "CB", "LEU", "A", 2, 20, "C")
        + pdb_atom(5, "CA", "VAL", "B", 2, 14, "C")
        + pdb_atom(6, "CB", "VAL", "B", 2, 10, "C")
        + "END\n"
    )
    contacts = analyze_subunit_interface(capture_snapshot(path), "0", "A", "B").contacts
    assert len(contacts) == 2
    assert all(not p.interactions for p in contacts)


@pytest.mark.parametrize(
    "model,a,b,cutoff",
    [("3", "A", "B", 5), ("1", "A", "A", 5), ("1", "A", "X", 5), ("1", "A", "B", float("nan")), ("1", "A", "B", 0)],
)
def test_invalid_scope_is_rejected(complex_snapshot, model, a, b, cutoff):
    with pytest.raises(ValueError):
        analyze_subunit_interface(complex_snapshot, model, a, b, cutoff)


@pytest.mark.parametrize(
    ("res_a", "atom_a", "element_a", "res_b", "atom_b", "element_b", "expected"),
    [
        ("ALA", "N", "N", "GLY", "N", "N", set()),
        ("PRO", "N", "N", "GLY", "O", "O", set()),
        ("ASN", "ND2", "N", "ALA", "N", "N", set()),
        ("GLN", "NE2", "N", "ALA", "O", "O", {"hbond_geometric"}),
        ("HIS", "NE2", "N", "ASP", "OD2", "O", {"hbond_geometric"}),
        ("LEU", "CD1", "C", "VAL", "CG1", "C", {"hydrophobic"}),
    ],
)
def test_geometry_requires_appropriate_atoms(tmp_path, res_a, atom_a, element_a, res_b, atom_b, element_b, expected):
    path = tmp_path / "pair.pdb"
    path.write_text(
        pdb_atom(1, atom_a, res_a, "A", 1, 0, element_a) + pdb_atom(2, atom_b, res_b, "B", 2, 3, element_b) + "END\n"
    )
    result = analyze_subunit_interface(capture_snapshot(path), "0", "A", "B")
    assert {e.interaction_type.value for e in result.contacts[0].interactions} == expected


def test_exclude_nucleic_acids_but_include_modified_amino_acids(tmp_path):
    path = tmp_path / "mixed.pdb"
    path.write_text(
        pdb_atom(1, "P", "DA", "A", 1, 0, "P")
        + pdb_atom(2, "CA", "MSE", "A", 2, 4, "C")
        + pdb_atom(3, "CA", "ALA", "B", 2, 3, "C")
        + "END\n"
    )
    result = analyze_subunit_interface(capture_snapshot(path), "0", "A", "B")
    assert result.excluded_a == 1 and result.observed_a == 1
    assert [p.residue_a.residue_name for p in result.contacts] == ["MSE"]


@pytest.mark.parametrize(
    "extra_atom",
    [
        "ATOM 12 C CA . VAL Y 8 44 ? 6.0 5.0 6.0 1.00 22.0 ? 104 VAL A CA 1",
        "ATOM 12 O O . ALA Y 8 44 ? 6.0 5.0 6.0 1.00 22.0 ? 100 ALA A O 1",
    ],
)
def test_ambiguous_mmcif_author_chain_is_rejected(tmp_path, extra_atom):
    from pathlib import Path

    fixture = Path(__file__).parents[2] / "fixtures" / "parsing" / "enriched.mmcif"
    text = fixture.read_text()
    extra = f"{extra_atom}\nATOM 13 C CA . VAL Z 9 45 ? 6.0 5.0 6.0 1.00 22.0 ? 105 VAL B CA 1\n"
    path = tmp_path / "mixed.cif"
    path.write_text(text.replace("ATOM 10 N", extra + "ATOM 10 N"))
    with pytest.raises(ValueError, match="ambiguous|multiple"):
        analyze_subunit_interface(capture_snapshot(path), "1", "A", "B")


def test_no_eligible_atoms_is_not_reported_as_empty_interface(tmp_path):
    path = tmp_path / "zero.pdb"
    path.write_text(pdb_atom(1, "N", "ALA", "A", 1, 0, "N", occ=0) + pdb_atom(2, "N", "ALA", "B", 1, 0, "N") + "END\n")
    with pytest.raises(ValueError, match="no observed"):
        analyze_subunit_interface(capture_snapshot(path), "0", "A", "B")


def test_unselected_mmcif_collision_does_not_reject_other_label_chains(tmp_path):
    from pathlib import Path

    fixture = Path(__file__).parents[2] / "fixtures" / "parsing" / "enriched.mmcif"
    text = fixture.read_text()
    extra = (
        "ATOM 12 C CA . VAL Y 8 44 ? 5.0 5.0 6.0 1.00 22.0 ? 104 VAL B CA 1\n"
        "ATOM 13 N N . GLY S 6 45 ? 6.0 5.0 6.0 1.00 22.0 ? 105 GLY C N 1\n"
        "ATOM 14 C CA . GLY Z 9 46 ? 6.0 5.0 6.0 1.00 22.0 ? 105 GLY C CA 1\n"
    )
    path = tmp_path / "collision.cif"
    path.write_text(text.replace("ATOM 10 N", extra + "ATOM 10 N"))
    result = analyze_subunit_interface(capture_snapshot(path), "1", "A", "B")
    assert result.contacts[0].distance_angstrom == 0
    assert result.contacts[0].atom_a == "CA"  # Higher-occupancy conformer, from X only.


def test_chooser_omits_duplicate_auth_chains_and_wrong_models(complex_snapshot):
    from dataclasses import replace

    from structlens.application.subunit_interface import eligible_interface_chains
    from structlens.core.parsing import InputSelection, load_structure_evidence

    selection = InputSelection(complex_snapshot.content_id, complex_snapshot.display_name, "pdb", "1")
    chains = load_structure_evidence(complex_snapshot, selection=selection).protein_structure.chains
    assert eligible_interface_chains(chains, "2") == ()
    assert eligible_interface_chains((*chains, replace(chains[0], label_chain_id="extra")), "1") == ("B",)


def test_pre_cancelled_interface_stops_without_parsing(complex_snapshot):
    from threading import Event

    from structlens.core.errors import AnalysisCancelledError

    stop = Event()
    stop.set()
    with pytest.raises(AnalysisCancelledError):
        analyze_subunit_interface(complex_snapshot, "1", "A", "B", cancel_event=stop)
