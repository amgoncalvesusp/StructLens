"""Tests for the pre-parser raw coordinate quality boundary."""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from structlens.core.parsing import coordinate_audit as coordinate_audit_module
from structlens.core.parsing import normalize as normalize_module
from structlens.core.parsing.coordinate_audit import (
    CoordinateAudit,
    CoordinateQualityError,
    audit_mmcif_source,
    audit_pdb_source,
)
from structlens.core.parsing.limits import ParseLimits, SnapshotLimitError
from structlens.core.parsing.normalize import load_structure_evidence
from structlens.core.quality.models import CoordinateResidue


def _pdb_line(
    record: str = "ATOM",
    serial: int = 1,
    name: str = "CA",
    residue: str = "ALA",
    chain: str = "A",
    sequence: int = 1,
    x: str = "1.000",
    y: str = "2.000",
    z: str = "3.000",
    occupancy: str = "1.00",
    b_factor: str = "10.00",
    element: str = "C",
    altloc: str = " ",
) -> str:
    return (
        f"{record:<6}{serial:5d} {name:>4}{altloc}{residue:>3} {chain:1}"
        f"{sequence:4d}    {x:>8}{y:>8}{z:>8}{occupancy:>6}{b_factor:>6}"
        f"          {element:>2}  "
    )


def _mmcif(rows: list[dict[str, str]]) -> dict[str, list[str]]:
    fields = {
        "group_PDB": "ATOM",
        "id": "1",
        "type_symbol": "C",
        "label_atom_id": "CA",
        "auth_atom_id": "CA",
        "label_alt_id": ".",
        "label_comp_id": "ALA",
        "auth_comp_id": "ALA",
        "label_asym_id": "X",
        "auth_asym_id": "A",
        "label_seq_id": "1",
        "auth_seq_id": "1",
        "pdbx_PDB_ins_code": "?",
        "Cartn_x": "1.0",
        "Cartn_y": "2.0",
        "Cartn_z": "3.0",
        "occupancy": "1.0",
        "B_iso_or_equiv": "10.0",
        "pdbx_PDB_model_num": "1",
        "pdbx_formal_charge": "?",
    }
    columns: dict[str, list[str]] = {}
    for row in rows:
        current = fields | row
        for key, value in current.items():
            columns.setdefault(f"_atom_site.{key}", []).append(value)
    return columns


def test_valid_pdb_audit_preserves_every_atom_and_nonpolymer_backbone_is_not_required() -> None:
    source = (
        _pdb_line(name="N", element="N")
        + "\n"
        + _pdb_line(serial=2, name="CA")
        + "\n"
        + _pdb_line(record="HETATM", serial=3, name="ZN", residue="ZN", element="ZN")
        + "\nEND\n"
    ).encode("ascii")

    audit = audit_pdb_source(source, source_id="valid.pdb")

    assert isinstance(audit, CoordinateAudit)
    assert audit.can_normalize is True
    assert audit.report.error_count == 0
    assert sum(len(residue.atoms) for residue in audit.residues) == 3
    assert audit.residues[-1].is_polymer is False
    assert all(
        not (item.code == "backbone.missing_atoms" and item.residue_id and item.residue_id.endswith(":ZN"))
        for item in audit.report.diagnostics
    )


def test_missing_element_is_blocking_and_explicit_unknown_is_only_a_warning() -> None:
    source = (_pdb_line(element=" ") + "\n" + _pdb_line(serial=2, element="Q") + "\nEND\n").encode(
        "ascii"
    )

    audit = audit_pdb_source(source, raise_on_error=False)

    assert [item.code for item in audit.report.diagnostics].count("element.missing") == 1
    assert [item.code for item in audit.report.diagnostics].count("element.unknown_radius") == 1
    assert audit.residues[0].atoms[0].element == ""
    assert audit.residues[0].atoms[0].name == "CA"
    assert audit.can_normalize is False


def test_invalid_pdb_numbers_are_preserved_and_raise_the_same_typed_report() -> None:
    source = (
        _pdb_line(x="notnum", occupancy="2.00", b_factor="-1.00") + "\nEND\n"
    ).encode("ascii")

    audit = audit_pdb_source(source, raise_on_error=False)
    assert isinstance(audit.residues[0], CoordinateResidue)
    assert audit.residues[0].atoms[0].coordinate[0] == "notnum"
    assert audit.residues[0].atoms[0].occupancy == 2.0
    assert audit.residues[0].atoms[0].b_factor == -1.0
    assert {
        item.code for item in audit.report.diagnostics
    } >= {
        "coordinate.invalid",
        "occupancy.invalid",
        "b_factor.invalid",
    }
    with pytest.raises(CoordinateQualityError) as caught:
        audit_pdb_source(source)
    assert caught.value.report.error_count == 3


def test_pdb_altlocs_are_distinct_but_duplicate_semantic_identity_blocks() -> None:
    source = (
        _pdb_line(altloc="A")
        + "\n"
        + _pdb_line(serial=2, altloc="B")
        + "\n"
        + _pdb_line(serial=3, altloc="A")
        + "\nEND\n"
    ).encode("ascii")

    audit = audit_pdb_source(source, raise_on_error=False)

    assert [atom.altloc for atom in audit.residues[0].atoms] == ["A", "B", "A"]
    assert "atom.duplicate_identity" in [item.code for item in audit.report.diagnostics]


def test_pdb_atom_serial_may_repeat_in_distinct_models() -> None:
    source = (
        "MODEL        1\n"
        + _pdb_line(serial=1)
        + "\nENDMDL\nMODEL        2\n"
        + _pdb_line(serial=1)
        + "\nENDMDL\nEND\n"
    ).encode("ascii")

    audit = audit_pdb_source(source, raise_on_error=False)

    assert [residue.residue_id.model_id for residue in audit.residues] == ["1", "2"]
    assert "atom.duplicate_identity" not in {item.code for item in audit.report.diagnostics}


def test_direct_pdb_audit_enforces_atom_limit_before_materializing_residues(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = (_pdb_line(serial=1) + "\n" + _pdb_line(serial=2) + "\nEND\n").encode("ascii")

    def forbidden(*_: object) -> tuple[CoordinateResidue, ...]:
        raise AssertionError("raw residues must not be materialized past a failed preflight")

    monkeypatch.setattr(coordinate_audit_module, "_pdb_residues", forbidden)

    with pytest.raises(SnapshotLimitError, match="atoms"):
        audit_pdb_source(source, limits=ParseLimits(max_atoms=1))


def test_direct_pdb_path_enforces_byte_limit_before_reading(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    path = tmp_path / "too-large.pdb"
    path.write_text(_pdb_line() + "\nEND\n", encoding="ascii")

    def forbidden(_: Path) -> bytes:
        raise AssertionError("oversized path must not be read into memory")

    monkeypatch.setattr(Path, "read_bytes", forbidden)

    with pytest.raises(SnapshotLimitError, match="raw_bytes"):
        audit_pdb_source(path, limits=ParseLimits(max_raw_bytes=1))


def test_mmcif_audit_uses_type_symbol_without_fabricating_element_and_keeps_rows() -> None:
    audit = audit_mmcif_source(
        _mmcif(
            [
                {"id": "1", "type_symbol": "C"},
                {"id": "2", "type_symbol": "?", "label_atom_id": "N", "auth_atom_id": "N"},
            ]
        ),
        source_id="valid.cif",
        raise_on_error=False,
    )

    assert sum(len(residue.atoms) for residue in audit.residues) == 2
    assert audit.residues[0].atoms[1].element == ""
    assert any(item.code == "element.missing" for item in audit.report.diagnostics)


def test_mmcif_invalid_occupancy_b_and_coordinates_are_reported_without_repair() -> None:
    audit = audit_mmcif_source(
        _mmcif(
            [
                {
                    "Cartn_x": "NaN",
                    "occupancy": "-0.2",
                    "B_iso_or_equiv": "inf",
                }
            ]
        ),
        raise_on_error=False,
    )

    atom = audit.residues[0].atoms[0]
    assert math.isnan(float(atom.coordinate[0]))
    assert atom.occupancy == -0.2
    assert math.isinf(float(atom.b_factor))
    assert {
        item.code for item in audit.report.diagnostics
    } >= {
        "coordinate.non_finite",
        "occupancy.invalid",
        "b_factor.invalid",
    }


def test_mmcif_hetatm_does_not_receive_polymer_backbone_diagnostic() -> None:
    audit = audit_mmcif_source(
        _mmcif(
            [
                {
                    "group_PDB": "HETATM",
                    "id": "1",
                    "label_comp_id": "ATP",
                    "auth_comp_id": "ATP",
                    "label_atom_id": "P",
                    "auth_atom_id": "P",
                    "type_symbol": "P",
                }
            ]
        )
    )
    assert audit.residues[0].is_polymer is False
    assert all(item.code != "backbone.missing_atoms" for item in audit.report.diagnostics)


def test_pdb_loader_rejects_raw_invalidity_before_biopython_parser(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    path = tmp_path / "invalid-before-parser.pdb"
    path.write_text(_pdb_line(occupancy="2.00") + "\nEND\n", encoding="ascii")

    class ForbiddenParser:
        def __init__(self, **_: object) -> None:
            raise AssertionError("PDBParser must not be constructed before coordinate QC")

    monkeypatch.setattr(normalize_module, "PDBParser", ForbiddenParser)

    with pytest.raises(CoordinateQualityError) as caught:
        load_structure_evidence(path)

    assert "occupancy.invalid" in {item.code for item in caught.value.report.diagnostics}


def test_mmcif_loader_rejects_raw_invalidity_before_biopython_parser(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    path = tmp_path / "invalid-before-parser.cif"
    path.write_text(
        """data_invalid
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA . ALA A 1 ? 1.0 2.0 3.0 1.20 10.0 1 ALA A CA 1
#
""",
        encoding="ascii",
    )

    class ForbiddenParser:
        def __init__(self, **_: object) -> None:
            raise AssertionError("MMCIFParser must not be constructed before coordinate QC")

    monkeypatch.setattr(normalize_module, "MMCIFParser", ForbiddenParser)

    with pytest.raises(CoordinateQualityError) as caught:
        load_structure_evidence(path)

    assert "occupancy.invalid" in {item.code for item in caught.value.report.diagnostics}
