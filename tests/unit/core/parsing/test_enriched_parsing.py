"""Evidence-preserving PDB/mmCIF parsing tests."""

from __future__ import annotations

import gzip
import hashlib
from pathlib import Path

import pytest

from structlens.core.models import ComponentKind
from structlens.core.parsing import (
    AltlocPolicy,
    AssemblyScope,
    InputSelection,
    ParseLimits,
    StructureFormat,
    StructureParseError,
    load_mmcif,
    load_pdb,
    load_structure,
    load_structure_evidence,
    load_structure_legacy,
)
from structlens.core.parsing import normalize as normalize_module
from structlens.core.parsing.snapshot import SourceSnapshot

_FIXTURES = Path(__file__).parents[3] / "fixtures" / "parsing"


def _selection(path: Path, fmt: StructureFormat, model: str = "1") -> InputSelection:
    return InputSelection(
        content_id=hashlib.sha256(path.read_bytes()).hexdigest(),
        display_name=path.name,
        format=fmt,
        model_id=model,
    )


_MMCIF_COLUMNS = (
    "_atom_site.group_PDB",
    "_atom_site.id",
    "_atom_site.type_symbol",
    "_atom_site.label_atom_id",
    "_atom_site.label_alt_id",
    "_atom_site.label_comp_id",
    "_atom_site.label_asym_id",
    "_atom_site.label_entity_id",
    "_atom_site.label_seq_id",
    "_atom_site.pdbx_PDB_ins_code",
    "_atom_site.Cartn_x",
    "_atom_site.Cartn_y",
    "_atom_site.Cartn_z",
    "_atom_site.occupancy",
    "_atom_site.B_iso_or_equiv",
    "_atom_site.pdbx_formal_charge",
    "_atom_site.auth_seq_id",
    "_atom_site.auth_comp_id",
    "_atom_site.auth_asym_id",
    "_atom_site.auth_atom_id",
    "_atom_site.pdbx_PDB_model_num",
)


def _write_mmcif_rows(path: Path, rows: list[str]) -> None:
    path.write_text(
        "data_rows\n#\nloop_\n"
        + "\n".join(_MMCIF_COLUMNS)
        + "\n"
        + "\n".join(rows)
        + "\n#\n",
        encoding="ascii",
    )


def test_enriched_pdb_retains_components_altloc_and_metadata() -> None:
    parsed = load_structure_evidence(
        _FIXTURES / "enriched.pdb",
        selection=_selection(_FIXTURES / "enriched.pdb", StructureFormat.PDB),
    )

    assert parsed.metadata.available_models == ("1", "2")
    assert parsed.metadata.selected_model == "1"
    assert parsed.metadata.experimental_method == "X-RAY DIFFRACTION"
    assert parsed.metadata.resolution_angstrom == pytest.approx(1.8)
    assert parsed.metadata.label_chain_ids == ()
    assert parsed.protein_structure.chains[0].sequence == "AGX"
    assert parsed.protein_structure.chains[0].label_chain_id is None
    assert parsed.protein_structure.chains[0].residue_records[0].atoms[1].coordinate == pytest.approx(
        (12.0, 10.0, 10.0)
    )
    atom = parsed.protein_structure.chains[0].residue_records[0].atoms[1]
    assert atom.b_factor == pytest.approx(21.0)
    assert atom.occupancy == pytest.approx(0.7)
    assert atom.source_serial == 3
    assert atom.source_atom_id == "3"
    assert parsed.components[0].metadata["altloc_inventory"] == {"CA": ("A", "B")}
    assert [atom.altloc for atom in parsed.components[0].atoms if atom.name == "CA"] == ["A", "B"]
    assert {component.kind for component in parsed.components} == {
        ComponentKind.POLYMER_RESIDUE,
        ComponentKind.WATER,
        ComponentKind.ION,
        ComponentKind.LIGAND,
        ComponentKind.OTHER,
    }


def test_pdb_requires_explicit_model_and_supports_author_chain_selection() -> None:
    source = _FIXTURES / "enriched.pdb"
    with pytest.raises(ValueError, match="model_id"):
        load_structure_evidence(source)

    content_id = hashlib.sha256(source.read_bytes()).hexdigest()
    parsed = load_structure_evidence(
        source,
        selection=InputSelection(
            content_id=content_id,
            display_name=source.name,
            format=StructureFormat.PDB,
            model_id="2",
            author_chain_ids=("A",),
        ),
    )
    assert parsed.metadata.selected_model == "2"
    assert parsed.protein_structure.chains[0].author_chain_id == "A"


@pytest.mark.parametrize("loader", [load_structure, load_pdb])
def test_legacy_pdb_public_loaders_keep_all_models_and_hetero(loader) -> None:
    structure = loader(_FIXTURES / "enriched.pdb")

    assert [(chain.model_id, chain.sequence) for chain in structure.chains] == [
        ("1", "AGXXXX"),
        ("2", "A"),
    ]


def test_legacy_mmcif_public_loader_keeps_non_water_hetero() -> None:
    structure = load_mmcif(_FIXTURES / "enriched.mmcif")

    assert [(chain.model_id, chain.sequence) for chain in structure.chains] == [
        ("1", "AXXXX"),
        ("2", "A"),
    ]


def test_legacy_loader_keeps_ligand_only_chain(tmp_path: Path) -> None:
    path = tmp_path / "ligand_only.pdb"
    path.write_text(
        "HETATM    1  P   ATP B 202      18.000  10.000  10.000  1.00 42.00           P  \nEND\n",
        encoding="ascii",
    )

    structure = load_structure(path)

    assert len(structure.chains) == 1
    assert structure.chains[0].chain_id == "B"
    assert structure.chains[0].sequence == "X"
    assert structure.chains[0].residue_records[0].residue_name == "ATP"


def test_legacy_pdb_preflight_rejects_duplicate_serial(tmp_path: Path) -> None:
    path = tmp_path / "duplicate.pdb"
    atom = "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 10.00           N  "
    path.write_text(f"{atom}\n{atom}\nEND\n", encoding="ascii")

    with pytest.raises(StructureParseError, match="duplicate PDB atom serial"):
        load_structure(path)


def test_legacy_pdb_preflight_rejects_non_preservable_serial(tmp_path: Path) -> None:
    path = tmp_path / "overflow_serial.pdb"
    path.write_text(
        "ATOM  *****  N   ALA A   1       0.000   0.000   0.000  1.00 10.00           N  \nEND\n",
        encoding="ascii",
    )

    with pytest.raises(StructureParseError, match="not preservable"):
        load_structure(path)


def test_pdb_model_serial_is_canonicalized_for_source_metadata(tmp_path: Path) -> None:
    lines = (_FIXTURES / "enriched.pdb").read_text(encoding="ascii").splitlines()
    lines[3] = lines[3][:10] + "0001" + lines[3][14:]
    ion_index = next(index for index, line in enumerate(lines) if line.startswith("HETATM    7"))
    lines[ion_index] = lines[ion_index][:78] + "+1"
    path = tmp_path / "padded_model.pdb"
    path.write_text("\n".join(lines) + "\n", encoding="ascii")

    parsed = load_structure_evidence(path, selection=_selection(path, StructureFormat.PDB))

    assert parsed.metadata.selected_model == "1"
    ion = next(component for component in parsed.components if component.kind is ComponentKind.ION)
    assert ion.atoms[0].source_atom_id == "7"
    assert ion.atoms[0].formal_charge == 1


def test_single_explicit_pdb_model_uses_the_source_model_serial(tmp_path: Path) -> None:
    path = tmp_path / "single_padded_model.pdb"
    path.write_text(
        "MODEL     0001\n"
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 10.00           N  \n"
        "HETATM    7  NA  NA  A 201      17.000  10.000  10.000  1.00 41.00          NA1+\n"
        "ENDMDL\nEND\n",
        encoding="ascii",
    )

    parsed = load_structure_evidence(path)

    assert parsed.metadata.available_models == ("1",)
    assert parsed.metadata.selected_model == "1"
    ion = next(component for component in parsed.components if component.kind is ComponentKind.ION)
    assert ion.model_id == "1"
    assert ion.atoms[0].source_atom_id == "7"
    assert ion.atoms[0].formal_charge == 1


def test_legacy_blank_chain_uses_stripped_identity(tmp_path: Path) -> None:
    path = tmp_path / "blank_chain.pdb"
    path.write_text(
        "ATOM      1  N   ALA     1       0.000   0.000   0.000  1.00 10.00           N  \nEND\n",
        encoding="ascii",
    )

    structure = load_structure_legacy(path)

    assert structure.chains[0].chain_id == ""
    assert structure.chains[0].residues[0].chain_id == ""


def test_pdb_default_selection_supports_blank_and_named_chains(tmp_path: Path) -> None:
    path = tmp_path / "multi_chain.pdb"
    path.write_text(
        "ATOM      1  N   ALA     1       0.000   0.000   0.000  1.00 10.00           N  \n"
        "ATOM      2  N   GLY A   2       1.000   0.000   0.000  1.00 10.00           N  \n"
        "END\n",
        encoding="ascii",
    )

    parsed = load_structure_evidence(path)

    assert parsed.selection.author_chain_ids == ("", "A")
    assert {chain.chain_id for chain in parsed.protein_structure.chains} == {"", "A"}


def test_pdb_label_numbering_is_explicitly_missing() -> None:
    path = _FIXTURES / "enriched.pdb"
    parsed = load_structure_evidence(path, selection=_selection(path, StructureFormat.PDB))

    assert all(record.numbering.label_seq_id is None for record in parsed.protein_structure.chains[0].residue_records)


def test_malformed_pdb_is_reported_at_the_parser_boundary(tmp_path: Path) -> None:
    path = tmp_path / "malformed.pdb"
    path.write_text("ATOM      1  N   ALA A   1       not coordinates\nEND\n", encoding="ascii")

    with pytest.raises(StructureParseError, match="unable to parse PDB source"):
        load_structure(path)


def test_legacy_wrappers_keep_all_models_and_non_water_hetero_records() -> None:
    structure = load_structure_legacy(_FIXTURES / "enriched.pdb")

    assert [(chain.model_id, chain.sequence) for chain in structure.chains] == [("1", "AGXXXX"), ("2", "A")]
    assert [record.residue_name for record in structure.chains[0].residue_records] == [
        "ALA",
        "GLY",
        "MSE",
        "NA",
        "ATP",
        "UNX",
    ]


def test_mmcif_preserves_label_entity_and_formal_charge() -> None:
    parsed = load_structure_evidence(
        _FIXTURES / "enriched.mmcif",
        selection=_selection(_FIXTURES / "enriched.mmcif", StructureFormat.MMCIF),
    )

    chain = parsed.protein_structure.chains[0]
    assert chain.author_chain_id == "A"
    assert chain.label_chain_id == "X"
    assert chain.entity_id == "1"
    assert chain.residue_records[0].numbering.label_seq_id == "42"
    ion = next(component for component in parsed.components if component.kind is ComponentKind.ION)
    assert ion.atoms[0].formal_charge == 1
    assert ion.atoms[0].source_atom_id == "7"
    ligand = next(component for component in parsed.components if component.kind is ComponentKind.LIGAND)
    assert ligand.author_chain_id == "A"
    assert ligand.label_chain_id == "L"
    assert ligand.component_id != parsed.components[0].component_id
    assert parsed.protein_structure.chains[0].label_chain_id == "X"


def test_mmcif_numeric_source_ids_are_canonicalized_only_for_lookup(tmp_path: Path) -> None:
    path = tmp_path / "padded_ids.mmcif"
    _write_mmcif_rows(
        path,
        ["HETATM 0007 NA NA . NA X 2 . ? 0 0 0 1 10 1 1 NA A NA 0001"],
    )

    parsed = load_structure_evidence(path)

    ion = parsed.components[0]
    assert parsed.metadata.selected_model == "1"
    assert ion.atoms[0].source_serial == 7
    assert ion.atoms[0].source_atom_id == "0007"
    assert ion.atoms[0].formal_charge == 1


def test_mmcif_missing_label_chain_is_not_fabricated_from_author_chain(tmp_path: Path) -> None:
    path = tmp_path / "missing_label.mmcif"
    _write_mmcif_rows(
        path,
        ["ATOM 1 C C . ALA ? 1 1 ? 0 0 0 1 10 ? 1 ALA A C 1"],
    )

    parsed = load_structure_evidence(path)

    assert parsed.components[0].author_chain_id == "A"
    assert parsed.components[0].label_chain_id is None
    assert parsed.protein_structure.chains[0].label_chain_id is None
    assert parsed.metadata.label_chain_ids == ()


def test_mmcif_label_filter_is_applied_per_residue() -> None:
    path = _FIXTURES / "enriched.mmcif"
    selection = InputSelection(
        content_id=hashlib.sha256(path.read_bytes()).hexdigest(),
        display_name=path.name,
        format=StructureFormat.MMCIF,
        model_id="1",
        author_chain_ids=("A",),
        label_chain_ids=("X",),
    )
    parsed = load_structure_evidence(path, selection=selection)
    ligand = next(component for component in parsed.components if component.kind is ComponentKind.LIGAND)
    assert ligand.label_chain_id == "L"
    assert ligand.metadata["selected_for_analysis"] is False
    assert ligand.metadata["retention_scope"] == "selected_model"
    assert all(
        component.metadata["retention_scope"] == "selected_model" for component in parsed.components
    )
    assert all(
        component.metadata["selected_for_analysis"] is True
        for component in parsed.components
        if component.kind is ComponentKind.POLYMER_RESIDUE
    )


def test_mmcif_components_are_retained_for_unselected_chain_and_labels(tmp_path: Path) -> None:
    path = tmp_path / "chain_locators.mmcif"
    _write_mmcif_rows(
        path,
        [
            "ATOM 1 C C . ALA X 1 1 ? 0 0 0 1 10 ? 1 ALA A C 1",
            "HETATM 2 P P . ATP L 4 . ? 1 0 0 1 20 ? 2 ATP A P 1",
            "HETATM 3 C C1 . ATP Y 6 . ? 2 0 0 1 21 ? 3 ATP B C1 1",
        ],
    )
    selection = InputSelection(
        content_id=hashlib.sha256(path.read_bytes()).hexdigest(),
        display_name=path.name,
        format=StructureFormat.MMCIF,
        model_id="1",
        author_chain_ids=("A",),
        label_chain_ids=("X",),
    )

    parsed = load_structure_evidence(path, selection=selection)

    assert len(parsed.protein_structure.chains) == 1
    assert parsed.protein_structure.chains[0].label_chain_id == "X"
    assert {(component.author_chain_id, component.label_chain_id) for component in parsed.components} == {
        ("A", "X"),
        ("A", "L"),
        ("B", "Y"),
    }
    assert all(
        component.metadata["selected_for_analysis"] is False
        for component in parsed.components
        if component.kind is not ComponentKind.POLYMER_RESIDUE
    )


def test_mmcif_polymer_segments_same_author_are_grouped_and_metadata_deduplicated(tmp_path: Path) -> None:
    path = tmp_path / "segments.mmcif"
    _write_mmcif_rows(
        path,
        [
            "ATOM 1 C C . ALA X 1 1 ? 0 0 0 1 10 ? 1 ALA A C 1",
            "ATOM 2 C C . GLY Y 2 2 ? 1 0 0 1 11 ? 2 GLY A C 1",
        ],
    )

    parsed = load_structure_evidence(path)

    assert [(chain.author_chain_id, chain.label_chain_id) for chain in parsed.protein_structure.chains] == [
        ("A", "X"),
        ("A", "Y"),
    ]
    assert parsed.metadata.analyzed_chain_ids == ("A",)
    assert parsed.metadata.author_chain_ids == ("A",)


def test_glycan_without_connectivity_evidence_is_other(tmp_path: Path) -> None:
    path = tmp_path / "glycan.mmcif"
    _write_mmcif_rows(
        path,
        ["HETATM 1 C C1 . NAG X 2 . ? 0 0 0 1 10 ? 1 NAG A C1 1"],
    )

    parsed = load_structure_evidence(path)

    assert parsed.components[0].kind is ComponentKind.OTHER
    assert not parsed.protein_structure.chains


def test_mmcif_limit_fails_before_materializing_the_biotables(monkeypatch: pytest.MonkeyPatch) -> None:
    path = _FIXTURES / "enriched.mmcif"

    def fail_if_materialized(*args, **kwargs):
        raise AssertionError("MMCIF2Dict should not be called when preflight rejects")

    monkeypatch.setattr(normalize_module, "MMCIF2Dict", fail_if_materialized)
    with pytest.raises(ValueError, match="atoms"):
        load_structure_evidence(path, limits=ParseLimits(max_atoms=1))


def test_altloc_all_is_rejected_until_conformer_aware_analysis() -> None:
    path = _FIXTURES / "enriched.pdb"
    base = _selection(path, StructureFormat.PDB)
    selection = InputSelection(
        content_id=base.content_id,
        display_name=base.display_name,
        format=base.format,
        model_id=base.model_id,
        altloc_policy=AltlocPolicy.ALL,
    )
    with pytest.raises(ValueError, match="conformer-aware"):
        load_structure_evidence(path, selection=selection)


def test_snapshot_bytes_are_consumed_after_path_replacement(tmp_path: Path) -> None:
    original = (_FIXTURES / "enriched.pdb").read_bytes()
    path = tmp_path / "input.pdb"
    path.write_bytes(original)
    snapshot = SourceSnapshot.from_path(path)
    path.write_bytes(b"not a PDB")

    parsed = load_structure_evidence(
        snapshot,
        selection=InputSelection(
            content_id=snapshot.content_id,
            display_name=snapshot.display_name,
            format=StructureFormat.PDB,
            model_id="1",
        ),
    )
    assert parsed.raw_source_hash == snapshot.raw_sha256


def test_gzip_and_parse_limits_are_applied(tmp_path: Path) -> None:
    path = tmp_path / "input.pdb.gz"
    with gzip.open(path, "wb") as handle:
        handle.write((_FIXTURES / "enriched.pdb").read_bytes())
    parsed = load_structure_evidence(
        path,
        selection=InputSelection(
            content_id=hashlib.sha256((_FIXTURES / "enriched.pdb").read_bytes()).hexdigest(),
            display_name=path.name,
            format=StructureFormat.PDB,
            model_id="1",
        ),
    )
    assert parsed.metadata.format is StructureFormat.PDB
    with pytest.raises(ValueError, match="atoms"):
        load_structure_evidence(path, limits=ParseLimits(max_atoms=1))


@pytest.mark.parametrize(
    ("field", "message"),
    [("max_models", "models"), ("max_chains", "chains"), ("max_residues", "residues")],
)
def test_parse_limits_apply_to_entrypoint(field: str, message: str) -> None:
    path = _FIXTURES / "enriched.pdb"
    limits = ParseLimits(**{field: 1})
    with pytest.raises(ValueError, match=message):
        load_structure_evidence(path, limits=limits)


def test_biological_assembly_is_rejected_until_implemented() -> None:
    path = _FIXTURES / "enriched.pdb"
    content_id = hashlib.sha256(path.read_bytes()).hexdigest()
    with pytest.raises(ValueError, match="biological assembly"):
        load_structure_evidence(
            path,
            selection=InputSelection(
                content_id=content_id,
                display_name=path.name,
                format=StructureFormat.PDB,
                model_id="1",
                assembly_scope=AssemblyScope.BIOLOGICAL_ASSEMBLY,
            ),
        )


def test_mmcif_tokenized_preflight_rejects_duplicate_site_id(tmp_path: Path) -> None:
    path = tmp_path / "duplicate.mmcif"
    path.write_text(
        "data_duplicate\n#\nloop_\n"
        "_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n"
        "_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n"
        "_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n"
        "_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n"
        "_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n_atom_site.auth_seq_id\n"
        "_atom_site.auth_comp_id\n_atom_site.auth_asym_id\n_atom_site.auth_atom_id\n"
        "_atom_site.pdbx_PDB_model_num\n"
        "ATOM 1 C CA . ALA X 1 1 ? 0 0 0 1 1 1 ALA A CA 1\n"
        "ATOM 1 C CB . ALA X 1 1 ? 0 0 0 1 1 1 ALA A CB 1\n#\n",
        encoding="utf-8",
    )
    with pytest.raises(StructureParseError, match="duplicate mmCIF atom_site id"):
        load_structure_evidence(path)


def test_mmcif_tokenized_preflight_rejects_incomplete_loop(tmp_path: Path) -> None:
    path = tmp_path / "malformed.mmcif"
    path.write_text("data_bad\nloop_\n_atom_site.id\n_atom_site.group_PDB\n1 ATOM 2\n", encoding="utf-8")
    with pytest.raises(StructureParseError, match="incomplete"):
        load_structure_evidence(path)
