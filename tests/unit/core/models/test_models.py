"""Contract tests for StructLens' foundational domain models.

These tests intentionally cover the identity and enum contracts before any
parsing or analysis service is introduced.
"""

from dataclasses import FrozenInstanceError, is_dataclass

import pytest

from structlens.core.models.components import ComponentKind, StructureComponent
from structlens.core.models.correspondence import (
    CorrespondenceStatus,
    ResidueCorrespondence,
)
from structlens.core.models.mutation import MutationEvent, MutationKind
from structlens.core.models.residue import ResidueId, ResidueNumbering
from structlens.core.models.settings import AlignmentMode, AnalysisSettings
from structlens.core.models.structure import AtomRecord, ProteinChain, ProteinStructure, ResidueRecord


def test_residue_id_distinguishes_insertion_codes() -> None:
    a = ResidueId("x", "1", "A", "100", None, "GLY")
    b = ResidueId("x", "1", "A", "100", "A", "GLY")

    assert a != b


def test_residue_id_is_immutable_and_hashable() -> None:
    residue_id = ResidueId("x", "1", "A", "100", None, "GLY")

    assert is_dataclass(residue_id)
    assert hash(residue_id) == hash(ResidueId("x", "1", "A", "100", None, "GLY"))
    with pytest.raises(FrozenInstanceError):
        residue_id.chain_id = "B"  # type: ignore[misc]


def test_residue_numbering_keeps_author_and_label_numbers() -> None:
    numbering = ResidueNumbering("100", "42", "A")

    assert numbering.auth_seq_id == "100"
    assert numbering.label_seq_id == "42"
    assert numbering.insertion_code == "A"


def test_correspondence_status_uses_enum() -> None:
    correspondence = ResidueCorrespondence(
        alignment_index=0,
        reference=ResidueId("x", "1", "A", "100", None, "GLY"),
        target=ResidueId("y", "1", "A", "101", None, "ALA"),
        reference_one_letter="G",
        target_one_letter="A",
        status=CorrespondenceStatus.SUBSTITUTION,
    )

    assert isinstance(correspondence.status, CorrespondenceStatus)
    assert correspondence.status.value == "substitution"


def test_correspondence_defaults_are_explicit_and_unit_named() -> None:
    correspondence = ResidueCorrespondence(
        alignment_index=4,
        reference=None,
        target=None,
        reference_one_letter=None,
        target_one_letter=None,
        status=CorrespondenceStatus.UNMAPPED,
    )

    assert correspondence.mapping_source == "unknown"
    assert correspondence.mapping_locked is False
    assert correspondence.ca_displacement_angstrom is None
    assert correspondence.is_outlier is False


def test_alignment_mode_values_are_stable() -> None:
    assert AlignmentMode.AUTO.value == "auto"
    assert AlignmentMode.SEQUENCE.value == "sequence"
    assert AlignmentMode.STRUCTURE.value == "structure"
    assert AlignmentMode.MANUAL.value == "manual"


def test_analysis_settings_exposes_auto_thresholds() -> None:
    settings = AnalysisSettings()

    assert settings.alignment_mode is AlignmentMode.AUTO
    assert settings.minimum_sequence_identity == pytest.approx(0.30)
    assert settings.minimum_sequence_coverage == pytest.approx(0.70)


def test_mutation_event_uses_mutation_kind_enum() -> None:
    event = MutationEvent(
        alignment_index=2,
        kind=MutationKind.SUBSTITUTION,
        reference=None,
        target=None,
        reference_aa="S",
        target_aa="T",
        reference_label="130",
        target_label="130",
        canonical_notation="S130T",
        blosum62_score=1,
        grantham_distance=58,
        physicochemical_class="conservative",
    )

    assert isinstance(event.kind, MutationKind)
    assert event.kind.value == "substitution"
    assert event.canonical_notation == "S130T"


def test_atom_record_normalizes_array_like_coordinate_to_immutable_tuple() -> None:
    atom = AtomRecord("CA", "C", [1, 2.5, 3], altloc="A", occupancy=0.5)

    assert atom.coordinate == (1.0, 2.5, 3.0)
    with pytest.raises(FrozenInstanceError):
        atom.name = "CB"  # type: ignore[misc]


def test_chain_keeps_legacy_residue_ids_and_rich_residue_records() -> None:
    residue_id = ResidueId("x", "1", "A", "100", None, "GLY")
    residue = ResidueRecord(
        residue_id=residue_id,
        numbering=ResidueNumbering("100", "100", None),
        residue_name="GLY",
        one_letter="G",
        atoms=(AtomRecord("CA", "C", (0, 0, 0)),),
    )
    chain = ProteinChain(
        "x",
        "1",
        "A",
        residues=(residue_id,),
        sequence="G",
        residue_records=(residue,),
    )

    assert chain.residues == (residue_id,)
    assert chain.residue_records[0].atoms[0].coordinate == (0.0, 0.0, 0.0)


def test_atom_record_retains_coordinate_evidence_and_source_identity() -> None:
    atom = AtomRecord(
        "CA",
        "C",
        (1, 2, 3),
        altloc="A",
        occupancy=0.75,
        b_factor=22.5,
        source_atom_id="atom-site-17",
        source_serial=17,
        formal_charge=0,
    )

    assert atom.b_factor == pytest.approx(22.5)
    assert atom.source_atom_id == "atom-site-17"
    assert atom.source_serial == 17
    assert atom.formal_charge == 0
    with pytest.raises(FrozenInstanceError):
        atom.b_factor = 1.0  # type: ignore[misc]


def test_protein_chain_retains_author_label_and_entity_ids() -> None:
    chain = ProteinChain(
        "structure",
        "1",
        "A",
        author_chain_id="auth-A",
        label_chain_id="label-A",
        entity_id="2",
    )

    assert chain.author_chain_id == "auth-A"
    assert chain.label_chain_id == "label-A"
    assert chain.entity_id == "2"


def test_structure_component_keeps_other_components_explicitly() -> None:
    atom = AtomRecord("C1", "C", (0, 0, 0), source_serial=8)
    metadata = {"residue_name": "UNX", "auth_chain_id": "A"}
    component = StructureComponent(
        component_id="1:A:8:UNX",
        kind=ComponentKind.OTHER,
        atoms=(atom,),
        metadata=metadata,
    )
    metadata["residue_name"] = "LIG"

    assert component.kind is ComponentKind.OTHER
    assert component.component_id == "1:A:8:UNX"
    assert component.atoms == (atom,)
    assert component.metadata["residue_name"] == "UNX"
    with pytest.raises(TypeError):
        component.metadata["new"] = "value"  # type: ignore[index]


def test_nested_structure_collections_are_tuples_and_validate_types() -> None:
    atom = AtomRecord("CA", "C", (0, 0, 0))
    residue_id = ResidueId("x", "1", "A", "1", None, "ALA")
    residue = ResidueRecord(
        residue_id=residue_id,
        numbering=ResidueNumbering("1", "1", None),
        residue_name="ALA",
        one_letter="A",
        atoms=[atom],  # type: ignore[arg-type]
    )
    assert residue.atoms == (atom,)
    with pytest.raises(TypeError):
        ResidueRecord(
            residue_id=residue_id,
            numbering=ResidueNumbering("1", "1", None),
            residue_name="ALA",
            one_letter="A",
            atoms=(object(),),  # type: ignore[arg-type]
        )

    chain = ProteinChain("x", "1", "A", residues=(residue_id,), residue_records=(residue,))
    structure = ProteinStructure("x", chains=[chain])  # type: ignore[arg-type]
    assert structure.chains == (chain,)
    with pytest.raises(TypeError):
        ProteinStructure("x", chains=(object(),))  # type: ignore[arg-type]


def test_nested_structure_metadata_is_deeply_immutable() -> None:
    metadata = {"source": {"tags": ["experimental"]}}
    chain = ProteinChain("x", "1", "A", metadata=metadata)
    structure = ProteinStructure("x", metadata=metadata)
    metadata["source"]["tags"].append("changed")

    assert chain.metadata["source"]["tags"] == ("experimental",)
    assert structure.metadata["source"]["tags"] == ("experimental",)
    with pytest.raises(TypeError):
        chain.metadata["source"]["tags"] += ("changed",)  # type: ignore[index]
    with pytest.raises(TypeError):
        structure.metadata["source"]["tags"] += ("changed",)  # type: ignore[index]


def test_source_serial_rejects_non_integer_non_string_values() -> None:
    with pytest.raises(ValueError):
        AtomRecord("CA", "C", (0, 0, 0), source_serial=1.5)  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        AtomRecord("CA", "C", (0, 0, 0), source_serial=b"1")  # type: ignore[arg-type]


def test_structure_models_validate_residue_identity_and_record_alignment() -> None:
    residue_id = ResidueId("structure", "1", "A", "10", "A", "ALA")
    residue = ResidueRecord(
        residue_id=residue_id,
        numbering=ResidueNumbering("10", "42", "A"),
        residue_name="ALA",
        one_letter="A",
    )
    chain = ProteinChain(
        "structure",
        "1",
        "A",
        residues=[residue_id],  # type: ignore[arg-type]
        residue_records=[residue],  # type: ignore[arg-type]
    )
    assert chain.residues == (residue_id,)
    assert chain.residue_records == (residue,)
    wrong_record = ResidueRecord(
        residue_id=ResidueId("structure", "1", "A", "11", None, "GLY"),
        numbering=ResidueNumbering("11", "43", None),
        residue_name="GLY",
        one_letter="G",
    )
    with pytest.raises(ValueError, match="residue IDs"):
        ProteinChain("structure", "1", "A", residues=(residue_id,), residue_records=(wrong_record,))
    with pytest.raises(ValueError, match="structure/model/chain"):
        ProteinChain(
            "other",
            "1",
            "A",
            residues=(residue_id,),
        )
    with pytest.raises(ValueError, match="residue_name"):
        ResidueRecord(
            residue_id=residue_id,
            numbering=ResidueNumbering("10", "42", "A"),
            residue_name="GLY",
            one_letter="G",
        )
    with pytest.raises(ValueError, match="numbering"):
        ResidueRecord(
            residue_id=residue_id,
            numbering=ResidueNumbering("11", "42", "A"),
            residue_name="ALA",
            one_letter="A",
        )


def test_pdb_blank_chain_identifiers_are_present_and_not_none() -> None:
    chain = ProteinChain("structure", "1", "", author_chain_id="", label_chain_id="")

    assert chain.chain_id == ""
    assert chain.author_chain_id == ""
    assert chain.label_chain_id == ""


def test_structure_metadata_rejects_non_json_like_values() -> None:
    class MutableValue:
        pass

    with pytest.raises(TypeError):
        ProteinChain("structure", "1", "A", metadata={"bad": MutableValue()})
    with pytest.raises(ValueError):
        ProteinChain("structure", "1", "A", metadata={"bad": float("nan")})
    with pytest.raises(TypeError):
        ProteinChain("structure", "1", "A", metadata={"bad": {"value"}})
