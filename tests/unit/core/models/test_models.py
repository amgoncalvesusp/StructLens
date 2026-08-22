"""Contract tests for StructLens' foundational domain models.

These tests intentionally cover the identity and enum contracts before any
parsing or analysis service is introduced.
"""

from dataclasses import FrozenInstanceError, is_dataclass

import pytest

from structlens.core.models.correspondence import (
    CorrespondenceStatus,
    ResidueCorrespondence,
)
from structlens.core.models.mutation import MutationEvent, MutationKind
from structlens.core.models.residue import ResidueId, ResidueNumbering
from structlens.core.models.settings import AlignmentMode, AnalysisSettings
from structlens.core.models.structure import AtomRecord, ProteinChain, ResidueRecord


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
