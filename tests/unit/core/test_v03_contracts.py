"""Focused tests for the immutable v0.3 scientific contracts."""

from dataclasses import FrozenInstanceError

import numpy as np
import pytest

from structlens.core.difference_maps import DistanceDifferenceMatrix, ResidueDisplacementVector
from structlens.core.evidence import (
    EvidenceCard,
    EvidenceQuality,
    InteractionEvidence,
    SequenceEvidence,
    SiteEvidence,
    StructureEvidence,
)
from structlens.core.interactions import (
    InteractionChange,
    InteractionDifference,
    InteractionRecord,
    InteractionType,
    ReferenceInteractionKey,
)
from structlens.core.models import ResidueId
from structlens.core.msa import AnalysisSequence, MSAColumn, MSAResidueCell, SequenceResidueRef
from structlens.core.sites import SiteDefinition, SiteDefinitionMode, SiteMetrics


def _residue(number: str = "10", *, structure_id: str = "ref") -> ResidueId:
    return ResidueId(structure_id, "1", "A", number, None, "ALA")


def test_sequence_contract_preserves_source_residue_mapping() -> None:
    residue_ref = SequenceResidueRef(sequence_index=0, one_letter="A", residue_id=_residue())
    sequence = AnalysisSequence(
        structure_id="ref",
        chain_id="A",
        sequence="A",
        residues=(residue_ref,),
        source="structure",
    )

    assert sequence.residues == (residue_ref,)
    assert sequence.residues[0].residue_id == _residue()
    with pytest.raises(FrozenInstanceError):
        residue_ref.sequence_index = 1  # type: ignore[misc]


@pytest.mark.parametrize("source", ["structure", "fasta", "derived"])
def test_analysis_sequence_accepts_only_declared_sources(source: str) -> None:
    sequence = AnalysisSequence("ref", "A", "A", (SequenceResidueRef(0, "A", _residue()),), source)

    assert sequence.source == source


def test_analysis_sequence_rejects_invalid_source_and_mapping() -> None:
    with pytest.raises(ValueError, match="source"):
        AnalysisSequence("ref", "A", "A", (SequenceResidueRef(0, "A", _residue()),), "database")
    with pytest.raises(ValueError, match="sequence"):
        AnalysisSequence("ref", "A", "AA", (SequenceResidueRef(0, "A", _residue()),), "structure")


def test_msa_contract_exposes_reference_relative_column_semantics() -> None:
    residue_ref = SequenceResidueRef(0, "A", _residue())
    cell = MSAResidueCell(structure_id="ref", alignment_column=3, residue=residue_ref, character="A")
    column = MSAColumn(
        index=3,
        reference_label="A10",
        reference_residue=_residue(),
        cells=(cell,),
        non_gap_count=1,
        gap_fraction=0.0,
        ambiguous_fraction=0.0,
        conservation_score=1.0,
        entropy_bits=0.0,
    )

    assert column.cells[0].residue == residue_ref
    assert column.reference_residue == _residue()
    assert column.alignment_index == 3


def test_msa_gap_cell_keeps_missing_residue_explicit() -> None:
    cell = MSAResidueCell("target", 4, None, "-")
    column = MSAColumn(4, "ins-4", None, (cell,), 0, 1.0, 0.0, None, None)

    assert cell.residue is None
    assert column.reference_residue is None
    with pytest.raises(ValueError, match="alignment_column"):
        MSAColumn(4, "ins-4", None, (MSAResidueCell("target", 5, None, "-"),), 0, 1.0, 0.0, None, None)


def test_msa_metrics_validate_counts_fractions_and_finite_values() -> None:
    cell = MSAResidueCell("ref", 0, SequenceResidueRef(0, "A", _residue()), "A")
    with pytest.raises(ValueError, match="non_gap_count"):
        MSAColumn(0, "A10", _residue(), (cell,), 2, 0.0, 0.0, 1.0, 0.0)
    with pytest.raises(ValueError, match="gap_fraction"):
        MSAColumn(0, "A10", _residue(), (cell,), 1, float("nan"), 0.0, 1.0, 0.0)


def test_interaction_contract_uses_exact_types_and_explicit_evidence_mode() -> None:
    record = InteractionRecord(
        structure_id="target",
        interaction_type=InteractionType.HBOND_GEOMETRIC,
        residue_a=_residue(structure_id="target"),
        residue_b=_residue("20", structure_id="target"),
        atom_a="N",
        atom_b="O",
        distance_angstrom=3.1,
        angle_degrees=145.0,
        ligand_or_metal_id=None,
        evidence_mode="heavy_atom_geometry",
    )
    key = ReferenceInteractionKey(InteractionType.HBOND_GEOMETRIC, "A:10", "A:20", None)
    difference = InteractionDifference(key, InteractionChange.GAINED, None, record)

    assert record.evidence_mode == "heavy_atom_geometry"
    assert difference.change is InteractionChange.GAINED
    assert key.reference_position_a == "A:10"


def test_interaction_contract_rejects_nonfinite_measurements_and_invalid_changes() -> None:
    with pytest.raises(ValueError, match="distance_angstrom"):
        InteractionRecord(
            "ref",
            InteractionType.SALT_BRIDGE,
            _residue(),
            _residue("20"),
            None,
            None,
            float("nan"),
            None,
            None,
            "distance_cutoff",
        )
    key = ReferenceInteractionKey(InteractionType.SALT_BRIDGE, "A:10", "A:20", None)
    with pytest.raises(ValueError, match="target_record"):
        InteractionDifference(key, InteractionChange.GAINED, None, None)


def test_site_definition_and_metrics_match_v03_schema() -> None:
    definition = SiteDefinition(
        site_id="active-site",
        name="Active site",
        mode=SiteDefinitionMode.KEY_RESIDUES,
        reference_residues=(_residue(),),
        center_residue=None,
        ligand_id=None,
        radius_angstrom=None,
    )
    metrics = SiteMetrics(
        site_id="active-site",
        structure_id="target",
        mapped_residue_count=1,
        coverage_fraction=1.0,
        global_frame_backbone_rmsd_angstrom=None,
        site_fitted_backbone_rmsd_angstrom=None,
        centroid_displacement_angstrom=0.8,
        radius_of_gyration_angstrom=2.1,
        atomic_envelope_volume_angstrom3=None,
        sasa_angstrom2=42.0,
        polar_residue_fraction=0.5,
        charged_residue_fraction=None,
    )

    assert definition.reference_residues == (_residue(),)
    assert metrics.global_frame_backbone_rmsd_angstrom is None
    assert metrics.atomic_envelope_volume_angstrom3 is None


def test_site_contract_requires_mode_specific_inputs_and_valid_fractions() -> None:
    with pytest.raises(ValueError, match="radius_angstrom"):
        SiteDefinition("near-ligand", "Near ligand", SiteDefinitionMode.LIGAND_RADIUS, (), None, "ATP", None)
    with pytest.raises(ValueError, match="coverage_fraction"):
        SiteMetrics("site", "target", 1, 1.1, None, None, None, None, None, None, None, None)


def test_distance_difference_matrix_freezes_all_arrays() -> None:
    reference = np.array([[0.0, 2.0], [2.0, 0.0]], dtype=np.float32)
    target = np.array([[0.0, 3.5], [3.5, 0.0]], dtype=np.float64)
    matrix = DistanceDifferenceMatrix(
        reference_positions=("A:10", "A:20"),
        reference_distances_angstrom=reference,
        target_distances_angstrom=target,
        delta_angstrom=target - reference,
        valid_mask=np.ones((2, 2), dtype=bool),
    )

    assert matrix.delta_angstrom.dtype == np.float64
    assert np.array_equal(matrix.delta_angstrom, [[0.0, 1.5], [1.5, 0.0]])
    assert not matrix.reference_distances_angstrom.flags.writeable
    assert not matrix.target_distances_angstrom.flags.writeable
    assert not matrix.delta_angstrom.flags.writeable
    assert not matrix.valid_mask.flags.writeable


def test_distance_difference_matrix_validates_shape_finiteness_and_delta() -> None:
    with pytest.raises(ValueError, match="square"):
        DistanceDifferenceMatrix(("A:10",), [[0.0, 1.0]], [[0.0, 1.0]], [[0.0, 0.0]], [[True, True]])
    with pytest.raises(ValueError, match="finite"):
        DistanceDifferenceMatrix(("A:10",), [[float("inf")]], [[0.0]], [[0.0]], [[True]])
    with pytest.raises(ValueError, match="target minus reference"):
        DistanceDifferenceMatrix(("A:10",), [[0.0]], [[1.0]], [[0.0]], [[True]])


def test_residue_displacement_vector_uses_reference_frame_coordinates() -> None:
    vector = ResidueDisplacementVector(
        reference_position="A:10",
        reference_residue=_residue(),
        target_residue=_residue(structure_id="target"),
        start_xyz=(1.0, 2.0, 3.0),
        end_xyz=(2.0, 4.0, 5.0),
        vector_xyz=(1.0, 2.0, 2.0),
        magnitude_angstrom=3.0,
    )

    assert vector.magnitude_angstrom == pytest.approx(3.0)
    with pytest.raises(ValueError, match="vector_xyz"):
        ResidueDisplacementVector("A:10", _residue(), _residue(), (0, 0, 0), (1, 1, 1), (2, 2, 2), 3.0)


def test_evidence_card_has_explicit_sections_quality_and_no_scores() -> None:
    card = EvidenceCard(
        reference_residue=_residue(),
        target_id="target",
        sequence=SequenceEvidence(),
        structure=StructureEvidence(),
        interactions=InteractionEvidence(),
        site=SiteEvidence(),
        quality=EvidenceQuality(),
    )

    assert card.reference_residue == _residue()
    assert card.target_id == "target"
    assert card.quality.overall_status == "unavailable"
    assert not hasattr(card, "impact_score")
    assert not hasattr(card, "damage_score")


def test_nested_evidence_collections_are_immutable_tuples() -> None:
    quality = EvidenceQuality(available_sections=["sequence"], warnings=["partial mapping"])
    interactions = InteractionEvidence(reference_interactions=[])

    assert quality.available_sections == ("sequence",)
    assert quality.warnings == ("partial mapping",)
    assert interactions.reference_interactions == ()
