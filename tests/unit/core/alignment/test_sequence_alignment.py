"""Behavioural tests for global sequence alignment and explicit mapping."""

from __future__ import annotations

import pytest

from structlens.core.alignment.sequence import GlobalSequenceAlignmentEngine
from structlens.core.mapping.sequence_mapper import SequenceResidueMapper
from structlens.core.metrics.sequence_metrics import calculate_sequence_metrics
from structlens.core.models import (
    CorrespondenceStatus,
    ProteinChain,
    ResidueId,
    SequenceAlignmentSettings,
)


def make_chain(structure_id: str, sequence: str, start: int = 1) -> ProteinChain:
    residues = tuple(
        ResidueId(structure_id, "1", "A", str(start + index), None, "UNK")
        for index, _ in enumerate(sequence)
    )
    return ProteinChain(structure_id, "1", "A", residues, sequence)


def test_aligns_identical_sequences_and_reports_complete_identity() -> None:
    reference = make_chain("reference", "ACDE")
    target = make_chain("target", "ACDE")

    result = GlobalSequenceAlignmentEngine().align(
        reference, target, SequenceAlignmentSettings()
    )
    correspondences = SequenceResidueMapper().build_correspondence(
        reference, target, result
    )
    metrics = calculate_sequence_metrics(correspondences)

    assert result.aligned_reference == "ACDE"
    assert result.aligned_target == "ACDE"
    assert [item.status for item in correspondences] == [
        CorrespondenceStatus.CONSERVED,
    ] * 4
    assert metrics.identity == pytest.approx(1.0)
    assert metrics.coverage == pytest.approx(1.0)


def test_marks_canonical_amino_acid_difference_as_substitution() -> None:
    reference = make_chain("reference", "ACDE")
    target = make_chain("target", "ATDE")
    engine = GlobalSequenceAlignmentEngine()
    result = engine.align(reference, target, SequenceAlignmentSettings())

    correspondences = SequenceResidueMapper().build_correspondence(
        reference, target, result
    )

    assert correspondences[1].status is CorrespondenceStatus.SUBSTITUTION
    assert correspondences[1].reference_one_letter == "C"
    assert correspondences[1].target_one_letter == "T"
    assert correspondences[1].sequence_score is not None


def test_represents_target_insertion_and_reference_deletion_explicitly() -> None:
    reference = make_chain("reference", "ACDE")
    target = make_chain("target", "ACXDE")
    result = GlobalSequenceAlignmentEngine().align(
        reference, target, SequenceAlignmentSettings()
    )

    correspondences = SequenceResidueMapper().build_correspondence(
        reference, target, result
    )

    assert any(
        item.status is CorrespondenceStatus.INSERTION for item in correspondences
    )
    insertion = next(
        item
        for item in correspondences
        if item.status is CorrespondenceStatus.INSERTION
    )
    assert insertion.reference is None
    assert insertion.target is not None


def test_maps_sequences_with_different_residue_numbering_by_alignment_position() -> (
    None
):
    reference = make_chain("reference", "ACDE", start=70)
    target = make_chain("target", "ACDE", start=68)
    result = GlobalSequenceAlignmentEngine().align(
        reference, target, SequenceAlignmentSettings()
    )

    correspondences = SequenceResidueMapper().build_correspondence(
        reference, target, result
    )

    assert [
        pair.reference.auth_seq_id for pair in correspondences if pair.reference
    ] == [
        "70",
        "71",
        "72",
        "73",
    ]
    assert [pair.target.auth_seq_id for pair in correspondences if pair.target] == [
        "68",
        "69",
        "70",
        "71",
    ]
    assert all(
        pair.status is CorrespondenceStatus.CONSERVED for pair in correspondences
    )
