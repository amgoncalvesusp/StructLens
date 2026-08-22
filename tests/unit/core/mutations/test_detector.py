"""Tests for residue-map based mutation descriptors."""

from structlens.core.models import (
    CorrespondenceStatus,
    MutationKind,
    ResidueCorrespondence,
    ResidueId,
)
from structlens.core.mutations.blosum import blosum62_score
from structlens.core.mutations.detector import detect_mutations
from structlens.core.mutations.grantham import grantham_distance
from structlens.core.mutations.physicochemistry import (
    amino_acid_category,
    classify_substitution,
)


def _residue(number: str, name: str) -> ResidueId:
    return ResidueId("structure", "1", "A", number, None, name)


def _correspondence(
    status: CorrespondenceStatus,
    reference: ResidueId | None,
    target: ResidueId | None,
    reference_aa: str | None,
    target_aa: str | None,
) -> ResidueCorrespondence:
    return ResidueCorrespondence(
        alignment_index=7,
        reference=reference,
        target=target,
        reference_one_letter=reference_aa,
        target_one_letter=target_aa,
        status=status,
    )


def test_detects_conserved_residue_without_mutation_descriptors() -> None:
    event = detect_mutations(
        [
            _correspondence(
                CorrespondenceStatus.CONSERVED,
                _residue("130", "SER"),
                _residue("128", "SER"),
                "S",
                "S",
            )
        ]
    )[0]

    assert event.kind is MutationKind.CONSERVED
    assert event.canonical_notation == "S130S"
    assert event.blosum62_score is None
    assert event.grantham_distance is None
    assert event.physicochemical_class is None


def test_detects_substitution_and_describes_it() -> None:
    event = detect_mutations(
        [
            _correspondence(
                CorrespondenceStatus.SUBSTITUTION,
                _residue("130", "SER"),
                _residue("128", "THR"),
                "S",
                "T",
            )
        ]
    )[0]

    assert event.kind is MutationKind.SUBSTITUTION
    assert event.canonical_notation == "S130T"
    assert event.blosum62_score == 1
    assert event.grantham_distance == 58
    assert event.physicochemical_class == "conservative"


def test_detects_unanchored_insertion_with_structured_notation() -> None:
    event = detect_mutations(
        [
            _correspondence(
                CorrespondenceStatus.INSERTION, None, _residue("146", "GLY"), None, "G"
            )
        ]
    )[0]

    assert event.kind is MutationKind.INSERTION
    assert event.reference is None
    assert event.canonical_notation == "insertion:G:146"
    assert event.blosum62_score is None


def test_detects_deletion_using_reference_residue_numbering() -> None:
    event = detect_mutations(
        [
            _correspondence(
                CorrespondenceStatus.DELETION, _residue("145", "GLY"), None, "G", None
            )
        ]
    )[0]

    assert event.kind is MutationKind.DELETION
    assert event.canonical_notation == "ΔGly145"
    assert event.target is None
    assert event.grantham_distance is None


def test_nonstandard_residue_has_no_substitution_scores() -> None:
    event = detect_mutations(
        [
            _correspondence(
                CorrespondenceStatus.SUBSTITUTION,
                _residue("145", "MSE"),
                _residue("143", "MET"),
                "X",
                "M",
            )
        ]
    )[0]

    assert event.kind is MutationKind.NONSTANDARD
    assert event.blosum62_score is None
    assert event.grantham_distance is None
    assert event.physicochemical_class is None


def test_blosum62_lookup_is_symmetric_and_uses_embedded_canonical_values() -> None:
    assert blosum62_score("W", "W") == 11
    assert blosum62_score("S", "T") == 1
    assert blosum62_score("T", "S") == 1
    assert blosum62_score("X", "S") is None


def test_grantham_lookup_is_symmetric_and_limited_to_canonical_substitutions() -> None:
    assert grantham_distance("S", "T") == 58
    assert grantham_distance("T", "S") == 58
    assert grantham_distance("W", "C") == 215
    assert grantham_distance("X", "S") is None


def test_physicochemical_categories_are_explicit_descriptors() -> None:
    assert amino_acid_category("S") == "polar"
    assert amino_acid_category("R") == "positively_charged"
    assert classify_substitution("S", "T") == "conservative"
    assert classify_substitution("D", "K") == "non_conservative"
    assert classify_substitution("X", "S") is None
