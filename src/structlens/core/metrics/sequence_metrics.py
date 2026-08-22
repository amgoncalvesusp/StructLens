"""Sequence identity, similarity, and coverage derived from correspondences."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from Bio.Align import substitution_matrices

from structlens.core.models import CorrespondenceStatus, ResidueCorrespondence

_BLOSUM62 = substitution_matrices.load("BLOSUM62")  # type: ignore[no-untyped-call]
_CANONICAL_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


@dataclass(frozen=True, slots=True)
class SequenceAlignmentMetrics:
    """Metrics over canonical residues, with explicit denominator semantics."""

    identity: float
    similarity: float
    coverage: float
    reference_coverage: float
    target_coverage: float
    aligned_canonical_residue_count: int


def calculate_sequence_metrics(
    correspondences: Sequence[ResidueCorrespondence],
) -> SequenceAlignmentMetrics:
    """Calculate mapped-alignment metrics without inferring from residue numbers."""

    reference_canonical = sum(
        _is_canonical(item.reference_one_letter) for item in correspondences
    )
    target_canonical = sum(
        _is_canonical(item.target_one_letter) for item in correspondences
    )
    aligned = [
        item
        for item in correspondences
        if _is_canonical(item.reference_one_letter)
        and _is_canonical(item.target_one_letter)
    ]
    aligned_count = len(aligned)
    identical = sum(item.status is CorrespondenceStatus.CONSERVED for item in aligned)
    similar = sum(
        _is_similar(item.reference_one_letter, item.target_one_letter)
        for item in aligned
    )
    reference_coverage = _ratio(aligned_count, reference_canonical)
    target_coverage = _ratio(aligned_count, target_canonical)
    return SequenceAlignmentMetrics(
        identity=_ratio(identical, aligned_count),
        similarity=_ratio(similar, aligned_count),
        coverage=reference_coverage,
        reference_coverage=reference_coverage,
        target_coverage=target_coverage,
        aligned_canonical_residue_count=aligned_count,
    )


def _is_canonical(letter: str | None) -> bool:
    return letter in _CANONICAL_AMINO_ACIDS


def _is_similar(reference: str | None, target: str | None) -> bool:
    if reference is None or target is None:
        return False
    try:
        return bool(_BLOSUM62[reference, target] > 0)
    except IndexError:
        return False


def _ratio(numerator: int, denominator: int) -> float:
    return numerator / denominator if denominator else 0.0


__all__ = ["SequenceAlignmentMetrics", "calculate_sequence_metrics"]
