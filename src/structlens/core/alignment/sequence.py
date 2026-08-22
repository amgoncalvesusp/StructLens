"""Global pairwise amino-acid alignment using explicit, serialized settings."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

from Bio import Align
from Bio.Align import substitution_matrices

from structlens.core.models import ProteinChain, SequenceAlignmentSettings


@dataclass(frozen=True, slots=True)
class SequenceAlignmentResult:
    """A reproducible global alignment represented as equal-length strings."""

    aligned_reference: str
    aligned_target: str
    score: float

    def __post_init__(self) -> None:
        if len(self.aligned_reference) != len(self.aligned_target):
            raise ValueError("aligned sequences must have the same length")


class GlobalSequenceAlignmentEngine:
    """Align normalized chain sequences with BLOSUM-like score matrices."""

    def align(
        self,
        reference: ProteinChain,
        target: ProteinChain,
        settings: SequenceAlignmentSettings,
    ) -> SequenceAlignmentResult:
        """Return the highest-scoring global alignment for the supplied chains."""

        aligner: Any = Align.PairwiseAligner(mode="global")  # type: ignore[no-untyped-call]
        aligner.substitution_matrix = substitution_matrices.load(  # type: ignore[no-untyped-call]
            settings.substitution_matrix
        )
        aligner.open_gap_score = settings.gap_open
        aligner.extend_gap_score = settings.gap_extend
        alignment = aligner.align(reference.sequence, target.sequence)[0]
        aligned_reference, aligned_target = _gapped_sequences(
            reference.sequence, target.sequence, alignment.coordinates
        )
        return SequenceAlignmentResult(
            aligned_reference=aligned_reference,
            aligned_target=aligned_target,
            score=float(alignment.score),
        )


def _gapped_sequences(
    reference: str, target: str, coordinates: Sequence[Sequence[int]]
) -> tuple[str, str]:
    """Expand Biopython coordinate blocks to portable gapped alignment strings."""

    # Biopython exposes a 2-by-N ndarray. Keeping the conversion here avoids
    # persisting a Biopython Alignment object as scientific project state.
    coordinate_rows = tuple(tuple(int(value) for value in row) for row in coordinates)
    reference_blocks: list[str] = []
    target_blocks: list[str] = []
    for index in range(len(coordinate_rows[0]) - 1):
        reference_start, reference_end = (
            coordinate_rows[0][index],
            coordinate_rows[0][index + 1],
        )
        target_start, target_end = (
            coordinate_rows[1][index],
            coordinate_rows[1][index + 1],
        )
        reference_span = reference[reference_start:reference_end]
        target_span = target[target_start:target_end]
        if len(reference_span) == len(target_span):
            reference_blocks.append(reference_span)
            target_blocks.append(target_span)
        elif reference_span:
            reference_blocks.append(reference_span)
            target_blocks.append("-" * len(reference_span))
        else:
            reference_blocks.append("-" * len(target_span))
            target_blocks.append(target_span)
    return "".join(reference_blocks), "".join(target_blocks)


__all__ = ["GlobalSequenceAlignmentEngine", "SequenceAlignmentResult"]
