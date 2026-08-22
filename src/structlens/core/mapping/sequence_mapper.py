"""Build explicit residue correspondences from an aligned pair of sequences."""

from __future__ import annotations

from Bio.Align import substitution_matrices

from structlens.core.alignment.sequence import SequenceAlignmentResult
from structlens.core.models import (
    CorrespondenceStatus,
    ProteinChain,
    ResidueCorrespondence,
    ResidueId,
)

_CANONICAL_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
_BLOSUM62 = substitution_matrices.load("BLOSUM62")  # type: ignore[no-untyped-call]


class SequenceResidueMapper:
    """Translate each aligned column into one inspectable correspondence."""

    def build_correspondence(
        self,
        reference: ProteinChain,
        target: ProteinChain,
        alignment: SequenceAlignmentResult,
    ) -> list[ResidueCorrespondence]:
        """Create map entries without relying on residue-number equality."""

        reference_ids = _residue_ids(reference)
        target_ids = _residue_ids(target)
        reference_index = 0
        target_index = 0
        correspondences: list[ResidueCorrespondence] = []
        for alignment_index, (reference_aa, target_aa) in enumerate(
            zip(alignment.aligned_reference, alignment.aligned_target, strict=True)
        ):
            reference_id = (
                reference_ids[reference_index] if reference_aa != "-" else None
            )
            target_id = target_ids[target_index] if target_aa != "-" else None
            if reference_aa != "-":
                reference_index += 1
            if target_aa != "-":
                target_index += 1
            correspondences.append(
                ResidueCorrespondence(
                    alignment_index=alignment_index,
                    reference=reference_id,
                    target=target_id,
                    reference_one_letter=_letter_or_none(reference_aa),
                    target_one_letter=_letter_or_none(target_aa),
                    status=_status(reference_aa, target_aa),
                    sequence_score=_substitution_score(reference_aa, target_aa),
                    mapping_source="sequence",
                )
            )
        if reference_index != len(reference_ids) or target_index != len(target_ids):
            raise ValueError("alignment and chain residue counts are inconsistent")
        return correspondences


def _residue_ids(chain: ProteinChain) -> tuple[ResidueId, ...]:
    """Read Task 1's optional rich records while preserving old chain API."""

    residue_records = getattr(chain, "residue_records", ())
    if residue_records:
        return tuple(record.residue_id for record in residue_records)
    return chain.residues


def _letter_or_none(letter: str) -> str | None:
    return None if letter == "-" else letter


def _status(reference_aa: str, target_aa: str) -> CorrespondenceStatus:
    if reference_aa == "-":
        return CorrespondenceStatus.INSERTION
    if target_aa == "-":
        return CorrespondenceStatus.DELETION
    if (
        reference_aa not in _CANONICAL_AMINO_ACIDS
        or target_aa not in _CANONICAL_AMINO_ACIDS
    ):
        return CorrespondenceStatus.NONSTANDARD
    if reference_aa == target_aa:
        return CorrespondenceStatus.CONSERVED
    return CorrespondenceStatus.SUBSTITUTION


def _substitution_score(reference_aa: str, target_aa: str) -> float | None:
    if reference_aa == "-" or target_aa == "-":
        return None
    try:
        return float(_BLOSUM62[reference_aa, target_aa])
    except IndexError:
        return None


__all__ = ["SequenceResidueMapper"]
