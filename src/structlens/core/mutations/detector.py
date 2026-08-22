"""Build mutation descriptors from the authoritative correspondence map."""

from __future__ import annotations

from collections.abc import Sequence

from structlens.core.models import (
    CorrespondenceStatus,
    MutationEvent,
    MutationKind,
    ResidueCorrespondence,
    ResidueId,
)

from .blosum import CANONICAL_AMINO_ACIDS, blosum62_score
from .grantham import grantham_distance
from .physicochemistry import classify_substitution

_CANONICAL_THREE_LETTER_NAMES = frozenset(
    {
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
    }
)


def detect_mutations(
    correspondences: Sequence[ResidueCorrespondence],
) -> list[MutationEvent]:
    """Describe mapped positions without revising the correspondence map.

    ``ResidueCorrespondence.status`` remains authoritative for mapped states.
    An explicit non-standard status, a non-canonical residue name, or a
    non-canonical one-letter symbol produces a ``NONSTANDARD`` descriptor and
    intentionally leaves canonical substitution scores unavailable.
    """

    return [
        _event_from_correspondence(correspondence)
        for correspondence in correspondences
        if correspondence.status is not CorrespondenceStatus.UNMAPPED
    ]


def _event_from_correspondence(
    correspondence: ResidueCorrespondence,
) -> MutationEvent:
    reference_aa = _normalized_symbol(correspondence.reference_one_letter)
    target_aa = _normalized_symbol(correspondence.target_one_letter)
    kind = _mutation_kind(correspondence, reference_aa, target_aa)
    reference_label = _residue_label(correspondence.reference)
    target_label = _residue_label(correspondence.target)

    if kind is MutationKind.SUBSTITUTION:
        return MutationEvent(
            alignment_index=correspondence.alignment_index,
            kind=kind,
            reference=correspondence.reference,
            target=correspondence.target,
            reference_aa=reference_aa,
            target_aa=target_aa,
            reference_label=reference_label,
            target_label=target_label,
            canonical_notation=f"{reference_aa}{reference_label}{target_aa}",
            blosum62_score=blosum62_score(reference_aa, target_aa),
            grantham_distance=grantham_distance(reference_aa, target_aa),
            physicochemical_class=classify_substitution(reference_aa, target_aa),
        )

    return MutationEvent(
        alignment_index=correspondence.alignment_index,
        kind=kind,
        reference=correspondence.reference,
        target=correspondence.target,
        reference_aa=reference_aa,
        target_aa=target_aa,
        reference_label=reference_label,
        target_label=target_label,
        canonical_notation=_notation(kind, correspondence, reference_aa, target_aa),
        blosum62_score=None,
        grantham_distance=None,
        physicochemical_class=None,
    )


def _mutation_kind(
    correspondence: ResidueCorrespondence,
    reference_aa: str | None,
    target_aa: str | None,
) -> MutationKind:
    if correspondence.status is CorrespondenceStatus.NONSTANDARD:
        return MutationKind.NONSTANDARD
    if not _is_canonical_residue(correspondence.reference, reference_aa):
        if correspondence.reference is not None:
            return MutationKind.NONSTANDARD
    if not _is_canonical_residue(correspondence.target, target_aa):
        if correspondence.target is not None:
            return MutationKind.NONSTANDARD
    return MutationKind(correspondence.status.value)


def _is_canonical_residue(residue: ResidueId | None, symbol: str | None) -> bool:
    if residue is None:
        return True
    return (
        residue.residue_name.upper() in _CANONICAL_THREE_LETTER_NAMES
        and symbol in CANONICAL_AMINO_ACIDS
    )


def _notation(
    kind: MutationKind,
    correspondence: ResidueCorrespondence,
    reference_aa: str | None,
    target_aa: str | None,
) -> str:
    if kind is MutationKind.CONSERVED:
        return (
            f"{reference_aa or '?'}{_residue_label(correspondence.reference)}"
            f"{target_aa or '?'}"
        )
    if kind is MutationKind.INSERTION:
        return f"insertion:{target_aa or '?'}:{_residue_label(correspondence.target)}"
    if kind is MutationKind.DELETION:
        residue_name = (
            correspondence.reference.residue_name.title()
            if correspondence.reference
            else "?"
        )
        return f"Δ{residue_name}{_residue_label(correspondence.reference)}"
    return _nonstandard_notation(correspondence, reference_aa, target_aa)


def _nonstandard_notation(
    correspondence: ResidueCorrespondence,
    reference_aa: str | None,
    target_aa: str | None,
) -> str:
    reference_name = (
        correspondence.reference.residue_name
        if correspondence.reference
        else reference_aa
    )
    target_name = (
        correspondence.target.residue_name if correspondence.target else target_aa
    )
    return (
        f"nonstandard:{reference_name or '?'}:"
        f"{_residue_label(correspondence.reference)}"
        f"->{target_name or '?'}:{_residue_label(correspondence.target)}"
    )


def _residue_label(residue: ResidueId | None) -> str:
    if residue is None:
        return "-"
    return f"{residue.auth_seq_id}{residue.insertion_code or ''}"


def _normalized_symbol(amino_acid: str | None) -> str | None:
    return amino_acid.upper() if amino_acid is not None else None


__all__ = ["detect_mutations"]
