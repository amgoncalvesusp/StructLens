"""Residue mutation descriptors produced from a correspondence map."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from enum import Enum

from .residue import ResidueId


class MutationKind(str, Enum):
    """Descriptive mutation categories; none imply functional effects."""

    CONSERVED = "conserved"
    SUBSTITUTION = "substitution"
    INSERTION = "insertion"
    DELETION = "deletion"
    NONSTANDARD = "nonstandard"


SEQUENCE_CHANGE_KINDS = (MutationKind.SUBSTITUTION, MutationKind.INSERTION, MutationKind.DELETION)


@dataclass(frozen=True, slots=True)
class MutationEvent:
    """A mutation descriptor attached to one aligned position."""

    alignment_index: int
    kind: MutationKind
    reference: ResidueId | None
    target: ResidueId | None
    reference_aa: str | None
    target_aa: str | None
    reference_label: str
    target_label: str
    canonical_notation: str
    blosum62_score: int | None
    grantham_distance: int | None
    physicochemical_class: str | None


def count_sequence_changes(events: Iterable[MutationEvent]) -> int:
    return sum(event.kind in SEQUENCE_CHANGE_KINDS for event in events)


__all__ = ["SEQUENCE_CHANGE_KINDS", "MutationEvent", "MutationKind", "count_sequence_changes"]
