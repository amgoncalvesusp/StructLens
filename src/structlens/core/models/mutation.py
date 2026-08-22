"""Residue mutation descriptors produced from a correspondence map."""

from __future__ import annotations

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


__all__ = ["MutationEvent", "MutationKind"]
