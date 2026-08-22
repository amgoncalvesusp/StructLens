"""Authoritative residue-to-residue correspondence records."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from .residue import ResidueId


class CorrespondenceStatus(str, Enum):
    """The state of one aligned residue position."""

    CONSERVED = "conserved"
    SUBSTITUTION = "substitution"
    INSERTION = "insertion"
    DELETION = "deletion"
    NONSTANDARD = "nonstandard"
    UNMAPPED = "unmapped"


@dataclass(slots=True)
class ResidueCorrespondence:
    """One authoritative aligned position across reference and target.

    The record is intentionally not frozen: alignment and geometry stages
    progressively populate metrics while preserving the same correspondence
    identity. Identifiers themselves are immutable value objects.
    """

    alignment_index: int
    reference: ResidueId | None
    target: ResidueId | None
    reference_one_letter: str | None
    target_one_letter: str | None
    status: CorrespondenceStatus
    sequence_score: float | None = None
    ca_displacement_angstrom: float | None = None
    backbone_rmsd_angstrom: float | None = None
    sidechain_rmsd_angstrom: float | None = None
    all_heavy_atom_rmsd_angstrom: float | None = None
    is_outlier: bool = False
    is_key_residue: bool = False
    mapping_source: str = "unknown"
    mapping_locked: bool = False

    def __post_init__(self) -> None:
        """Normalize enum-compatible input and reject impossible indices."""

        if self.alignment_index < 0:
            raise ValueError("alignment_index must be non-negative")
        if not isinstance(self.status, CorrespondenceStatus):
            self.status = CorrespondenceStatus(self.status)


__all__ = ["CorrespondenceStatus", "ResidueCorrespondence"]
