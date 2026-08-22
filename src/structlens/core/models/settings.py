"""Serializable settings shared by alignment and analysis workflows."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class AlignmentMode(str, Enum):
    """How residue correspondence should be established."""

    AUTO = "auto"
    SEQUENCE = "sequence"
    STRUCTURE = "structure"
    MANUAL = "manual"


@dataclass(frozen=True, slots=True)
class SequenceAlignmentSettings:
    """Explicit pairwise sequence-alignment parameters."""

    substitution_matrix: str = "BLOSUM62"
    gap_open: float = -10.0
    gap_extend: float = -0.5


@dataclass(frozen=True, slots=True)
class StructuralAlignmentSettings:
    """Parameters reserved for the external structural alignment adapter."""

    executable: str = "USalign"
    timeout_seconds: float = 120.0


@dataclass(frozen=True, slots=True)
class AnalysisSettings:
    """Analysis policy, including the explicit AUTO branch thresholds."""

    alignment_mode: AlignmentMode = AlignmentMode.AUTO
    minimum_sequence_identity: float = 0.30
    minimum_sequence_coverage: float = 0.70
    substitution_matrix: str = "BLOSUM62"
    gap_open: float = -10.0
    gap_extend: float = -0.5
    refined_rmsd: bool = False
    refinement_cutoff_angstrom: float = 2.0
    refinement_max_iterations: int = 10
    usalign_executable: str | None = None

    def __post_init__(self) -> None:
        for name, value in (
            ("minimum_sequence_identity", self.minimum_sequence_identity),
            ("minimum_sequence_coverage", self.minimum_sequence_coverage),
        ):
            if not 0.0 <= value <= 1.0:
                raise ValueError(f"{name} must be between 0 and 1")
        if not isinstance(self.alignment_mode, AlignmentMode):
            object.__setattr__(
                self, "alignment_mode", AlignmentMode(self.alignment_mode)
            )
        if self.refinement_cutoff_angstrom <= 0:
            raise ValueError("refinement_cutoff_angstrom must be positive")
        if self.refinement_max_iterations < 1:
            raise ValueError("refinement_max_iterations must be at least 1")


__all__ = [
    "AlignmentMode",
    "AnalysisSettings",
    "SequenceAlignmentSettings",
    "StructuralAlignmentSettings",
]
