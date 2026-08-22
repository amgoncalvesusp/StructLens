"""Stable contracts for PyMOL-independent alignment and mapping engines."""

from __future__ import annotations

from typing import Protocol

from structlens.core.models import (
    ProteinChain,
    ResidueCorrespondence,
    SequenceAlignmentSettings,
)

from .sequence import SequenceAlignmentResult


class SequenceAlignmentEngine(Protocol):
    """Produces a global amino-acid alignment for two normalized chains."""

    def align(
        self,
        reference: ProteinChain,
        target: ProteinChain,
        settings: SequenceAlignmentSettings,
    ) -> SequenceAlignmentResult: ...


class ResidueMapper(Protocol):
    """Converts an alignment into the authoritative correspondence map."""

    def build_correspondence(
        self,
        reference: ProteinChain,
        target: ProteinChain,
        alignment: SequenceAlignmentResult,
    ) -> list[ResidueCorrespondence]: ...


__all__ = ["ResidueMapper", "SequenceAlignmentEngine"]
