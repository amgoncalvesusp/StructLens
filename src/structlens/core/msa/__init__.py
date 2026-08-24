"""Immutable multiple-sequence-alignment domain contracts."""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Literal, Protocol

from structlens.core.models import ResidueId

MSASource = Literal["structure", "fasta", "derived"]


@dataclass(frozen=True, slots=True)
class SequenceResidueRef:
    sequence_index: int
    one_letter: str
    residue_id: ResidueId | None = None

    def __post_init__(self) -> None:
        if self.sequence_index < 0:
            raise ValueError("sequence_index must be non-negative")
        if len(self.one_letter) != 1:
            raise ValueError("one_letter must contain one character")


@dataclass(frozen=True, slots=True, init=False)
class AnalysisSequence:
    structure_id: str
    chain_id: str | None
    sequence: str
    residues: tuple[SequenceResidueRef, ...]
    source: MSASource

    def __init__(
        self,
        structure_id: str,
        chain_id: str | None = None,
        sequence: str = "",
        residues: tuple[SequenceResidueRef, ...] = (),
        source: MSASource = "derived",
        **legacy: object,
    ) -> None:
        if "sequence_id" in legacy and not structure_id:
            structure_id = str(legacy["sequence_id"])
        if "is_reference" in legacy or "display_name" in legacy or "source_uri" in legacy:
            # These v0.2-era presentation fields are intentionally ignored;
            # source identity is represented by structure_id/chain_id/source.
            pass
        if legacy:
            unknown = sorted(set(legacy) - {"sequence_id", "is_reference", "display_name", "source_uri"})
            if unknown:
                raise TypeError(f"unexpected AnalysisSequence fields: {', '.join(unknown)}")
        if not structure_id:
            raise ValueError("structure_id must not be empty")
        if source not in {"structure", "fasta", "derived"}:
            raise ValueError("source must be structure, fasta, or derived")
        normalized = tuple(residues)
        if len(sequence) != len(normalized):
            raise ValueError("sequence length must match residues")
        if tuple(item.sequence_index for item in normalized) != tuple(range(len(normalized))):
            raise ValueError("residue sequence_index values must be contiguous and ordered")
        object.__setattr__(self, "structure_id", structure_id)
        object.__setattr__(self, "chain_id", chain_id)
        object.__setattr__(self, "sequence", sequence)
        object.__setattr__(self, "residues", normalized)
        object.__setattr__(self, "source", source)

    @property
    def sequence_id(self) -> str:
        return self.structure_id


@dataclass(frozen=True, slots=True, init=False)
class MSAResidueCell:
    structure_id: str
    alignment_column: int
    residue: SequenceResidueRef | None
    character: str

    def __init__(
        self,
        structure_id: str,
        alignment_column: int,
        residue: SequenceResidueRef | None,
        character: str,
        **legacy: object,
    ) -> None:
        if legacy:
            raise TypeError(f"unexpected MSAResidueCell fields: {', '.join(sorted(legacy))}")
        if not structure_id:
            raise ValueError("structure_id must not be empty")
        if alignment_column < 0:
            raise ValueError("alignment_column must be non-negative")
        if len(character) != 1:
            raise ValueError("character must contain one symbol")
        is_gap = character in {"-", "."}
        if is_gap and residue is not None:
            raise ValueError("gap cells must not contain a residue")
        if not is_gap and residue is None:
            raise ValueError("non-gap cells require a residue")
        object.__setattr__(self, "structure_id", structure_id)
        object.__setattr__(self, "alignment_column", alignment_column)
        object.__setattr__(self, "residue", residue)
        object.__setattr__(self, "character", character)

    @property
    def symbol(self) -> str:
        return self.character

    @property
    def residue_ref(self) -> SequenceResidueRef | None:
        return self.residue


@dataclass(frozen=True, slots=True, init=False)
class MSAColumn:
    index: int
    reference_label: str
    reference_residue: ResidueId | None
    cells: tuple[MSAResidueCell, ...]
    non_gap_count: int
    gap_fraction: float
    ambiguous_fraction: float
    conservation_score: float | None
    entropy_bits: float | None

    def __init__(
        self,
        index: int,
        reference_label: str,
        reference_residue: ResidueId | None,
        cells: tuple[MSAResidueCell, ...],
        non_gap_count: int,
        gap_fraction: float,
        ambiguous_fraction: float,
        conservation_score: float | None,
        entropy_bits: float | None,
        **legacy: object,
    ) -> None:
        if legacy:
            raise TypeError(f"unexpected MSAColumn fields: {', '.join(sorted(legacy))}")
        if index < 0:
            raise ValueError("index must be non-negative")
        normalized = tuple(cells)
        if any(cell.alignment_column != index for cell in normalized):
            raise ValueError("alignment_column must match column index")
        if non_gap_count != sum(cell.character not in {"-", "."} for cell in normalized):
            raise ValueError("non_gap_count does not match cells")
        for name, value in (("gap_fraction", gap_fraction), ("ambiguous_fraction", ambiguous_fraction), ("conservation_score", conservation_score)):
            if value is not None and (not math.isfinite(value) or not 0.0 <= value <= 1.0):
                raise ValueError(f"{name} must be finite and between 0 and 1")
        if entropy_bits is not None and (not math.isfinite(entropy_bits) or entropy_bits < 0.0):
            raise ValueError("entropy_bits must be finite and non-negative")
        object.__setattr__(self, "index", index)
        object.__setattr__(self, "reference_label", reference_label)
        object.__setattr__(self, "reference_residue", reference_residue)
        object.__setattr__(self, "cells", normalized)
        object.__setattr__(self, "non_gap_count", non_gap_count)
        object.__setattr__(self, "gap_fraction", gap_fraction)
        object.__setattr__(self, "ambiguous_fraction", ambiguous_fraction)
        object.__setattr__(self, "conservation_score", conservation_score)
        object.__setattr__(self, "entropy_bits", entropy_bits)

    @property
    def alignment_index(self) -> int:
        return self.index


@dataclass(frozen=True, slots=True)
class MSASettings:
    """Reproducible settings for a multiple-sequence alignment run."""

    algorithm: str = "muscle5"
    mode: str = "align"


@dataclass(frozen=True, slots=True)
class MultipleSequenceAlignment:
    """Rows and reference-aware columns produced by an alignment engine."""

    sequences: tuple[AnalysisSequence, ...]
    aligned_rows: tuple[tuple[str, str], ...]
    columns: tuple[MSAColumn, ...]
    reference_structure_id: str | None = None
    algorithm: str = "fallback"
    provenance: tuple[str, ...] = field(default_factory=tuple)


class MultipleSequenceAlignmentEngine(Protocol):
    def align(
        self,
        sequences: Sequence[AnalysisSequence],
        settings: MSASettings,
    ) -> MultipleSequenceAlignment:
        ...


__all__ = [
    "AnalysisSequence",
    "MSAColumn",
    "MSAResidueCell",
    "MSASettings",
    "MSASource",
    "MultipleSequenceAlignment",
    "MultipleSequenceAlignmentEngine",
    "SequenceResidueRef",
]
