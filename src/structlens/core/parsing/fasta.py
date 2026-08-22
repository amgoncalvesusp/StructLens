"""FASTA parsing that preserves record boundaries and source descriptions."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO


@dataclass(frozen=True, slots=True)
class ProteinSequence:
    """One amino-acid sequence record from a FASTA file."""

    identifier: str
    description: str
    sequence: str


def load_fasta(path: Path) -> list[ProteinSequence]:
    """Read all FASTA records from ``path`` without joining separate entries."""

    return [
        ProteinSequence(
            identifier=record.id,
            description=record.description,
            sequence=str(record.seq).upper(),
        )
        for record in SeqIO.parse(path, "fasta")  # type: ignore[no-untyped-call]
    ]


__all__ = ["ProteinSequence", "load_fasta"]
