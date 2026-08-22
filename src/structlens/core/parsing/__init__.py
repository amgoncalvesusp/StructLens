"""Input parsing and normalization for structure and sequence files."""

from .fasta import ProteinSequence, load_fasta
from .normalize import load_structure

__all__ = ["ProteinSequence", "load_fasta", "load_structure"]
