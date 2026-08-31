"""PDB parser compatibility entry points."""

from __future__ import annotations

from pathlib import Path

from structlens.core.models import ProteinStructure

from .normalize import load_structure_legacy


def load_pdb(path: Path) -> ProteinStructure:
    """Load a PDB (optionally gzip-compressed) into normalized records."""

    return load_structure_legacy(path)


__all__ = ["load_pdb"]
