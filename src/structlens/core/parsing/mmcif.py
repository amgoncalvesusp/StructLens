"""mmCIF parser compatibility entry point."""

from __future__ import annotations

from pathlib import Path

from structlens.core.models import ProteinStructure

from .normalize import load_structure_legacy


def load_mmcif(path: Path) -> ProteinStructure:
    """Load an mmCIF (optionally gzip-compressed) into normalized records."""

    return load_structure_legacy(path)


__all__ = ["load_mmcif"]
