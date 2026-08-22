"""PDB reader entry point."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from Bio.PDB.PDBParser import PDBParser

from structlens.core.models import ProteinStructure

from .normalize import normalize_biopython_structure


def load_pdb(path: Path) -> ProteinStructure:
    """Load a PDB (optionally gzip-compressed) into normalized domain records."""

    parser = PDBParser(QUIET=True)  # type: ignore[no-untyped-call]
    with _open_text(path) as handle:
        structure = parser.get_structure(_structure_id(path), handle)  # type: ignore[no-untyped-call]
    return normalize_biopython_structure(structure, path)


def _open_text(path: Path) -> Any:
    if path.suffix.lower() == ".gz":
        import gzip

        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def _structure_id(path: Path) -> str:
    name = path.name
    for suffix in (".gz", ".pdb", ".ent"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
    return name


__all__ = ["load_pdb"]
