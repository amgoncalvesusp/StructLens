"""mmCIF reader entry point with author/label numbering retention."""

from __future__ import annotations

import gzip
import tempfile
from pathlib import Path

from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.MMCIFParser import MMCIFParser

from structlens.core.models import ProteinStructure

from .normalize import normalize_biopython_structure


def load_mmcif(path: Path) -> ProteinStructure:
    """Load an mmCIF file while retaining both auth and label sequence IDs."""

    parser = MMCIFParser(QUIET=True)  # type: ignore[no-untyped-call]
    if path.suffix.lower() != ".gz":
        structure = parser.get_structure(_structure_id(path), str(path))  # type: ignore[no-untyped-call]
        return normalize_biopython_structure(structure, path, _label_numbering(path))
    temporary_path: Path | None = None
    try:
        with (
            gzip.open(path, "rb") as source,
            tempfile.NamedTemporaryFile(suffix=".cif", delete=False) as destination,
        ):
            destination.write(source.read())
            temporary_path = Path(destination.name)
        structure = parser.get_structure(_structure_id(path), str(temporary_path))  # type: ignore[no-untyped-call]
        return normalize_biopython_structure(
            structure, path, _label_numbering(temporary_path)
        )
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def _label_numbering(
    path: Path,
) -> dict[tuple[str, str, str, str | None, str], str | None]:
    cif = MMCIF2Dict(str(path))  # type: ignore[no-untyped-call]
    values = zip(
        cif.get("_atom_site.pdbx_PDB_model_num", []),
        cif.get("_atom_site.auth_asym_id", []),
        cif.get("_atom_site.auth_seq_id", []),
        cif.get("_atom_site.pdbx_PDB_ins_code", []),
        cif.get("_atom_site.auth_comp_id", []),
        cif.get("_atom_site.label_seq_id", []),
        strict=True,
    )
    result: dict[tuple[str, str, str, str | None, str], str | None] = {}
    for model, chain, auth_seq, insertion, residue_name, label_seq in values:
        result[
            (
                str(model),
                str(chain),
                str(auth_seq),
                _none_if_unknown(str(insertion)),
                str(residue_name).upper(),
            )
        ] = _none_if_unknown(str(label_seq))
    return result


def _none_if_unknown(value: str) -> str | None:
    return None if value in {".", "?", ""} else value


def _structure_id(path: Path) -> str:
    name = path.name
    for suffix in (".gz", ".mmcif", ".cif"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
    return name


__all__ = ["load_mmcif"]
