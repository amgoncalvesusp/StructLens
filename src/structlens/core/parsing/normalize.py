"""Convert Biopython structures into immutable StructLens domain records."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

from Bio.Data.IUPACData import protein_letters_3to1

from structlens.core.models import (
    AtomRecord,
    ProteinChain,
    ProteinStructure,
    ResidueId,
    ResidueNumbering,
    ResidueRecord,
)

_WATER_NAMES = frozenset({"HOH", "WAT", "DOD"})
_ONE_LETTER_BY_THREE_LETTER = {
    name.upper(): letter for name, letter in protein_letters_3to1.items()
}
_NumberingMap = Mapping[tuple[str, str, str, str | None, str], str | None]


def load_structure(path: Path) -> ProteinStructure:
    """Load PDB or mmCIF input selected by filename, including PDB gzip input."""

    suffixes = tuple(suffix.lower() for suffix in path.suffixes)
    if ".pdb" in suffixes or ".ent" in suffixes:
        from .pdb import load_pdb

        return load_pdb(path)
    if ".cif" in suffixes or ".mmcif" in suffixes:
        from .mmcif import load_mmcif

        return load_mmcif(path)
    raise ValueError(f"unsupported structure format: {path.name}")


def normalize_biopython_structure(
    structure: Any,
    path: Path,
    label_numbering: _NumberingMap | None = None,
) -> ProteinStructure:
    """Normalize ordering, identifiers, numbering, and selected alternate atoms."""

    chains: list[ProteinChain] = []
    for model in structure:
        model_id = str(getattr(model, "serial_num", model.id))
        for source_chain in model:
            records = tuple(
                record
                for residue in source_chain
                if (
                    record := _normalize_residue(
                        residue,
                        str(structure.id),
                        model_id,
                        str(source_chain.id),
                        label_numbering,
                    )
                )
                is not None
            )
            if records:
                chains.append(
                    ProteinChain(
                        structure_id=str(structure.id),
                        model_id=model_id,
                        chain_id=str(source_chain.id),
                        residues=tuple(record.residue_id for record in records),
                        sequence="".join(
                            record.one_letter or "X" for record in records
                        ),
                        residue_records=records,
                        source_path=str(path),
                    )
                )
    return ProteinStructure(
        structure_id=str(structure.id), chains=tuple(chains), source_path=str(path)
    )


def _normalize_residue(
    residue: Any,
    structure_id: str,
    model_id: str,
    chain_id: str,
    label_numbering: _NumberingMap | None,
) -> ResidueRecord | None:
    residue_name = str(residue.get_resname()).upper()
    if residue_name in _WATER_NAMES:
        return None
    _, raw_auth_seq, raw_insertion_code = residue.id
    auth_seq_id = str(raw_auth_seq)
    insertion_code = _none_if_blank(str(raw_insertion_code))
    residue_id = ResidueId(
        structure_id,
        model_id,
        chain_id,
        auth_seq_id,
        insertion_code,
        residue_name,
    )
    label_seq_id = None
    if label_numbering is not None:
        label_seq_id = label_numbering.get(
            (model_id, chain_id, auth_seq_id, insertion_code, residue_name)
        )
    one_letter = _ONE_LETTER_BY_THREE_LETTER.get(residue_name)
    return ResidueRecord(
        residue_id=residue_id,
        numbering=ResidueNumbering(auth_seq_id, label_seq_id, insertion_code),
        residue_name=residue_name,
        one_letter=one_letter,
        atoms=tuple(_atom_record(atom) for atom in _selected_atoms(residue)),
        is_standard=one_letter is not None,
    )


def _selected_atoms(residue: Any) -> list[Any]:
    selected: list[Any] = []
    for atom in residue.child_dict.values():
        if atom.is_disordered() == 2:
            choices = atom.disordered_get_list()
            selected.append(max(choices, key=_alternate_location_sort_key))
        else:
            selected.append(atom)
    return selected


def _alternate_location_sort_key(atom: Any) -> tuple[float, int, str]:
    occupancy = atom.get_occupancy()
    altloc = str(atom.get_altloc()).strip()
    return (float(occupancy) if occupancy is not None else -1.0, altloc == "", altloc)


def _atom_record(atom: Any) -> AtomRecord:
    return AtomRecord(
        name=str(atom.get_name()).strip(),
        element=str(atom.element).strip().upper(),
        coordinate=tuple(float(value) for value in atom.get_coord()),
        altloc=_none_if_blank(str(atom.get_altloc())),
        occupancy=(
            float(atom.get_occupancy()) if atom.get_occupancy() is not None else None
        ),
    )


def _none_if_blank(value: str) -> str | None:
    cleaned = value.strip()
    return None if cleaned in {"", ".", "?"} else cleaned


__all__ = ["load_structure", "normalize_biopython_structure"]
