"""Explicit residue chemistry group tables used by the detector."""

from __future__ import annotations

from structlens.core.models import AtomRecord, ResidueRecord

DONOR_ATOMS = frozenset({"N", "NE", "NH1", "NH2", "ND1", "ND2", "NE1", "NE2", "NZ", "OG", "OG1", "OH"})
ACCEPTOR_ATOMS = frozenset({"O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH", "N", "ND1", "NE2"})
HYDROPHOBIC_ELEMENTS = frozenset({"C", "S"})
AROMATIC_RESIDUES = frozenset({"PHE", "TYR", "TRP", "HIS"})
CATIONIC_RESIDUES = frozenset({"ARG", "LYS", "HIS"})
CHARGED_RESIDUES = frozenset({"ARG", "LYS", "HIS", "ASP", "GLU"})
METAL_ELEMENTS = frozenset({"ZN", "FE", "MG", "MN", "CA", "CU", "CO", "NI"})


def atom_is_donor(atom: AtomRecord) -> bool:
    return atom.name.upper().strip() in DONOR_ATOMS and atom.element.upper() in {"N", "O", "S"}


def atom_is_acceptor(atom: AtomRecord) -> bool:
    return atom.name.upper().strip() in ACCEPTOR_ATOMS and atom.element.upper() in {"N", "O", "S"}


def residue_is_hydrophobic(residue: ResidueRecord) -> bool:
    backbone = {"N", "CA", "C", "O", "OXT"}
    return residue.residue_name.upper() in {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"} and any(
        atom.name.upper() not in backbone and atom.element.upper() in HYDROPHOBIC_ELEMENTS for atom in residue.atoms
    )


def residue_is_aromatic(residue: ResidueRecord) -> bool:
    return residue.residue_name.upper() in AROMATIC_RESIDUES


def residue_is_cationic(residue: ResidueRecord) -> bool:
    return residue.residue_name.upper() in CATIONIC_RESIDUES


__all__ = [
    "ACCEPTOR_ATOMS",
    "AROMATIC_RESIDUES",
    "CATIONIC_RESIDUES",
    "CHARGED_RESIDUES",
    "DONOR_ATOMS",
    "HYDROPHOBIC_ELEMENTS",
    "METAL_ELEMENTS",
    "atom_is_acceptor",
    "atom_is_donor",
    "residue_is_aromatic",
    "residue_is_cationic",
    "residue_is_hydrophobic",
]
