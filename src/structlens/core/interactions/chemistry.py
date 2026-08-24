"""Explicit residue chemistry and geometry helpers used by the detector.

The detector deliberately relies on typed atom groups rather than matching
arbitrary atom-name substrings.  Aromatic and cationic geometries are kept in
this module so all consumers use the same scientific definitions.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from structlens.core.models import AtomRecord, ResidueRecord

DONOR_ATOMS = frozenset({"N", "NE", "NH1", "NH2", "ND1", "ND2", "NE1", "NE2", "NZ", "OG", "OG1", "OH"})
ACCEPTOR_ATOMS = frozenset({"O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH", "N", "ND1", "NE2"})
HYDROPHOBIC_ELEMENTS = frozenset({"C", "S"})
AROMATIC_RESIDUES = frozenset({"PHE", "TYR", "TRP", "HIS"})
CATIONIC_RESIDUES = frozenset({"ARG", "LYS", "HIS"})
CHARGED_RESIDUES = frozenset({"ARG", "LYS", "HIS", "ASP", "GLU"})
METAL_ELEMENTS = frozenset({"ZN", "FE", "MG", "MN", "CA", "CU", "CO", "NI"})

# Each tuple is one chemically connected aromatic ring.  Tryptophan has two
# rings, and both are valid interaction surfaces.  The atom names are exact
# PDB/mmCIF names; no substring matching is used.
AROMATIC_RING_ATOMS: dict[str, tuple[tuple[str, ...], ...]] = {
    "PHE": (("CG", "CD1", "CE1", "CZ", "CE2", "CD2"),),
    "TYR": (("CG", "CD1", "CE1", "CZ", "CE2", "CD2"),),
    "TRP": (
        ("CG", "CD1", "NE1", "CE2", "CD2"),
        ("CD2", "CE2", "CE3", "CZ3", "CH2", "CZ2"),
    ),
    "HIS": (("CG", "ND1", "CE1", "NE2", "CD2"),),
}

# Explicit cationic groups.  Histidine is included as a potentially
# protonated cation, consistent with the existing descriptive residue table;
# callers still receive geometry evidence rather than a protonation claim.
CATIONIC_GROUP_ATOMS: dict[str, tuple[tuple[str, ...], ...]] = {
    "ARG": (("NE", "CZ", "NH1", "NH2"),),
    "LYS": (("NZ",),),
    "HIS": (("ND1", "NE2"),),
}


@dataclass(frozen=True, slots=True)
class AromaticRingGeometry:
    """Centroid, normal and source atoms for one aromatic ring."""

    residue_name: str
    atom_names: tuple[str, ...]
    centroid: tuple[float, float, float]
    normal: tuple[float, float, float]


@dataclass(frozen=True, slots=True)
class CationicGroupGeometry:
    """Centroid and source atoms for one explicit cationic group."""

    residue_name: str
    atom_names: tuple[str, ...]
    centroid: tuple[float, float, float]


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


def aromatic_ring_atom_names(residue_name: str) -> tuple[tuple[str, ...], ...]:
    """Return exact atom-name definitions for the residue's aromatic rings."""

    return AROMATIC_RING_ATOMS.get(residue_name.upper().strip(), ())


def cationic_group_atom_names(residue_name: str) -> tuple[tuple[str, ...], ...]:
    """Return exact atom-name definitions for explicit cationic groups."""

    return CATIONIC_GROUP_ATOMS.get(residue_name.upper().strip(), ())


def aromatic_ring_geometries(residue: ResidueRecord) -> tuple[AromaticRingGeometry, ...]:
    """Build available aromatic ring geometries from complete atom groups.

    A ring with missing atoms is omitted, rather than approximated.  This
    preserves the distinction between unavailable geometry and a zero-valued
    measurement.
    """

    atoms = {atom.name.upper().strip(): atom for atom in residue.atoms}
    geometries: list[AromaticRingGeometry] = []
    for names in aromatic_ring_atom_names(residue.residue_name):
        if not all(name in atoms for name in names):
            continue
        coordinates = tuple(_coordinate(atoms[name]) for name in names)
        centroid = _centroid(coordinates)
        normal = _plane_normal(coordinates, centroid)
        if normal is None:
            continue
        geometries.append(AromaticRingGeometry(residue.residue_name.upper(), names, centroid, normal))
    return tuple(geometries)


def cationic_group_geometries(residue: ResidueRecord) -> tuple[CationicGroupGeometry, ...]:
    """Build available explicit cationic-group geometries.

    All atoms in a group are required.  In particular, a partially parsed
    guanidinium group must not be silently treated as a fully observed cation.
    """

    atoms = {atom.name.upper().strip(): atom for atom in residue.atoms}
    geometries: list[CationicGroupGeometry] = []
    for names in cationic_group_atom_names(residue.residue_name):
        if not all(name in atoms for name in names):
            continue
        coordinates = tuple(_coordinate(atoms[name]) for name in names)
        geometries.append(CationicGroupGeometry(residue.residue_name.upper(), names, _centroid(coordinates)))
    return tuple(geometries)


def _centroid(coordinates: tuple[tuple[float, float, float], ...]) -> tuple[float, float, float]:
    if not coordinates:
        raise ValueError("coordinates must not be empty")
    count = float(len(coordinates))
    return tuple(sum(point[index] for point in coordinates) / count for index in range(3))  # type: ignore[return-value]


def _coordinate(atom: AtomRecord) -> tuple[float, float, float]:
    values = tuple(float(value) for value in atom.coordinate)
    if len(values) != 3:
        raise ValueError("atom coordinate must contain exactly three values")
    return values[0], values[1], values[2]


def _plane_normal(
    coordinates: tuple[tuple[float, float, float], ...],
    centroid: tuple[float, float, float],
) -> tuple[float, float, float] | None:
    centered = tuple(tuple(point[index] - centroid[index] for index in range(3)) for point in coordinates)
    best: tuple[float, float, float] | None = None
    best_norm = 0.0
    for first_index, first in enumerate(centered):
        for second in centered[first_index + 1 :]:
            cross = (
                first[1] * second[2] - first[2] * second[1],
                first[2] * second[0] - first[0] * second[2],
                first[0] * second[1] - first[1] * second[0],
            )
            norm = math.sqrt(sum(value * value for value in cross))
            if norm > best_norm:
                best = cross
                best_norm = norm
    if best is None or best_norm <= 1e-12:
        return None
    return tuple(value / best_norm for value in best)  # type: ignore[return-value]


__all__ = [
    "ACCEPTOR_ATOMS",
    "AROMATIC_RESIDUES",
    "AROMATIC_RING_ATOMS",
    "AromaticRingGeometry",
    "CATIONIC_RESIDUES",
    "CATIONIC_GROUP_ATOMS",
    "CationicGroupGeometry",
    "CHARGED_RESIDUES",
    "DONOR_ATOMS",
    "HYDROPHOBIC_ELEMENTS",
    "METAL_ELEMENTS",
    "atom_is_acceptor",
    "atom_is_donor",
    "aromatic_ring_atom_names",
    "aromatic_ring_geometries",
    "cationic_group_atom_names",
    "cationic_group_geometries",
    "residue_is_aromatic",
    "residue_is_cationic",
    "residue_is_hydrophobic",
]
