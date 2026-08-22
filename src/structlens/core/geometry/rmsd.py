"""Strict coordinate and symmetry-aware residue RMSD calculations."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .symmetry import symmetry_permutations

_BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O"})


def _coordinate_array(value: ArrayLike, *, name: str) -> NDArray[np.float64]:
    coordinates = np.asarray(value, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1:] != (3,):
        raise ValueError(f"{name} must have shape (n, 3)")
    if coordinates.shape[0] == 0:
        raise ValueError(f"{name} must contain at least one coordinate")
    if not np.isfinite(coordinates).all():
        raise ValueError(f"{name} must contain only finite coordinates")
    return coordinates


def rmsd(reference_coordinates: ArrayLike, target_coordinates: ArrayLike) -> float:
    """Return RMSD in Å for coordinates that are already in one frame."""

    reference = _coordinate_array(reference_coordinates, name="reference_coordinates")
    target = _coordinate_array(target_coordinates, name="target_coordinates")
    if reference.shape != target.shape:
        raise ValueError(
            "reference_coordinates and target_coordinates must have the same shape"
        )
    return float(np.sqrt(np.mean(np.sum((reference - target) ** 2, axis=1))))


def _atom_mapping(
    atoms: Mapping[str, ArrayLike] | Iterable[Any],
) -> dict[str, NDArray[np.float64]]:
    """Normalize mappings and future parser atom records without owning a model."""

    if isinstance(atoms, Mapping):
        items: list[tuple[Any, Any]] = list(atoms.items())
    else:
        items = []
        for atom in atoms:
            name = getattr(atom, "atom_name", getattr(atom, "name", None))
            coordinate = getattr(atom, "coordinates", getattr(atom, "coordinate", None))
            if name is None or coordinate is None:
                raise TypeError(
                    "atoms must be a name-to-coordinate mapping or atom records"
                )
            items.append((name, coordinate))
    normalized: dict[str, NDArray[np.float64]] = {}
    for name, coordinate in items:
        atom_name = str(name).strip().upper()
        point = np.asarray(coordinate, dtype=float)
        if not atom_name or point.shape != (3,) or not np.isfinite(point).all():
            raise ValueError(
                "atom names and coordinates must be finite three-dimensional values"
            )
        if atom_name in normalized:
            raise ValueError(f"duplicate atom name: {atom_name}")
        normalized[atom_name] = point
    return normalized


def _is_heavy_atom(name: str) -> bool:
    return not name.lstrip("0123456789").startswith("H")


def _symmetry_rmsd(
    reference: Mapping[str, NDArray[np.float64]],
    target: Mapping[str, NDArray[np.float64]],
    residue_name: str,
    names: frozenset[str],
) -> float:
    reference_coordinates = np.array([reference[name] for name in sorted(names)])
    return min(
        rmsd(
            reference_coordinates,
            np.array([target[permutation[name]] for name in sorted(names)]),
        )
        for permutation in symmetry_permutations(residue_name, names)
    )


@dataclass(frozen=True, slots=True)
class ResidueRmsdMetrics:
    """Independent per-residue metrics in Å; unavailable metrics are ``None``."""

    backbone_rmsd_angstrom: float | None
    sidechain_rmsd_angstrom: float | None
    all_heavy_atom_rmsd_angstrom: float | None


def residue_rmsds(
    reference_atoms: Mapping[str, ArrayLike] | Iterable[Any],
    target_atoms: Mapping[str, ArrayLike] | Iterable[Any],
    residue_name: str,
) -> ResidueRmsdMetrics:
    """Calculate backbone, side-chain and all-heavy-atom RMSDs in Å.

    A metric is unavailable when either residue lacks any atom in that metric's
    complete atom set.  The function therefore never drops a missing atom to
    fabricate a lower RMSD.  Symmetric names use only chemistry-valid swaps.
    """

    reference = _atom_mapping(reference_atoms)
    target = _atom_mapping(target_atoms)
    reference_heavy = frozenset(name for name in reference if _is_heavy_atom(name))
    target_heavy = frozenset(name for name in target if _is_heavy_atom(name))
    backbone = (
        _symmetry_rmsd(reference, target, residue_name, _BACKBONE_ATOMS)
        if _BACKBONE_ATOMS.issubset(reference) and _BACKBONE_ATOMS.issubset(target)
        else None
    )
    reference_sidechain = reference_heavy - _BACKBONE_ATOMS
    target_sidechain = target_heavy - _BACKBONE_ATOMS
    sidechain = (
        _symmetry_rmsd(reference, target, residue_name, reference_sidechain)
        if reference_sidechain and reference_sidechain == target_sidechain
        else None
    )
    all_heavy = (
        _symmetry_rmsd(reference, target, residue_name, reference_heavy)
        if reference_heavy and reference_heavy == target_heavy
        else None
    )
    return ResidueRmsdMetrics(backbone, sidechain, all_heavy)


__all__ = ["ResidueRmsdMetrics", "residue_rmsds", "rmsd"]
