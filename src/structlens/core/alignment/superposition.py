"""Reproducible strict structural superposition."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from structlens.core.geometry.kabsch import apply_transform, kabsch
from structlens.core.geometry.rmsd import rmsd


@dataclass(frozen=True, slots=True)
class SuperpositionResult:
    """The rigid transform that maps target coordinates onto reference coordinates."""

    rotation: NDArray[np.float64]
    translation: NDArray[np.float64]
    strict_rmsd_angstrom: float
    atom_count: int
    residue_count: int


def superpose(
    reference_coordinates: ArrayLike,
    target_coordinates: ArrayLike,
    *,
    residue_count: int | None = None,
) -> SuperpositionResult:
    """Fit target coordinates to reference and report the unreduced strict RMSD."""

    reference = np.asarray(reference_coordinates, dtype=float)
    target = np.asarray(target_coordinates, dtype=float)
    rotation, translation = kabsch(reference, target)
    atom_count = int(reference.shape[0])
    if residue_count is None:
        residue_count = atom_count
    if residue_count < 0:
        raise ValueError("residue_count must be non-negative")
    fitted_target = apply_transform(target, rotation, translation)
    rotation = rotation.copy()
    translation = translation.copy()
    rotation.setflags(write=False)
    translation.setflags(write=False)
    return SuperpositionResult(
        rotation=rotation,
        translation=translation,
        strict_rmsd_angstrom=rmsd(reference, fitted_target),
        atom_count=atom_count,
        residue_count=residue_count,
    )


fit_superposition = superpose


__all__ = ["SuperpositionResult", "fit_superposition", "superpose"]
