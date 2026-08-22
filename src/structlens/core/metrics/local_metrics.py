"""Local neighborhood and residue-level structural metrics."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from structlens.core.geometry.rmsd import rmsd


def neighborhood_indices(
    reference_ca_coordinates: np.ndarray,
    center_index: int,
    radius_angstrom: float = 5.0,
) -> tuple[int, ...]:
    coordinates = np.asarray(reference_ca_coordinates, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError("reference_ca_coordinates must have shape (n, 3)")
    if not 0 <= center_index < len(coordinates):
        raise IndexError("center_index is outside the coordinate array")
    if radius_angstrom <= 0:
        raise ValueError("radius_angstrom must be positive")
    distances = np.linalg.norm(coordinates - coordinates[center_index], axis=1)
    return tuple(int(index) for index, distance in enumerate(distances) if distance <= radius_angstrom)


def local_rmsd(
    reference_ca_coordinates: np.ndarray,
    target_ca_coordinates: np.ndarray,
    indices: Sequence[int],
) -> float | None:
    reference = np.asarray(reference_ca_coordinates, dtype=float)
    target = np.asarray(target_ca_coordinates, dtype=float)
    selected = tuple(indices)
    if not selected:
        return None
    if reference.shape != target.shape:
        raise ValueError("reference and target coordinates must have the same shape")
    return rmsd(reference[list(selected)], target[list(selected)])


__all__ = ["local_rmsd", "neighborhood_indices"]
