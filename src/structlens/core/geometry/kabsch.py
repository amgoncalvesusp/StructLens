"""Rigid least-squares fitting with the Kabsch algorithm."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

Coordinates = NDArray[np.float64]


def _coordinates(value: ArrayLike, *, name: str) -> Coordinates:
    coordinates = np.asarray(value, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1:] != (3,):
        raise ValueError(f"{name} must have shape (n, 3)")
    if coordinates.shape[0] == 0:
        raise ValueError(f"{name} must contain at least one coordinate")
    if not np.isfinite(coordinates).all():
        raise ValueError(f"{name} must contain only finite coordinates")
    return coordinates


def kabsch(
    reference_coordinates: ArrayLike, target_coordinates: ArrayLike
) -> tuple[Coordinates, Coordinates]:
    """Fit ``target`` onto ``reference`` and return rotation and translation.

    Coordinates are row vectors: fitted coordinates are computed as
    ``target_coordinates @ rotation + translation``.  Reflection is explicitly
    disallowed, so the returned rotation always has determinant +1.
    """

    reference = _coordinates(reference_coordinates, name="reference_coordinates")
    target = _coordinates(target_coordinates, name="target_coordinates")
    if reference.shape != target.shape:
        raise ValueError(
            "reference_coordinates and target_coordinates must have the same shape"
        )

    reference_centroid = reference.mean(axis=0)
    target_centroid = target.mean(axis=0)
    covariance = (target - target_centroid).T @ (reference - reference_centroid)
    left, _, right_transposed = np.linalg.svd(covariance)
    rotation = left @ right_transposed
    if np.linalg.det(rotation) < 0.0:
        left[:, -1] *= -1.0
        rotation = left @ right_transposed
    translation = reference_centroid - target_centroid @ rotation
    return rotation, translation


def apply_transform(
    coordinates: ArrayLike, rotation: ArrayLike, translation: ArrayLike
) -> Coordinates:
    """Apply a row-vector rigid transform returned by :func:`kabsch`."""

    points = _coordinates(coordinates, name="coordinates")
    matrix = np.asarray(rotation, dtype=float)
    offset = np.asarray(translation, dtype=float)
    if matrix.shape != (3, 3):
        raise ValueError("rotation must have shape (3, 3)")
    if offset.shape != (3,):
        raise ValueError("translation must have shape (3,)")
    if not np.isfinite(matrix).all() or not np.isfinite(offset).all():
        raise ValueError("rotation and translation must be finite")
    return np.asarray(points @ matrix + offset, dtype=np.float64)


__all__ = ["Coordinates", "apply_transform", "kabsch"]
