"""Distance-difference and displacement-vector contracts."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from structlens.core.models import ResidueId


def _matrix(value: ArrayLike, name: str, size: int | None = None) -> np.ndarray:
    matrix = np.asarray(value, dtype=np.float64)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"{name} must be square")
    if size is not None and matrix.shape != (size, size):
        raise ValueError(f"{name} shape must match reference_positions")
    if not np.isfinite(matrix).all():
        raise ValueError(f"{name} must be finite")
    frozen = np.array(matrix, dtype=np.float64, copy=True)
    frozen.setflags(write=False)
    return frozen


@dataclass(frozen=True, slots=True)
class DistanceDifferenceMatrix:
    reference_positions: tuple[str, ...]
    reference_distances_angstrom: ArrayLike
    target_distances_angstrom: ArrayLike
    delta_angstrom: ArrayLike
    valid_mask: ArrayLike

    def __post_init__(self) -> None:
        positions = tuple(self.reference_positions)
        size = len(positions)
        reference = _matrix(self.reference_distances_angstrom, "reference_distances_angstrom", size)
        target = _matrix(self.target_distances_angstrom, "target_distances_angstrom", size)
        delta = _matrix(self.delta_angstrom, "delta_angstrom", size)
        mask = np.asarray(self.valid_mask, dtype=bool)
        if mask.shape != (size, size):
            raise ValueError("valid_mask shape must match reference_positions")
        expected = target - reference
        if not np.allclose(delta[mask], expected[mask], rtol=1e-10, atol=1e-10):
            raise ValueError("delta_angstrom must equal target minus reference")
        if not np.allclose(reference, reference.T) or not np.allclose(target, target.T) or not np.allclose(delta, delta.T):
            raise ValueError("distance matrices must be symmetric")
        mask = np.array(mask, dtype=bool, copy=True)
        mask.setflags(write=False)
        object.__setattr__(self, "reference_positions", positions)
        object.__setattr__(self, "reference_distances_angstrom", reference)
        object.__setattr__(self, "target_distances_angstrom", target)
        object.__setattr__(self, "delta_angstrom", delta)
        object.__setattr__(self, "valid_mask", mask)

    @property
    def matrix_angstrom(self) -> np.ndarray:
        return np.asarray(self.delta_angstrom, dtype=np.float64)


@dataclass(frozen=True, slots=True)
class ResidueDisplacementVector:
    reference_position: str
    reference_residue: ResidueId
    target_residue: ResidueId
    start_xyz: tuple[float, float, float]
    end_xyz: tuple[float, float, float]
    vector_xyz: tuple[float, float, float]
    magnitude_angstrom: float

    def __post_init__(self) -> None:
        if not self.reference_position:
            raise ValueError("reference_position must not be empty")
        for name in ("start_xyz", "end_xyz", "vector_xyz"):
            value = tuple(float(item) for item in getattr(self, name))
            if len(value) != 3 or not all(math.isfinite(item) for item in value):
                raise ValueError(f"{name} must contain three finite values")
            object.__setattr__(self, name, value)
        expected = tuple(end - start for start, end in zip(self.start_xyz, self.end_xyz, strict=True))
        if not np.allclose(self.vector_xyz, expected, rtol=1e-9, atol=1e-9):
            raise ValueError("vector_xyz must equal end_xyz minus start_xyz")
        magnitude = math.sqrt(sum(item * item for item in self.vector_xyz))
        if not math.isfinite(self.magnitude_angstrom) or self.magnitude_angstrom < 0.0:
            raise ValueError("magnitude_angstrom must be finite and non-negative")
        if not math.isclose(self.magnitude_angstrom, magnitude, rel_tol=1e-9, abs_tol=1e-9):
            raise ValueError("magnitude_angstrom must match vector_xyz")

    @property
    def vector(self) -> tuple[float, float, float]:
        return self.vector_xyz


DisplacementVector = ResidueDisplacementVector

__all__ = ["DistanceDifferenceMatrix", "DisplacementVector", "ResidueDisplacementVector"]
