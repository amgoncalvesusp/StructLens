"""Vectorized internal-distance differences and displacement vectors."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np

from structlens.core.difference_maps import DistanceDifferenceMatrix, ResidueDisplacementVector
from structlens.core.models import ResidueId


def _distance_matrix(coordinates: np.ndarray) -> np.ndarray:
    delta = coordinates[:, None, :] - coordinates[None, :, :]
    result = np.sqrt(np.einsum("ijk,ijk->ij", delta, delta, dtype=np.float64))
    return np.asarray(result, dtype=np.float64)


def calculate_distance_difference(
    reference_positions: Sequence[str],
    reference_ca: Mapping[str, Sequence[float]],
    target_ca: Mapping[str, Sequence[float]],
) -> DistanceDifferenceMatrix:
    labels = tuple(reference_positions)
    valid = np.asarray([label in reference_ca and label in target_ca for label in labels], dtype=bool)
    reference = np.zeros((len(labels), len(labels)), dtype=np.float64)
    target = np.full_like(reference, np.nan)
    if valid.any():
        ref_coords = np.asarray([reference_ca[label] for label, ok in zip(labels, valid, strict=True) if ok], dtype=np.float64)
        tar_coords = np.asarray([target_ca[label] for label, ok in zip(labels, valid, strict=True) if ok], dtype=np.float64)
        indices = np.flatnonzero(valid)
        reference[np.ix_(indices, indices)] = _distance_matrix(ref_coords)
        target[np.ix_(indices, indices)] = _distance_matrix(tar_coords)
    valid_mask = valid[:, None] & valid[None, :]
    delta = np.zeros_like(reference)
    delta[valid_mask] = target[valid_mask] - reference[valid_mask]
    return DistanceDifferenceMatrix(labels, reference, target, delta, valid_mask)


def build_displacement_vectors(
    positions: Sequence[str],
    reference_residues: Mapping[str, ResidueId],
    target_residues: Mapping[str, ResidueId],
    reference_ca: Mapping[str, Sequence[float]],
    target_ca_in_reference_frame: Mapping[str, Sequence[float]],
    *,
    minimum_magnitude_angstrom: float = 0.5,
    maximum_vectors: int = 100,
) -> tuple[ResidueDisplacementVector, ...]:
    values: list[ResidueDisplacementVector] = []
    for position in positions:
        if position not in reference_ca or position not in target_ca_in_reference_frame or position not in reference_residues or position not in target_residues:
            continue
        raw_start = tuple(float(item) for item in reference_ca[position])
        raw_end = tuple(float(item) for item in target_ca_in_reference_frame[position])
        if len(raw_start) != 3 or len(raw_end) != 3:
            raise ValueError("C-alpha coordinates must contain three values")
        start = (raw_start[0], raw_start[1], raw_start[2])
        end = (raw_end[0], raw_end[1], raw_end[2])
        vector = (end[0] - start[0], end[1] - start[1], end[2] - start[2])
        magnitude = float(np.linalg.norm(vector))
        if magnitude >= minimum_magnitude_angstrom:
            values.append(ResidueDisplacementVector(position, reference_residues[position], target_residues[position], start, end, vector, magnitude))
    values.sort(key=lambda item: (-item.magnitude_angstrom, item.reference_position))
    return tuple(values[:maximum_vectors])


__all__ = ["build_displacement_vectors", "calculate_distance_difference"]
