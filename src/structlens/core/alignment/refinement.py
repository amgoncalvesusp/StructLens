"""Deterministic transparent refinement of a strict superposition."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from structlens.core.geometry.kabsch import apply_transform

from .superposition import SuperpositionResult, superpose


@dataclass(frozen=True, slots=True)
class RefinementResult:
    refined_superposition: SuperpositionResult
    included_alignment_indices: tuple[int, ...]
    excluded_alignment_indices: tuple[int, ...]
    cutoff_angstrom: float
    cycles: int


def refine_superposition(
    reference_coordinates: np.ndarray,
    target_coordinates: np.ndarray,
    *,
    alignment_indices: Sequence[int] | None = None,
    cutoff_angstrom: float = 2.0,
    max_iterations: int = 10,
) -> RefinementResult:
    if cutoff_angstrom <= 0:
        raise ValueError("cutoff_angstrom must be positive")
    if max_iterations < 1:
        raise ValueError("max_iterations must be at least 1")
    reference = np.asarray(reference_coordinates, dtype=float)
    target = np.asarray(target_coordinates, dtype=float)
    if reference.shape != target.shape or reference.ndim != 2 or reference.shape[1] != 3:
        raise ValueError("reference and target coordinates must both have shape (n, 3)")
    if alignment_indices is None:
        indices = tuple(range(len(reference)))
    else:
        indices = tuple(alignment_indices)
        if len(indices) != len(reference):
            raise ValueError("alignment_indices must match the coordinate count")
    keep = np.ones(len(reference), dtype=bool)
    cycles = 0
    result = superpose(reference, target, residue_count=len(reference))
    for cycle in range(1, max_iterations + 1):
        cycles = cycle
        result = superpose(reference[keep], target[keep], residue_count=int(keep.sum()))
        fitted = apply_transform(target, result.rotation, result.translation)
        distances = np.linalg.norm(reference - fitted, axis=1)
        new_keep = distances <= cutoff_angstrom
        if new_keep.sum() == 0 or np.array_equal(new_keep, keep):
            break
        keep = new_keep
    included = tuple(indices[i] for i, value in enumerate(keep) if value)
    excluded = tuple(indices[i] for i, value in enumerate(keep) if not value)
    return RefinementResult(result, included, excluded, cutoff_angstrom, cycles)


__all__ = ["RefinementResult", "refine_superposition"]
