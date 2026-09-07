"""Aggregate structural metrics over aligned Cα coordinates."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from structlens.core.alignment.superposition import SuperpositionResult, superpose


@dataclass(frozen=True, slots=True)
class StructuralMetrics:
    strict_ca_rmsd_angstrom: float
    mapped_residue_count: int
    rotation: np.ndarray[Any, Any]
    translation: np.ndarray[Any, Any]


def calculate_structural_metrics(
    reference_ca_coordinates: np.ndarray[Any, Any],
    target_ca_coordinates: np.ndarray[Any, Any],
) -> StructuralMetrics:
    result: SuperpositionResult = superpose(reference_ca_coordinates, target_ca_coordinates)
    return StructuralMetrics(
        result.strict_rmsd_angstrom,
        result.residue_count,
        result.rotation,
        result.translation,
    )


__all__ = ["StructuralMetrics", "calculate_structural_metrics"]
