"""Per-residue distances after a global structural superposition."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike


def ca_displacement(
    reference_ca_coordinate: ArrayLike, target_ca_coordinate: ArrayLike
) -> float:
    """Return the Euclidean Cα displacement in Å for one aligned residue.

    Inputs must already share a coordinate frame.  This deliberately does not
    call the value an RMSD: one atom pair only defines a distance.
    """

    reference = np.asarray(reference_ca_coordinate, dtype=float)
    target = np.asarray(target_ca_coordinate, dtype=float)
    if reference.shape != (3,) or target.shape != (3,):
        raise ValueError("Cα coordinates must each have shape (3,)")
    if not np.isfinite(reference).all() or not np.isfinite(target).all():
        raise ValueError("Cα coordinates must be finite")
    return float(np.linalg.norm(reference - target))


__all__ = ["ca_displacement"]
