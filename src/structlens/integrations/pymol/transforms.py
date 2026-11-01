"""PyMOL transform calls through an injected command proxy."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np


def apply_transform(command: object, object_name: str, rotation: Sequence[Sequence[float]] | np.ndarray, translation: Sequence[float]) -> None:
    matrix = np.asarray(rotation, dtype=float)
    offset = tuple(float(value) for value in translation)
    if matrix.shape != (3, 3) or len(offset) != 3:
        raise ValueError("rotation must be 3x3 and translation must contain 3 values")
    # StructLens applies row vectors; PyMOL accepts a row-major homogeneous matrix acting on columns.
    transform = np.eye(4)
    transform[:3, :3] = matrix.T
    transform[:3, 3] = offset
    transform_fn = getattr(command, "transform_object", None)
    if transform_fn is None:
        raise RuntimeError("PyMOL command proxy does not expose transform_object")
    transform_fn(object_name, transform.reshape(-1).tolist(), 1, 0, "", 1)


__all__ = ["apply_transform"]
