from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_allclose

from structlens.core.geometry.kabsch import apply_transform, kabsch


def test_kabsch_recovers_known_rigid_transform() -> None:
    target = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]]
    )
    rotation = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    translation = np.array([4.0, -2.0, 1.5])
    reference = target @ rotation + translation

    fitted_rotation, fitted_translation = kabsch(reference, target)

    assert_allclose(
        apply_transform(target, fitted_rotation, fitted_translation), reference
    )
    assert_allclose(fitted_rotation, rotation, atol=1e-12)
    assert_allclose(fitted_translation, translation, atol=1e-12)


def test_kabsch_returns_a_proper_rotation_for_reflected_coordinates() -> None:
    target = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]]
    )
    reflected_reference = target * np.array([-1.0, 1.0, 1.0])

    rotation, _ = kabsch(reflected_reference, target)

    assert np.linalg.det(rotation) == pytest.approx(1.0)


def test_kabsch_rejects_mismatched_coordinate_sets() -> None:
    with pytest.raises(ValueError, match="same shape"):
        kabsch(np.zeros((3, 3)), np.zeros((2, 3)))
