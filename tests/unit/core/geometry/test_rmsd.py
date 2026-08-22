from __future__ import annotations

import numpy as np
import pytest

from structlens.core.geometry.rmsd import rmsd


def test_rmsd_is_zero_for_identical_coordinates() -> None:
    coordinates = np.array([[0.0, 1.0, 2.0], [2.0, 3.0, 4.0]])

    assert rmsd(coordinates, coordinates.copy()) == pytest.approx(0.0)


def test_rmsd_rejects_mismatched_coordinate_sets() -> None:
    with pytest.raises(ValueError, match="same shape"):
        rmsd(np.zeros((2, 3)), np.zeros((3, 3)))
