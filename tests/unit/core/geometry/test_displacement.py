from __future__ import annotations

import numpy as np
import pytest

from structlens.core.geometry.displacement import ca_displacement


def test_ca_displacement_is_a_distance_not_a_single_atom_rmsd() -> None:
    result = ca_displacement(np.array([1.0, 2.0, 3.0]), np.array([1.0, 2.0, 5.5]))

    assert result == pytest.approx(2.5)
