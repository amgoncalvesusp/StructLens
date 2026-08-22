from __future__ import annotations

import numpy as np
import pytest

from structlens.core.alignment.superposition import superpose


def test_superposition_reports_strict_rmsd_without_internal_deformation() -> None:
    target = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]]
    )
    rotation = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    reference = target @ rotation + np.array([2.0, 3.0, -1.0])

    result = superpose(reference, target, residue_count=4)

    assert result.strict_rmsd_angstrom == pytest.approx(0.0, abs=1e-12)
    assert result.atom_count == 4
    assert result.residue_count == 4
