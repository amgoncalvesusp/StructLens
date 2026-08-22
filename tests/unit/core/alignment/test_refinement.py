from __future__ import annotations

import numpy as np

from structlens.core.alignment.refinement import refine_superposition


def test_refinement_excludes_extreme_outlier_and_preserves_indices() -> None:
    reference = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [20, 0, 0]], dtype=float)
    target = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [20, 0, 0]], dtype=float)
    target[-1] += np.array([0, 10, 0])

    result = refine_superposition(reference, target, alignment_indices=(10, 11, 12, 13), cutoff_angstrom=1.0)

    assert result.excluded_alignment_indices == (13,)
    assert result.refined_superposition.residue_count == 3
    assert result.cutoff_angstrom == 1.0
