from __future__ import annotations

import numpy as np
import pytest

from structlens.core.geometry.rmsd import residue_rmsds, rmsd


def _asp_atoms() -> dict[str, np.ndarray]:
    return {
        "N": np.array([0.0, 0.0, 0.0]),
        "CA": np.array([1.0, 0.0, 0.0]),
        "C": np.array([2.0, 0.0, 0.0]),
        "O": np.array([3.0, 0.0, 0.0]),
        "CB": np.array([1.0, 1.0, 0.0]),
        "CG": np.array([1.0, 2.0, 0.0]),
        "OD1": np.array([0.0, 3.0, 0.0]),
        "OD2": np.array([2.0, 3.0, 0.0]),
    }


def test_symmetric_sidechain_names_produce_minimum_chemically_valid_rmsd() -> None:
    reference = _asp_atoms()
    target = _asp_atoms()
    target["OD1"], target["OD2"] = target["OD2"], target["OD1"]
    naive = rmsd(
        np.array([reference["OD1"], reference["OD2"]]),
        np.array([target["OD1"], target["OD2"]]),
    )

    metrics = residue_rmsds(reference, target, "ASP")

    assert naive > 0.0
    assert metrics.backbone_rmsd_angstrom == pytest.approx(0.0)
    assert metrics.sidechain_rmsd_angstrom == pytest.approx(0.0)
    assert metrics.all_heavy_atom_rmsd_angstrom == pytest.approx(0.0)


def test_missing_required_backbone_atom_produces_none_not_a_substitute() -> None:
    reference = _asp_atoms()
    target = _asp_atoms()
    del target["O"]

    metrics = residue_rmsds(reference, target, "ASP")

    assert metrics.backbone_rmsd_angstrom is None
    assert metrics.all_heavy_atom_rmsd_angstrom is None
