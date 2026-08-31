from __future__ import annotations

import pytest

from structlens.core.pockets import (
    POCKET_RADII_VERSION,
    POCKET_VDW_RADII_ANGSTROM,
    vdw_radius_angstrom,
)


def test_pocket_radii_lookup_is_normalized_and_conservative() -> None:
    assert POCKET_RADII_VERSION == "structlens-pocket-vdw-bondi-polymer-se-1"
    assert vdw_radius_angstrom(" c ") == pytest.approx(1.70)
    assert vdw_radius_angstrom("N") == pytest.approx(1.55)
    assert vdw_radius_angstrom("SE") == pytest.approx(1.90)
    assert vdw_radius_angstrom(None) is None
    assert vdw_radius_angstrom("ZN") is None


def test_pocket_radii_table_is_immutable() -> None:
    with pytest.raises(TypeError):
        POCKET_VDW_RADII_ANGSTROM["C"] = 9.99  # type: ignore[index]
