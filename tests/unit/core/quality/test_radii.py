from __future__ import annotations

from types import MappingProxyType

from structlens.core.quality.clashes import (
    VAN_DER_WAALS_RADII_ANGSTROM as CLASH_ALIAS,
)
from structlens.core.quality.clashes import (
    VDW_RADII as CLASH_RADII,
)
from structlens.core.quality.clashes import (
    VDW_RADII_VERSION as CLASH_RADII_VERSION,
)
from structlens.core.quality.coordinates import vdw_radius_angstrom
from structlens.core.quality.radii import (
    VAN_DER_WAALS_RADII_ANGSTROM,
    VDW_RADII,
    VDW_RADII_VERSION,
)


def test_polymer_radius_table_is_one_immutable_versioned_source() -> None:
    assert isinstance(VDW_RADII, MappingProxyType)
    assert VDW_RADII is CLASH_RADII
    assert VDW_RADII is VAN_DER_WAALS_RADII_ANGSTROM
    assert VDW_RADII is CLASH_ALIAS
    assert VDW_RADII_VERSION == CLASH_RADII_VERSION
    assert VDW_RADII_VERSION.startswith("structlens-vdw-bondi-polymer-se-")


def test_table_contains_bondi_polymer_elements_and_explicit_selenium_extension() -> None:
    expected = {"H", "D", "T", "C", "N", "O", "F", "P", "S", "CL", "BR", "I", "SE"}
    assert set(VDW_RADII) == expected
    assert VDW_RADII["SE"] == 1.90
    assert vdw_radius_angstrom(" se ") == 1.90


def test_unsupported_or_metal_elements_are_unknown_instead_of_getting_ad_hoc_radii() -> None:
    assert vdw_radius_angstrom("Q") is None
    assert vdw_radius_angstrom("ZN") is None
    assert vdw_radius_angstrom(None) is None
