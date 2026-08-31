"""Versioned van der Waals radii for pocket geometry.

The table is intentionally conservative and restricted to common biomolecular
polymer elements plus selenium. Metals remain unknown until a validated,
context-aware model exists.
"""

from __future__ import annotations

from collections.abc import Mapping
from types import MappingProxyType

POCKET_RADII_VERSION = "structlens-pocket-vdw-bondi-polymer-se-1"

POCKET_VDW_RADII_ANGSTROM: Mapping[str, float] = MappingProxyType(
    {
        "H": 1.20,
        "D": 1.20,
        "T": 1.20,
        "C": 1.70,
        "N": 1.55,
        "O": 1.52,
        "F": 1.47,
        "P": 1.80,
        "S": 1.80,
        "CL": 1.75,
        "BR": 1.85,
        "I": 1.98,
        "SE": 1.90,
    }
)


def vdw_radius_angstrom(element: str | None) -> float | None:
    if element is None:
        return None
    return POCKET_VDW_RADII_ANGSTROM.get(str(element).strip().upper())


__all__ = [
    "POCKET_RADII_VERSION",
    "POCKET_VDW_RADII_ANGSTROM",
    "vdw_radius_angstrom",
]
