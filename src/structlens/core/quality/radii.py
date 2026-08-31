"""Versioned van der Waals radii used by structural quality screens.

The table is intentionally limited to elements with a conservative Bondi-style
radius for biomolecular polymer screening, plus an explicit selenium extension
for selenocysteine and other selenium-containing residues.  Metals and other
elements are left unknown until a validated, context-aware model is added;
callers must report those atoms instead of silently assigning a radius.
"""

from __future__ import annotations

from collections.abc import Mapping
from types import MappingProxyType

# The version is part of every radius-dependent report.  A new table or
# interpretation must receive a new version so results cannot be compared
# silently across scientific conventions.
VDW_RADII_VERSION = "structlens-vdw-bondi-polymer-se-1"

# Values are in Å.  Bondi-style values cover the common biomolecular polymer
# elements; selenium is an explicit extension used by selenium-containing
# residues.  Do not add a metal without a documented, context-appropriate
# radius model.
VDW_RADII: Mapping[str, float] = MappingProxyType(
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

# Compatibility name retained for callers that use the unabbreviated form.
VAN_DER_WAALS_RADII_ANGSTROM = VDW_RADII


def vdw_radius_angstrom(element: str | None) -> float | None:
    """Return a supported radius in Å, or ``None`` for an unknown element."""

    if element is None:
        return None
    return VDW_RADII.get(str(element).strip().upper())


__all__ = [
    "VAN_DER_WAALS_RADII_ANGSTROM",
    "VDW_RADII",
    "VDW_RADII_VERSION",
    "vdw_radius_angstrom",
]
