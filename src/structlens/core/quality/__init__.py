"""Typed structural coordinate-quality checks."""

from .clashes import (
    CLASH_ALGORITHM_VERSION,
    VDW_RADII,
    VDW_RADII_VERSION,
    ClashAtom,
    ClashDiagnostic,
    ClashObservation,
    ClashScreeningResult,
    screen_heavy_atom_overlaps,
)
from .coordinates import (
    check_coordinate_quality,
    check_residue_quality,
    diagnostics_for_atoms,
    vdw_radius_angstrom,
)
from .geometry import screen_ca_pseudo_geometry, screen_chain_geometry, screen_peptide_continuity
from .models import CoordinateAtom, CoordinateQCSettings, CoordinateResidue, StructureQualityReport

__all__ = [
    "CLASH_ALGORITHM_VERSION",
    "ClashAtom",
    "ClashDiagnostic",
    "ClashObservation",
    "ClashScreeningResult",
    "CoordinateAtom",
    "CoordinateQCSettings",
    "CoordinateResidue",
    "StructureQualityReport",
    "VDW_RADII",
    "VDW_RADII_VERSION",
    "check_coordinate_quality",
    "check_residue_quality",
    "diagnostics_for_atoms",
    "screen_ca_pseudo_geometry",
    "screen_chain_geometry",
    "screen_heavy_atom_overlaps",
    "screen_peptide_continuity",
    "vdw_radius_angstrom",
]
