"""Active-site definition and descriptive geometry metric contracts."""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum

from structlens.core.models import ResidueId


class SiteDefinitionMode(str, Enum):
    KEY_RESIDUES = "key_residues"
    LIGAND_RADIUS = "ligand_radius"
    RESIDUE_RADIUS = "residue_radius"


SiteDefinitionKind = SiteDefinitionMode
SiteType = SiteDefinitionMode


def _fraction(value: float | None, name: str) -> None:
    if value is not None and (not math.isfinite(value) or not 0.0 <= value <= 1.0):
        raise ValueError(f"{name} must be finite and between 0 and 1")


def _distance(value: float | None, name: str) -> None:
    if value is not None and (not math.isfinite(value) or value < 0.0):
        raise ValueError(f"{name} must be finite and non-negative")


@dataclass(frozen=True, slots=True, init=False)
class SiteDefinition:
    site_id: str
    name: str
    mode: SiteDefinitionMode
    reference_residues: tuple[ResidueId, ...]
    center_residue: ResidueId | None
    ligand_id: str | None
    radius_angstrom: float | None

    def __init__(
        self,
        site_id: str,
        name: str | None = None,
        mode: SiteDefinitionMode | None = None,
        reference_residues: tuple[ResidueId, ...] = (),
        center_residue: ResidueId | None = None,
        ligand_id: str | None = None,
        radius_angstrom: float | None = None,
        **legacy: object,
    ) -> None:
        if name is None:
            name = site_id
        if mode is None and "kind" in legacy:
            mode = SiteDefinitionMode(legacy.pop("kind"))
        if "key_residues" in legacy:
            reference_residues = tuple(legacy.pop("key_residues"))  # type: ignore[arg-type]
        if legacy:
            raise TypeError(f"unexpected SiteDefinition fields: {', '.join(sorted(legacy))}")
        if not site_id or not name:
            raise ValueError("site_id and name must not be empty")
        if mode is None:
            raise ValueError("mode must be specified")
        mode = SiteDefinitionMode(mode)
        residues = tuple(reference_residues)
        if any(not isinstance(item, ResidueId) for item in residues):
            raise TypeError("reference_residues must contain ResidueId values")
        if mode is SiteDefinitionMode.KEY_RESIDUES and not residues:
            raise ValueError("key_residues sites require reference_residues")
        if mode is SiteDefinitionMode.LIGAND_RADIUS and not ligand_id:
            raise ValueError("ligand_radius sites require ligand_id")
        if mode in {SiteDefinitionMode.LIGAND_RADIUS, SiteDefinitionMode.RESIDUE_RADIUS} and radius_angstrom is None:
            raise ValueError("radius_angstrom is required for radius sites")
        _distance(radius_angstrom, "radius_angstrom")
        object.__setattr__(self, "site_id", site_id)
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "mode", mode)
        object.__setattr__(self, "reference_residues", residues)
        object.__setattr__(self, "center_residue", center_residue)
        object.__setattr__(self, "ligand_id", ligand_id)
        object.__setattr__(self, "radius_angstrom", radius_angstrom)

    @property
    def kind(self) -> SiteDefinitionMode:
        return self.mode

    @property
    def key_residues(self) -> tuple[ResidueId, ...]:
        return self.reference_residues


@dataclass(frozen=True, slots=True)
class SiteMetrics:
    site_id: str
    structure_id: str
    mapped_residue_count: int
    coverage_fraction: float
    global_frame_backbone_rmsd_angstrom: float | None = None
    site_fitted_backbone_rmsd_angstrom: float | None = None
    centroid_displacement_angstrom: float | None = None
    radius_of_gyration_angstrom: float | None = None
    atomic_envelope_volume_angstrom3: float | None = None
    sasa_angstrom2: float | None = None
    polar_residue_fraction: float | None = None
    charged_residue_fraction: float | None = None

    def __post_init__(self) -> None:
        if not self.site_id or not self.structure_id:
            raise ValueError("site_id and structure_id must not be empty")
        if self.mapped_residue_count < 0:
            raise ValueError("mapped_residue_count must be non-negative")
        _fraction(self.coverage_fraction, "coverage_fraction")
        for name in (
            "global_frame_backbone_rmsd_angstrom",
            "site_fitted_backbone_rmsd_angstrom",
            "centroid_displacement_angstrom",
            "radius_of_gyration_angstrom",
            "atomic_envelope_volume_angstrom3",
            "sasa_angstrom2",
        ):
            _distance(getattr(self, name), name)
        for name in ("polar_residue_fraction", "charged_residue_fraction"):
            _fraction(getattr(self, name), name)

    @property
    def global_frame_rmsd_angstrom(self) -> float | None:
        return self.global_frame_backbone_rmsd_angstrom

    @property
    def site_fitted_rmsd_angstrom(self) -> float | None:
        return self.site_fitted_backbone_rmsd_angstrom


__all__ = ["SiteDefinition", "SiteDefinitionKind", "SiteDefinitionMode", "SiteMetrics", "SiteType"]
