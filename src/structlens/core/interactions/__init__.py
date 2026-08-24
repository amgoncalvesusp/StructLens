"""Reference-normalized structural interaction contracts."""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum

from structlens.core.models import ResidueId


class InteractionType(str, Enum):
    HBOND_GEOMETRIC = "hbond_geometric"
    SALT_BRIDGE = "salt_bridge"
    HYDROPHOBIC = "hydrophobic"
    PI_STACKING = "pi_stacking"
    CATION_PI = "cation_pi"
    METAL_CONTACT = "metal_contact"
    HYDROGEN_BOND = "hbond_geometric"


class InteractionChange(str, Enum):
    CONSERVED = "conserved"
    GAINED = "gained"
    LOST = "lost"
    TARGET_ONLY_UNMAPPED = "target_only_unmapped"


InteractionChangeKind = InteractionChange
InteractionDifferenceKind = InteractionChange
InteractionDifferenceStatus = InteractionChange


def _finite(value: float, name: str, *, nonnegative: bool = False) -> None:
    if not math.isfinite(value) or (nonnegative and value < 0.0):
        raise ValueError(f"{name} must be finite" + (" and non-negative" if nonnegative else ""))


@dataclass(frozen=True, slots=True)
class InteractionRecord:
    structure_id: str
    interaction_type: InteractionType
    residue_a: ResidueId
    residue_b: ResidueId | None
    atom_a: str | None
    atom_b: str | None
    distance_angstrom: float
    angle_degrees: float | None = None
    ligand_or_metal_id: str | None = None
    evidence_mode: str = "heavy_atom_geometry"

    def __post_init__(self) -> None:
        if not self.structure_id:
            raise ValueError("structure_id must not be empty")
        if not isinstance(self.interaction_type, InteractionType):
            object.__setattr__(self, "interaction_type", InteractionType(self.interaction_type))
        _finite(self.distance_angstrom, "distance_angstrom", nonnegative=True)
        if self.angle_degrees is not None:
            _finite(self.angle_degrees, "angle_degrees")
        if not self.evidence_mode:
            raise ValueError("evidence_mode must not be empty")

    @property
    def is_putative(self) -> bool:
        return self.evidence_mode == "heavy_atom_geometry"


@dataclass(frozen=True, slots=True)
class ReferenceInteractionKey:
    interaction_type: InteractionType
    reference_position_a: str
    reference_position_b: str | None = None
    external_partner_id: str | None = None

    def __post_init__(self) -> None:
        if not self.reference_position_a:
            raise ValueError("reference_position_a must not be empty")
        if not isinstance(self.interaction_type, InteractionType):
            object.__setattr__(self, "interaction_type", InteractionType(self.interaction_type))


@dataclass(frozen=True, slots=True)
class InteractionDifference:
    key: ReferenceInteractionKey
    change: InteractionChange
    reference_record: InteractionRecord | None = None
    target_record: InteractionRecord | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.change, InteractionChange):
            object.__setattr__(self, "change", InteractionChange(self.change))
        if self.change is InteractionChange.CONSERVED and (self.reference_record is None or self.target_record is None):
            raise ValueError("conserved interaction differences require both records")
        if self.change is InteractionChange.GAINED and self.target_record is None:
            raise ValueError("gained interaction differences require target_record")
        if self.change is InteractionChange.LOST and self.reference_record is None:
            raise ValueError("lost interaction differences require reference_record")

    @property
    def reference(self) -> InteractionRecord | None:
        return self.reference_record

    @property
    def target(self) -> InteractionRecord | None:
        return self.target_record

    @property
    def kind(self) -> InteractionChange:
        return self.change


InteractionKey = ReferenceInteractionKey

__all__ = [
    "InteractionChange",
    "InteractionChangeKind",
    "InteractionDifference",
    "InteractionDifferenceKind",
    "InteractionDifferenceStatus",
    "InteractionRecord",
    "InteractionType",
    "ReferenceInteractionKey",
    "InteractionKey",
]
