"""Pure visualization state; scientific analysis objects remain untouched."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from enum import Enum
from typing import Protocol, TypeVar

from structlens.core.models import CorrespondenceStatus


class VisualizableCorrespondence(Protocol):
    status: CorrespondenceStatus
    is_key_residue: bool
    is_outlier: bool
    ca_displacement_angstrom: float | None


_Correspondence = TypeVar("_Correspondence", bound=VisualizableCorrespondence)


class HighlightFilter(str, Enum):
    ALL = "all"
    MUTATIONS = "mutations"
    CONSERVED = "conserved"
    KEY = "key"
    MUTATED_KEY = "mutated_key"
    OUTLIERS = "outliers"
    INSERTIONS_DELETIONS = "insertions_deletions"
    DISPLACEMENT = "displacement"


class ColorMode(str, Enum):
    REFERENCE_TARGET = "reference_target"
    MUTATION_STATUS = "mutation_status"
    CA_DISPLACEMENT = "ca_displacement"
    BACKBONE_RMSD = "backbone_rmsd"


class Representation(str, Enum):
    STICKS = "sticks"
    SPHERES = "spheres"
    STICKS_SPHERES = "sticks_spheres"
    CARTOON_STICKS = "cartoon_sticks"
    SURFACE = "surface"
    LABELS = "labels"


@dataclass(frozen=True, slots=True)
class VisualizationState:
    highlight_filter: HighlightFilter = HighlightFilter.ALL
    color_mode: ColorMode = ColorMode.MUTATION_STATUS
    representation: Representation = Representation.STICKS
    show_labels: bool = False
    show_reference: bool = True
    show_target: bool = True
    local_radius_angstrom: float = 5.0
    preset: str = "Minimal"


class VisualizationRenderer:
    def filtered_correspondences(
        self,
        correspondences: Sequence[_Correspondence],
        state: VisualizationState,
    ) -> tuple[_Correspondence, ...]:
        return tuple(item for item in correspondences if _matches(item, state.highlight_filter))

    def apply_preset(self, preset: str) -> VisualizationState:
        presets = {
            "Minimal": VisualizationState(preset="Minimal"),
            "Publication": VisualizationState(
                preset="Publication",
                representation=Representation.CARTOON_STICKS,
                color_mode=ColorMode.REFERENCE_TARGET,
            ),
            "Mutation focus": VisualizationState(
                preset="Mutation focus",
                highlight_filter=HighlightFilter.MUTATIONS,
                representation=Representation.STICKS_SPHERES,
            ),
            "Structural deviation": VisualizationState(
                preset="Structural deviation",
                highlight_filter=HighlightFilter.DISPLACEMENT,
                color_mode=ColorMode.CA_DISPLACEMENT,
            ),
            "Active site": VisualizationState(preset="Active site", highlight_filter=HighlightFilter.KEY),
            "Presentation": VisualizationState(
                preset="Presentation",
                representation=Representation.CARTOON_STICKS,
                show_labels=True,
            ),
        }
        if preset not in presets:
            raise ValueError(f"Unknown StructLens visualization preset: {preset}")
        return replace(presets[preset])


def visualization_state_from_mapping(payload: object) -> VisualizationState:
    """Parse persisted visualization controls into a closed typed value."""

    if payload is None:
        payload = {}
    if not isinstance(payload, Mapping):
        raise TypeError("visualization state must be an object")

    def boolean(key: str, default: bool) -> bool:
        value = payload.get(key, default)
        if not isinstance(value, bool):
            raise TypeError(f"{key} must be a boolean")
        return value

    radius_value = payload.get("local_radius_angstrom", 5.0)
    if isinstance(radius_value, bool):
        raise TypeError("local_radius_angstrom must be numeric")
    radius = float(radius_value)
    if not math.isfinite(radius) or not 0.1 <= radius <= 20.0:
        raise ValueError("local_radius_angstrom must be finite and between 0.1 and 20.0")
    preset = str(payload.get("preset", "Minimal"))
    allowed_presets = {
        "Minimal",
        "Publication",
        "Mutation focus",
        "Structural deviation",
        "Active site",
        "Presentation",
    }
    if preset not in allowed_presets:
        raise ValueError("preset is not supported")
    return VisualizationState(
        highlight_filter=HighlightFilter(str(payload.get("highlight_filter", HighlightFilter.ALL.value))),
        color_mode=ColorMode(str(payload.get("color_mode", ColorMode.MUTATION_STATUS.value))),
        representation=Representation(str(payload.get("representation", Representation.STICKS.value))),
        show_labels=boolean("show_labels", False),
        show_reference=boolean("show_reference", True),
        show_target=boolean("show_target", True),
        local_radius_angstrom=radius,
        preset=preset,
    )


def _matches(item: VisualizableCorrespondence, filter_value: HighlightFilter) -> bool:
    if filter_value is HighlightFilter.ALL:
        return True
    if filter_value is HighlightFilter.MUTATIONS:
        return item.status not in {
            CorrespondenceStatus.CONSERVED,
            CorrespondenceStatus.UNMAPPED,
        }
    if filter_value is HighlightFilter.CONSERVED:
        return item.status is CorrespondenceStatus.CONSERVED
    if filter_value is HighlightFilter.KEY:
        return item.is_key_residue
    if filter_value is HighlightFilter.MUTATED_KEY:
        return item.is_key_residue and item.status is not CorrespondenceStatus.CONSERVED
    if filter_value is HighlightFilter.OUTLIERS:
        return item.is_outlier
    if filter_value is HighlightFilter.INSERTIONS_DELETIONS:
        return item.status in {
            CorrespondenceStatus.INSERTION,
            CorrespondenceStatus.DELETION,
        }
    if filter_value is HighlightFilter.DISPLACEMENT:
        return item.ca_displacement_angstrom is not None
    return False


__all__ = [
    "ColorMode",
    "HighlightFilter",
    "Representation",
    "VisualizationRenderer",
    "VisualizationState",
    "visualization_state_from_mapping",
]
