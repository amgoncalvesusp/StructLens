"""Pure visualization state; scientific analysis objects remain untouched."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, replace
from enum import Enum

from structlens.core.models import CorrespondenceStatus, ResidueCorrespondence


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
        correspondences: Sequence[ResidueCorrespondence],
        state: VisualizationState,
    ) -> tuple[ResidueCorrespondence, ...]:
        return tuple(
            item for item in correspondences if _matches(item, state.highlight_filter)
        )

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
            "Active site": VisualizationState(
                preset="Active site", highlight_filter=HighlightFilter.KEY
            ),
            "Presentation": VisualizationState(
                preset="Presentation",
                representation=Representation.CARTOON_STICKS,
                show_labels=True,
            ),
        }
        if preset not in presets:
            raise ValueError(f"Unknown StructLens visualization preset: {preset}")
        return replace(presets[preset])


def _matches(item: ResidueCorrespondence, filter_value: HighlightFilter) -> bool:
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
]
