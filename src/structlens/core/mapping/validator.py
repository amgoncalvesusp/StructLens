"""Validation helpers for explicit residue maps and AUTO policy."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, replace

from structlens.core.errors import MappingError
from structlens.core.models import (
    AlignmentMode,
    AnalysisSettings,
    ResidueCorrespondence,
)


@dataclass(frozen=True, slots=True)
class MappingDecision:
    mode: AlignmentMode
    reason: str


def choose_auto_mapping(
    identity: float,
    coverage: float,
    settings: AnalysisSettings | None = None,
) -> MappingDecision:
    settings = settings or AnalysisSettings()
    if (
        identity >= settings.minimum_sequence_identity
        and coverage >= settings.minimum_sequence_coverage
    ):
        return MappingDecision(
            AlignmentMode.SEQUENCE,
            (
                f"AUTO selected sequence-guided mapping: identity={identity:.3f} "
                f">= {settings.minimum_sequence_identity:.3f}, coverage={coverage:.3f} "
                f">= {settings.minimum_sequence_coverage:.3f}"
            ),
        )
    return MappingDecision(
        AlignmentMode.STRUCTURE,
        (
            f"AUTO selected structure-guided mapping: identity={identity:.3f} or "
            f"coverage={coverage:.3f} below configured thresholds"
        ),
    )


def validate_correspondence(
    correspondences: Iterable[ResidueCorrespondence],
    *,
    displacement_cutoff_angstrom: float | None = None,
) -> list[ResidueCorrespondence]:
    """Validate uniqueness and mark geometric outliers without deleting rows."""

    values = list(correspondences)
    reference_seen: set[object] = set()
    target_seen: set[object] = set()
    for item in values:
        if (
            item.reference is not None
            and item.reference in reference_seen
            and not item.mapping_locked
        ):
            raise MappingError(
                "A reference residue appears more than once in the correspondence map"
            )
        if (
            item.target is not None
            and item.target in target_seen
            and not item.mapping_locked
        ):
            raise MappingError(
                "A target residue appears more than once in the correspondence map"
            )
        if item.reference is not None:
            reference_seen.add(item.reference)
        if item.target is not None:
            target_seen.add(item.target)
    if displacement_cutoff_angstrom is None:
        return values
    if displacement_cutoff_angstrom <= 0:
        raise ValueError("displacement_cutoff_angstrom must be positive")
    return [
        replace(
            item,
            is_outlier=(
                item.ca_displacement_angstrom is not None
                and item.ca_displacement_angstrom > displacement_cutoff_angstrom
            ),
        )
        for item in values
    ]


__all__ = ["MappingDecision", "choose_auto_mapping", "validate_correspondence"]
