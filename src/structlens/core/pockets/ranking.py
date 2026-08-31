"""Transparent, deterministic ranking of detected pocket candidates."""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

from .models import PocketCandidate


@dataclass(frozen=True, slots=True)
class PocketRankingComponents:
    """Displayed geometric components used to order a candidate.

    These are intentionally not collapsed into a probability.  Consumers can
    show each component independently and preserve the method's audit trail.
    """

    sphere_count: int
    maximum_alpha_radius_angstrom: float
    total_alpha_radius_angstrom: float
    lining_residue_count: int

    def __post_init__(self) -> None:
        if self.sphere_count <= 0 or self.lining_residue_count < 0:
            raise ValueError("ranking counts must be non-negative and sphere_count must be positive")
        for value, name in (
            (self.maximum_alpha_radius_angstrom, "maximum_alpha_radius_angstrom"),
            (self.total_alpha_radius_angstrom, "total_alpha_radius_angstrom"),
        ):
            if not math.isfinite(float(value)) or float(value) <= 0.0:
                raise ValueError(f"{name} must be finite and positive")

    def to_json(self) -> dict[str, int | float]:
        return {
            "sphere_count": self.sphere_count,
            "maximum_alpha_radius_angstrom": self.maximum_alpha_radius_angstrom,
            "total_alpha_radius_angstrom": self.total_alpha_radius_angstrom,
            "lining_residue_count": self.lining_residue_count,
        }


@dataclass(frozen=True, slots=True)
class RankedPocketCandidate:
    """A candidate plus the independent components shown in the UI."""

    candidate: PocketCandidate
    components: PocketRankingComponents


def score_pocket_candidate(candidate: PocketCandidate) -> PocketRankingComponents:
    """Calculate explicit geometric ranking components for one candidate."""

    if not isinstance(candidate, PocketCandidate):
        raise TypeError("candidate must be PocketCandidate")
    radii = tuple(float(sphere.radius_angstrom) for sphere in candidate.alpha_spheres)
    return PocketRankingComponents(
        sphere_count=len(radii),
        maximum_alpha_radius_angstrom=max(radii),
        total_alpha_radius_angstrom=sum(radii),
        lining_residue_count=len(candidate.lining_residues),
    )


def _rank_key(candidate: PocketCandidate) -> tuple[int, float, float, int, str]:
    components = score_pocket_candidate(candidate)
    return (
        -components.sphere_count,
        -components.maximum_alpha_radius_angstrom,
        -components.total_alpha_radius_angstrom,
        -components.lining_residue_count,
        candidate.candidate_id,
    )


def rank_pocket_candidates(
    candidates: Sequence[PocketCandidate],
    *,
    maximum_candidates: int | None = None,
) -> tuple[PocketCandidate, ...]:
    """Return candidates in stable display order, applying an optional cap."""

    values = tuple(candidates)
    if any(not isinstance(item, PocketCandidate) for item in values):
        raise TypeError("candidates must contain PocketCandidate values")
    if maximum_candidates is not None:
        if (
            isinstance(maximum_candidates, bool)
            or not isinstance(maximum_candidates, int)
            or maximum_candidates <= 0
        ):
            raise ValueError("maximum_candidates must be a positive integer")
    ranked = tuple(sorted(values, key=_rank_key))
    return ranked if maximum_candidates is None else ranked[:maximum_candidates]


def rank_pocket_candidates_with_components(
    candidates: Sequence[PocketCandidate],
    *,
    maximum_candidates: int | None = None,
) -> tuple[RankedPocketCandidate, ...]:
    """Rank candidates while retaining the independent displayed components."""

    ranked = rank_pocket_candidates(candidates, maximum_candidates=maximum_candidates)
    return tuple(RankedPocketCandidate(item, score_pocket_candidate(item)) for item in ranked)


__all__ = [
    "PocketRankingComponents",
    "RankedPocketCandidate",
    "rank_pocket_candidates",
    "rank_pocket_candidates_with_components",
    "score_pocket_candidate",
]
