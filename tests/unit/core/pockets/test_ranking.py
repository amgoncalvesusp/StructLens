from __future__ import annotations

from dataclasses import FrozenInstanceError

import pytest

from structlens.core.models import ResidueId
from structlens.core.pockets import AlphaSphere, PocketCandidate
from structlens.core.pockets.ranking import (
    PocketRankingComponents,
    rank_pocket_candidates,
    rank_pocket_candidates_with_components,
    score_pocket_candidate,
)


def _residue(auth_seq_id: str) -> ResidueId:
    return ResidueId("reference", "1", "A", auth_seq_id, None, "ALA")


def _sphere(
    label: str,
    center: tuple[float, float, float],
    radius: float,
    residues: tuple[str, ...],
) -> AlphaSphere:
    return AlphaSphere(
        center_xyz=center,
        radius_angstrom=radius,
        touching_atom_ids=tuple(f"{label}-atom-{index}" for index in range(1, 5)),
        lining_residues=tuple(_residue(value) for value in residues),
        source_simplex_atom_ids=tuple(f"{label}-atom-{index}" for index in range(1, 5)),
    )


def test_rank_pocket_candidates_orders_by_explicit_geometry_components() -> None:
    largest_radius = PocketCandidate(
        alpha_spheres=(
            _sphere("large-a", (0.0, 0.0, 0.0), 4.4, ("1", "2", "3", "4")),
            _sphere("large-b", (0.8, 0.0, 0.0), 4.0, ("2", "3", "4", "5")),
        )
    )
    richest_cluster = PocketCandidate(
        alpha_spheres=(
            _sphere("rich-a", (10.0, 0.0, 0.0), 3.2, ("10", "11", "12", "13")),
            _sphere("rich-b", (10.9, 0.0, 0.0), 3.1, ("11", "12", "13", "14")),
            _sphere("rich-c", (11.8, 0.0, 0.0), 3.0, ("12", "13", "14", "15")),
        )
    )
    smallest = PocketCandidate(alpha_spheres=(_sphere("small", (20.0, 0.0, 0.0), 3.0, ("20", "21", "22", "23")),))

    ranked = rank_pocket_candidates((largest_radius, smallest, richest_cluster), maximum_candidates=3)

    assert ranked == (richest_cluster, largest_radius, smallest)


def test_rank_pocket_candidates_respects_the_display_cap() -> None:
    candidates = tuple(
        PocketCandidate(alpha_spheres=(_sphere(f"candidate-{index}", (float(index), 0.0, 0.0), 3.0 + index, (str(index), "9", "8", "7")),))
        for index in range(4)
    )

    ranked = rank_pocket_candidates(candidates, maximum_candidates=2)

    assert len(ranked) == 2
    assert ranked[0].alpha_spheres[0].radius_angstrom > ranked[1].alpha_spheres[0].radius_angstrom


def test_rank_pocket_candidates_is_permutation_invariant_and_has_no_cap_by_default() -> None:
    candidates = tuple(
        PocketCandidate(
            alpha_spheres=(_sphere(f"candidate-{index}", (float(index), 0.0, 0.0), 3.0 + index, (str(index),)),)
        )
        for index in range(3)
    )

    first = rank_pocket_candidates(candidates)
    second = rank_pocket_candidates(tuple(reversed(candidates)))

    assert first == second
    assert len(first) == 3


def test_rank_pocket_candidates_exposes_independent_components() -> None:
    candidate = PocketCandidate(
        alpha_spheres=(
            _sphere("a", (0.0, 0.0, 0.0), 3.0, ("1", "2")),
            _sphere("b", (1.0, 0.0, 0.0), 4.0, ("2", "3")),
        )
    )

    components = score_pocket_candidate(candidate)
    ranked = rank_pocket_candidates_with_components((candidate,))

    assert components.to_json() == {
        "sphere_count": 2,
        "maximum_alpha_radius_angstrom": 4.0,
        "total_alpha_radius_angstrom": 7.0,
        "lining_residue_count": 3,
    }
    assert ranked[0].candidate == candidate
    assert ranked[0].components == components
    with pytest.raises(FrozenInstanceError):
        components.sphere_count = 3  # type: ignore[misc]


@pytest.mark.parametrize(
    ("kwargs", "message"),
    (
        ({"sphere_count": 0, "maximum_alpha_radius_angstrom": 1.0, "total_alpha_radius_angstrom": 1.0, "lining_residue_count": 0}, "sphere_count"),
        ({"sphere_count": 1, "maximum_alpha_radius_angstrom": 1.0, "total_alpha_radius_angstrom": 1.0, "lining_residue_count": -1}, "non-negative"),
        ({"sphere_count": 1, "maximum_alpha_radius_angstrom": float("nan"), "total_alpha_radius_angstrom": 1.0, "lining_residue_count": 0}, "finite"),
        ({"sphere_count": 1, "maximum_alpha_radius_angstrom": 1.0, "total_alpha_radius_angstrom": 0.0, "lining_residue_count": 0}, "positive"),
    ),
)
def test_ranking_components_reject_invalid_values(kwargs: dict[str, object], message: str) -> None:
    with pytest.raises(ValueError, match=message):
        PocketRankingComponents(**kwargs)  # type: ignore[arg-type]


def test_rank_pocket_candidates_rejects_invalid_candidate_and_caps() -> None:
    with pytest.raises(TypeError, match="PocketCandidate"):
        rank_pocket_candidates((object(),))  # type: ignore[arg-type]
    for cap in (0, -1, True, 1.5):
        with pytest.raises(ValueError, match="positive integer"):
            rank_pocket_candidates((), maximum_candidates=cap)  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="PocketCandidate"):
        score_pocket_candidate(object())  # type: ignore[arg-type]


def test_rank_pocket_candidates_handles_empty_input() -> None:
    assert rank_pocket_candidates(()) == ()
    assert rank_pocket_candidates_with_components(()) == ()
