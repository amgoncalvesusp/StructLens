from __future__ import annotations

import json
from dataclasses import FrozenInstanceError

import pytest

from structlens.core.models import ResidueId
from structlens.core.pockets import (
    AlphaSphere,
    PocketCandidate,
    PocketDetectionSettings,
    PocketGeometrySettings,
)


def _residue(chain_id: str, auth_seq_id: str) -> ResidueId:
    return ResidueId("reference", "1", chain_id, auth_seq_id, None, "ALA")


def test_pocket_geometry_settings_validate_units_and_ranges() -> None:
    settings = PocketGeometrySettings(
        minimum_alpha_sphere_radius_angstrom=3.0,
        maximum_alpha_sphere_radius_angstrom=6.0,
        probe_radius_angstrom=1.4,
        lining_contact_slack_angstrom=0.5,
    )

    assert settings.to_json() == {
        "minimum_alpha_sphere_radius_angstrom": 3.0,
        "maximum_alpha_sphere_radius_angstrom": 6.0,
        "probe_radius_angstrom": 1.4,
        "lining_contact_slack_angstrom": 0.5,
    }

    with pytest.raises(ValueError, match="minimum_alpha_sphere_radius_angstrom"):
        PocketGeometrySettings(minimum_alpha_sphere_radius_angstrom=-1.0)
    with pytest.raises(ValueError, match="must not exceed"):
        PocketGeometrySettings(
            minimum_alpha_sphere_radius_angstrom=5.0,
            maximum_alpha_sphere_radius_angstrom=4.0,
        )


def test_alpha_sphere_canonicalizes_ordering_and_serialization() -> None:
    sphere = AlphaSphere(
        center_xyz=(1.0, 2.0, 3.0),
        radius_angstrom=4.0,
        touching_atom_ids=("atom-3", "atom-1", "atom-4", "atom-2"),
        lining_residues=(_residue("B", "20"), _residue("A", "10"), _residue("A", "11")),
        source_simplex_atom_ids=("atom-4", "atom-2", "atom-1", "atom-3"),
    )

    payload = sphere.to_json()

    assert payload["touching_atom_ids"] == ["atom-1", "atom-2", "atom-3", "atom-4"]
    assert payload["lining_residues"][0]["chain_id"] == "A"
    assert payload["source_simplex_atom_ids"] == ["atom-1", "atom-2", "atom-3", "atom-4"]
    assert len(payload["sphere_id"]) == 64
    assert json.loads(sphere.canonical_json_bytes()) == payload
    with pytest.raises(FrozenInstanceError):
        sphere.radius_angstrom = 5.0  # type: ignore[misc]


def test_pocket_candidate_identity_is_deterministic_under_input_reordering() -> None:
    sphere_a = AlphaSphere(
        center_xyz=(0.0, 0.0, 0.0),
        radius_angstrom=3.0,
        touching_atom_ids=("atom-1", "atom-2", "atom-3", "atom-4"),
        lining_residues=(_residue("A", "10"), _residue("A", "11")),
        source_simplex_atom_ids=("atom-1", "atom-2", "atom-3", "atom-4"),
    )
    sphere_b = AlphaSphere(
        center_xyz=(2.0, 0.0, 0.0),
        radius_angstrom=3.5,
        touching_atom_ids=("atom-5", "atom-6", "atom-7", "atom-8"),
        lining_residues=(_residue("A", "11"), _residue("B", "20")),
        source_simplex_atom_ids=("atom-5", "atom-6", "atom-7", "atom-8"),
    )

    first = PocketCandidate(alpha_spheres=(sphere_b, sphere_a))
    second = PocketCandidate(alpha_spheres=(sphere_a, sphere_b))

    assert first.candidate_id == second.candidate_id
    assert first.canonical_json_bytes() == second.canonical_json_bytes()
    assert first.lining_residues == (_residue("A", "10"), _residue("A", "11"), _residue("B", "20"))
    assert first.to_json()["centroid_xyz"] == [1.0, 0.0, 0.0]


@pytest.mark.parametrize(
    ("kwargs", "message"),
    (
        ({"minimum_alpha_sphere_radius_angstrom": float("nan")}, "finite"),
        ({"lining_contact_slack_angstrom": -0.1}, "non-negative"),
    ),
)
def test_pocket_geometry_settings_reject_nonfinite_and_negative_values(
    kwargs: dict[str, float],
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        PocketGeometrySettings(**kwargs)


def test_pocket_detection_settings_validate_caps_and_nested_geometry() -> None:
    settings = PocketDetectionSettings(
        geometry=PocketGeometrySettings(minimum_alpha_sphere_radius_angstrom=3.0),
        cluster_distance_padding_angstrom=1.25,
        solvent_grid_spacing_angstrom=0.75,
        solvent_boundary_margin_angstrom=2.5,
        minimum_cluster_size=2,
        maximum_candidates=4,
        max_atom_count=200,
        max_estimated_simplices=5000,
        max_solvent_grid_cells=6000,
        max_clearance_atom_checks=7000,
        max_solvent_raster_cells=8000,
    )

    assert settings.to_json() == {
        "geometry": {
            "minimum_alpha_sphere_radius_angstrom": 3.0,
            "maximum_alpha_sphere_radius_angstrom": 6.2,
            "probe_radius_angstrom": 1.4,
            "lining_contact_slack_angstrom": 0.5,
        },
        "cluster_distance_padding_angstrom": 1.25,
        "solvent_grid_spacing_angstrom": 0.75,
        "solvent_boundary_margin_angstrom": 2.5,
        "minimum_cluster_size": 2,
        "maximum_candidates": 4,
        "max_atom_count": 200,
        "max_estimated_simplices": 5000,
        "max_solvent_grid_cells": 6000,
        "max_clearance_atom_checks": 7000,
        "max_solvent_raster_cells": 8000,
    }


@pytest.mark.parametrize(
    ("kwargs", "message"),
    (
        ({"cluster_distance_padding_angstrom": -0.1}, "non-negative"),
        ({"solvent_grid_spacing_angstrom": 0.0}, "positive"),
        ({"minimum_cluster_size": 0}, "positive integer"),
        ({"maximum_candidates": 0}, "positive integer"),
        ({"max_atom_count": 0}, "positive integer"),
        ({"max_estimated_simplices": 0}, "positive integer"),
        ({"max_solvent_grid_cells": 0}, "positive integer"),
        ({"max_clearance_atom_checks": 0}, "positive integer"),
        ({"max_solvent_raster_cells": 0}, "positive integer"),
    ),
)
def test_pocket_detection_settings_reject_invalid_ranges(
    kwargs: dict[str, object],
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        PocketDetectionSettings(**kwargs)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    (
        ({"center_xyz": ("bad", 2.0, 3.0)}, "numeric"),
        ({"center_xyz": (1.0, 2.0)}, "three finite"),
        ({"touching_atom_ids": ("atom-1", "atom-1", "atom-2", "atom-3")}, "duplicates"),
        ({"touching_atom_ids": ("atom-1", "atom-2", "atom-3")}, "exactly 4"),
        ({"touching_atom_ids": ("atom-1", " ", "atom-3", "atom-4")}, "non-empty atom IDs"),
        ({"source_simplex_atom_ids": ("atom-1", "atom-2", "atom-3")}, "exactly 4"),
        ({"lining_residues": (object(),)}, "ResidueId"),
    ),
)
def test_alpha_sphere_rejects_invalid_geometry_contracts(
    kwargs: dict[str, object],
    message: str,
) -> None:
    defaults: dict[str, object] = {
        "center_xyz": (1.0, 2.0, 3.0),
        "radius_angstrom": 4.0,
        "touching_atom_ids": ("atom-1", "atom-2", "atom-3", "atom-4"),
        "lining_residues": (_residue("A", "10A"),),
        "source_simplex_atom_ids": ("atom-1", "atom-2", "atom-3", "atom-4"),
    }
    defaults.update(kwargs)
    with pytest.raises(ValueError, match=message):
        AlphaSphere(**defaults)  # type: ignore[arg-type]


def test_pocket_candidate_rejects_empty_or_untyped_sphere_collections() -> None:
    with pytest.raises(ValueError, match="AlphaSphere"):
        PocketCandidate(alpha_spheres=(object(),))  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="AlphaSphere"):
        PocketCandidate(alpha_spheres=())


def test_pocket_candidate_lineage_is_validated_without_changing_geometry_identity() -> None:
    sphere = AlphaSphere(
        center_xyz=(0.0, 0.0, 0.0),
        radius_angstrom=3.0,
        touching_atom_ids=("atom-1", "atom-2", "atom-3", "atom-4"),
        lining_residues=(_residue("A", "10"),),
        source_simplex_atom_ids=("atom-1", "atom-2", "atom-3", "atom-4"),
    )
    geometric = PocketCandidate((sphere,))
    bound = PocketCandidate(
        (sphere,),
        source_content_id="a" * 64,
        selection_id="b" * 64,
    )

    assert bound.candidate_id == geometric.candidate_id
    assert bound.to_json()["lineage"] == {
        "source_content_id": "a" * 64,
        "selection_id": "b" * 64,
    }
    with pytest.raises(ValueError, match="together"):
        PocketCandidate((sphere,), source_content_id="a" * 64)
    with pytest.raises(ValueError, match="SHA-256"):
        PocketCandidate((sphere,), source_content_id="bad", selection_id="b" * 64)
