from __future__ import annotations

from dataclasses import FrozenInstanceError

import pytest

from structlens.core.models import ResidueId
from structlens.core.pockets import AlphaSphere, PocketDetectionSettings, PocketGeometrySettings
from structlens.core.pockets.clustering import (
    PocketClusteringResult,
    cluster_alpha_spheres,
)
from structlens.core.pockets.delaunay import PocketAtom


def _residue(chain_id: str, auth_seq_id: str) -> ResidueId:
    return ResidueId("reference", "1", chain_id, auth_seq_id, None, "ALA")


def _atom(
    atom_id: str,
    coordinate: tuple[float, float, float],
    auth_seq_id: str,
) -> PocketAtom:
    return PocketAtom(atom_id, _residue("A", auth_seq_id), coordinate, 1.8)


def _dense_shell(*, center_x: float = 0.0, radius: float = 4.2, open_face: bool = False) -> tuple[PocketAtom, ...]:
    atoms: list[PocketAtom] = []
    index = 1
    for x in (-radius, 0.0, radius):
        for y in (-radius, 0.0, radius):
            for z in (-radius, 0.0, radius):
                if max(abs(x), abs(y), abs(z)) != radius:
                    continue
                if open_face and x == radius:
                    continue
                atoms.append(_atom(f"a-{center_x}-{index}", (center_x + x, y, z), str(index)))
                index += 1
    return tuple(atoms)


def _sphere(
    sphere_id: int,
    center: tuple[float, float, float],
    radius: float,
    residue_offset: int,
) -> AlphaSphere:
    return AlphaSphere(
        center_xyz=center,
        radius_angstrom=radius,
        touching_atom_ids=tuple(f"atom-{sphere_id}-{index}" for index in range(1, 5)),
        lining_residues=tuple(_residue("A", str(residue_offset + index)) for index in range(4)),
        source_simplex_atom_ids=tuple(f"atom-{sphere_id}-{index}" for index in range(1, 5)),
    )


def test_cluster_alpha_spheres_separates_two_buried_cavities() -> None:
    atoms = _dense_shell(center_x=-10.0) + _dense_shell(center_x=10.0)
    spheres = (
        _sphere(1, (-10.0, 0.0, 0.0), 4.2, 0),
        _sphere(2, (-9.2, 0.0, 0.0), 4.0, 10),
        _sphere(3, (10.0, 0.0, 0.0), 3.4, 20),
        _sphere(4, (10.8, 0.2, 0.0), 3.2, 30),
    )

    result = cluster_alpha_spheres(spheres, atoms, PocketDetectionSettings())

    assert result.diagnostics == ()
    assert len(result.candidates) == 2
    assert result.candidates[0].centroid_xyz[0] == pytest.approx(-9.6, abs=0.3)
    assert result.candidates[1].centroid_xyz[0] == pytest.approx(10.4, abs=0.3)


def test_cluster_alpha_spheres_rejects_solvent_exposed_clusters() -> None:
    spheres = (
        _sphere(1, (0.0, 0.0, 0.0), 4.0, 0),
        _sphere(2, (0.8, 0.0, 0.0), 4.0, 10),
    )
    result = cluster_alpha_spheres(
        spheres,
        _dense_shell(open_face=True),
        PocketDetectionSettings(solvent_grid_spacing_angstrom=0.75),
    )

    assert result.candidates == ()
    assert [item.code for item in result.diagnostics] == ["pocket.detect.solvent_exposed_cluster"]


def test_cluster_alpha_spheres_ignores_clusters_below_the_minimum_size() -> None:
    result = cluster_alpha_spheres(
        (_sphere(1, (0.0, 0.0, 0.0), 4.0, 0),),
        _dense_shell(),
        PocketDetectionSettings(minimum_cluster_size=2),
    )

    assert result.candidates == ()
    assert result.diagnostics == ()


def test_cluster_alpha_spheres_reports_a_configured_grid_resource_limit() -> None:
    result = cluster_alpha_spheres(
        (
            _sphere(1, (0.0, 0.0, 0.0), 4.0, 0),
            _sphere(2, (0.8, 0.0, 0.0), 4.0, 10),
        ),
        _dense_shell(),
        PocketDetectionSettings(max_solvent_grid_cells=10),
    )

    assert result.candidates == ()
    assert [item.code for item in result.diagnostics] == ["pocket.detect.resource_limit"]


def test_cluster_alpha_spheres_reports_excessive_raster_work_below_grid_cap() -> None:
    result = cluster_alpha_spheres(
        (
            _sphere(1, (0.0, 0.0, 0.0), 4.0, 0),
            _sphere(2, (0.8, 0.0, 0.0), 4.0, 10),
        ),
        _dense_shell(),
        PocketDetectionSettings(
            max_solvent_grid_cells=1_000_000,
            max_solvent_raster_cells=10,
        ),
    )

    assert result.candidates == ()
    assert [item.code for item in result.diagnostics] == ["pocket.detect.resource_limit"]


def test_cluster_exposure_checks_the_whole_candidate_region() -> None:
    result = cluster_alpha_spheres(
        (_sphere(1, (-1.25, 0.0, 0.0), 1.0, 0),),
        _dense_shell(open_face=True),
        PocketDetectionSettings(
            minimum_cluster_size=1,
            solvent_grid_spacing_angstrom=0.5,
        ),
    )

    assert result.candidates == ()
    assert [item.code for item in result.diagnostics] == ["pocket.detect.solvent_exposed_cluster"]


def test_cluster_alpha_spheres_empty_input_is_an_immutable_empty_result() -> None:
    result = cluster_alpha_spheres((), (), PocketDetectionSettings())

    assert result.candidates == ()
    assert result.diagnostics == ()
    with pytest.raises(FrozenInstanceError):
        result.candidates = ()  # type: ignore[misc]


def test_cluster_alpha_spheres_accepts_a_buried_candidate_without_atom_context() -> None:
    sphere = _sphere(1, (0.0, 0.0, 0.0), 4.0, 0)

    result = cluster_alpha_spheres((sphere,), (), PocketDetectionSettings(minimum_cluster_size=1))

    assert len(result.candidates) == 1
    assert result.diagnostics == ()


@pytest.mark.parametrize(
    ("spheres", "settings", "message"),
    (
        ((object(),), PocketDetectionSettings(), "AlphaSphere"),
        (("bad",), PocketDetectionSettings(), "AlphaSphere"),
        (("sentinel",), object(), "PocketDetectionSettings"),
    ),
)
def test_cluster_alpha_spheres_rejects_invalid_contracts(
    spheres: tuple[object, ...],
    settings: object,
    message: str,
) -> None:
    with pytest.raises((TypeError, ValueError), match=message):
        cluster_alpha_spheres(spheres, (), settings)  # type: ignore[arg-type]


def test_cluster_result_rejects_untyped_payloads() -> None:
    with pytest.raises(TypeError, match="PocketCandidate"):
        PocketClusteringResult(candidates=(object(),))  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="Diagnostic"):
        PocketClusteringResult(diagnostics=(object(),))  # type: ignore[arg-type]


def test_cluster_alpha_spheres_uses_nested_probe_setting_and_nonzero_padding() -> None:
    settings = PocketDetectionSettings(
        geometry=PocketGeometrySettings(probe_radius_angstrom=1.8),
        cluster_distance_padding_angstrom=0.2,
    )
    spheres = (
        _sphere(1, (-8.0, 0.0, 0.0), 2.0, 0),
        _sphere(2, (-4.1, 0.0, 0.0), 2.0, 10),
    )

    result = cluster_alpha_spheres(spheres, (), settings)

    assert len(result.candidates) == 1
    assert len(result.candidates[0].alpha_spheres) == 2


def test_cluster_alpha_spheres_rejects_nonfinite_atom_coordinates_before_grid_work() -> None:
    atom = _atom("invalid", (0.0, 0.0, 0.0), "1")
    object.__setattr__(atom, "coordinate", (float("nan"), 0.0, 0.0))

    with pytest.raises(ValueError, match="finite"):
        cluster_alpha_spheres((_sphere(1, (0.0, 0.0, 0.0), 4.0, 0),), (atom,), PocketDetectionSettings())


def test_cluster_alpha_spheres_orders_candidates_by_centroid_not_hash() -> None:
    spheres = (
        _sphere(1, (12.0, 0.0, 0.0), 2.0, 0),
        _sphere(2, (-12.0, 0.0, 0.0), 2.0, 10),
    )

    result = cluster_alpha_spheres(spheres, (), PocketDetectionSettings(minimum_cluster_size=1))

    assert tuple(item.centroid_xyz[0] for item in result.candidates) == (-12.0, 12.0)
