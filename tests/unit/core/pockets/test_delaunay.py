from __future__ import annotations

import math
from dataclasses import FrozenInstanceError

import numpy as np
import pytest
from scipy.spatial import QhullError

from structlens.core.evidence import Availability, Diagnostic
from structlens.core.models import ResidueId
from structlens.core.pockets import PocketDetectionSettings, PocketGeometrySettings
from structlens.core.pockets.delaunay import (
    PocketAtom,
    alpha_sphere_from_simplex,
    detect_alpha_spheres,
)


def _atom(
    atom_id: str,
    coordinate: tuple[float, float, float],
    auth_seq_id: str,
    *,
    radius: float = 1.8,
) -> PocketAtom:
    return PocketAtom(
        atom_id=atom_id,
        residue_id=ResidueId("reference", "1", "A", auth_seq_id, None, "ALA"),
        coordinate=coordinate,
        radius_angstrom=radius,
    )


def _regular_tetrahedron(radius: float = 4.0) -> tuple[PocketAtom, ...]:
    scale = radius / math.sqrt(3.0)
    vertices = (
        (scale, scale, scale),
        (-scale, -scale, scale),
        (-scale, scale, -scale),
        (scale, -scale, -scale),
    )
    return tuple(_atom(f"atom-{index}", vertex, str(index)) for index, vertex in enumerate(vertices, start=1))


def _dense_shell(*, radius: float = 4.2, open_face: bool = False) -> tuple[PocketAtom, ...]:
    atoms: list[PocketAtom] = []
    index = 1
    for x in (-radius, 0.0, radius):
        for y in (-radius, 0.0, radius):
            for z in (-radius, 0.0, radius):
                if max(abs(x), abs(y), abs(z)) != radius:
                    continue
                if open_face and x == radius:
                    continue
                atoms.append(_atom(f"shell-{index}", (x, y, z), str(index)))
                index += 1
    return tuple(atoms)


def test_alpha_sphere_from_simplex_rejects_nonvertex_intrusion() -> None:
    simplex = _regular_tetrahedron()
    intruder = _atom("intruder", (0.0, 0.0, 2.0), "99", radius=1.8)
    settings = PocketGeometrySettings(
        minimum_alpha_sphere_radius_angstrom=2.0,
        maximum_alpha_sphere_radius_angstrom=3.0,
    )

    accepted = alpha_sphere_from_simplex(simplex, simplex, settings)
    rejected = alpha_sphere_from_simplex(simplex, simplex + (intruder,), settings)

    assert accepted is not None
    assert accepted.center_xyz == pytest.approx((0.0, 0.0, 0.0), abs=1e-12)
    assert accepted.radius_angstrom == pytest.approx(2.2)
    assert rejected is None


def test_alpha_sphere_filters_on_vdw_clearance_not_atom_center_radius() -> None:
    simplex = _regular_tetrahedron()

    result = alpha_sphere_from_simplex(
        simplex,
        simplex,
        PocketGeometrySettings(
            minimum_alpha_sphere_radius_angstrom=2.8,
            maximum_alpha_sphere_radius_angstrom=6.2,
        ),
    )

    assert result is None


def test_surface_tangent_nonvertex_atom_is_retained_as_lining_evidence() -> None:
    simplex = _regular_tetrahedron()
    tangent = _atom("tangent", (4.2, 0.0, 0.0), "99", radius=2.0)

    result = alpha_sphere_from_simplex(
        simplex,
        simplex + (tangent,),
        PocketGeometrySettings(
            minimum_alpha_sphere_radius_angstrom=2.0,
            maximum_alpha_sphere_radius_angstrom=3.0,
            lining_contact_slack_angstrom=0.1,
        ),
    )

    assert result is not None
    assert tangent.residue_id in result.lining_residues


def test_pocket_atom_normalizes_known_radius_identity_and_is_immutable() -> None:
    atom = PocketAtom(
        " atom-1 ",
        ResidueId("reference", "1", "A", "7", None, "ALA"),
        (1, 2.0, 3.0),
        1.8,
        " c ",
    )

    assert atom.atom_id == "atom-1"
    assert atom.coordinate == (1.0, 2.0, 3.0)
    assert atom.element == "C"
    assert atom.has_known_radius is True
    assert atom.to_json()["radius_angstrom"] == 1.8
    assert atom.to_json()["residue_id"] == {
        "structure_id": "reference",
        "model_id": "1",
        "chain_id": "A",
        "auth_seq_id": "7",
        "insertion_code": None,
        "residue_name": "ALA",
    }
    with pytest.raises(FrozenInstanceError):
        atom.atom_id = "changed"  # type: ignore[misc]


@pytest.mark.parametrize(
    ("kwargs", "error", "message"),
    (
        ({"atom_id": " "}, ValueError, "atom_id"),
        ({"coordinate": (1.0, 2.0)}, ValueError, "three finite"),
        ({"coordinate": (1.0, 2.0, math.inf)}, ValueError, "three finite"),
        ({"radius_angstrom": 0.0}, ValueError, "finite and positive"),
    ),
)
def test_pocket_atom_rejects_malformed_values(
    kwargs: dict[str, object],
    error: type[Exception],
    message: str,
) -> None:
    defaults: dict[str, object] = {
        "atom_id": "atom-1",
        "residue_id": ResidueId("reference", "1", "A", "1", None, "ALA"),
        "coordinate": (1.0, 2.0, 3.0),
        "radius_angstrom": 1.8,
    }
    defaults.update(kwargs)

    with pytest.raises(error, match=message):
        PocketAtom(**defaults)  # type: ignore[arg-type]


def test_pocket_atom_marks_an_element_without_a_validated_radius_as_unknown() -> None:
    atom = PocketAtom(
        "unknown",
        ResidueId("reference", "1", "A", "1", None, "ALA"),
        (0.0, 0.0, 0.0),
        1.8,
        "XX",
    )

    assert atom.has_known_radius is False
    result = detect_alpha_spheres(
        _regular_tetrahedron() + (atom,),
        PocketDetectionSettings(
            geometry=PocketGeometrySettings(
                minimum_alpha_sphere_radius_angstrom=2.0,
                maximum_alpha_sphere_radius_angstrom=3.0,
            )
        ),
    )
    assert result.spheres
    assert [item.code for item in result.diagnostics] == ["pocket.detect.unknown_radius"]


def test_alpha_sphere_from_simplex_rejects_malformed_simplex_inputs() -> None:
    simplex = _regular_tetrahedron()

    with pytest.raises(ValueError, match="exactly four"):
        alpha_sphere_from_simplex(simplex[:3], simplex, PocketGeometrySettings())
    with pytest.raises(TypeError, match="PocketAtom"):
        alpha_sphere_from_simplex((object(),) * 4, simplex, PocketGeometrySettings())  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="distinct"):
        alpha_sphere_from_simplex((simplex[0], simplex[0], simplex[2], simplex[3]), simplex, PocketGeometrySettings())
    with pytest.raises(TypeError, match="PocketGeometrySettings"):
        alpha_sphere_from_simplex(simplex, simplex, object())  # type: ignore[arg-type]


def test_detect_alpha_spheres_rejects_untyped_atom_iterables() -> None:
    with pytest.raises(TypeError, match="PocketAtom"):
        detect_alpha_spheres([object()])  # type: ignore[list-item]


def test_detect_alpha_spheres_enforces_atom_cap_before_tessellation() -> None:
    atoms = _regular_tetrahedron() + (_atom("extra", (20.0, 20.0, 20.0), "5"),)

    result = detect_alpha_spheres(atoms, PocketDetectionSettings(max_atom_count=4))

    assert result.spheres == ()
    assert [item.code for item in result.diagnostics] == ["pocket.detect.resource_limit"]


def test_detect_alpha_spheres_reports_duplicate_atom_identity() -> None:
    duplicate = _atom("atom-1", (20.0, 20.0, 20.0), "5")

    result = detect_alpha_spheres(_regular_tetrahedron() + (duplicate,))

    assert [item.code for item in result.diagnostics] == ["pocket.detect.duplicate_atom_id"]
    assert result.availability is Availability.INVALID_INPUT
    assert result.to_json()["availability"] == "invalid_input"


def test_detect_alpha_spheres_honors_cancellation_during_generation() -> None:
    calls = 0

    def cancel_after_tessellation() -> bool:
        nonlocal calls
        calls += 1
        return calls >= 2

    result = detect_alpha_spheres(
        _regular_tetrahedron(),
        PocketDetectionSettings(
            geometry=PocketGeometrySettings(
                minimum_alpha_sphere_radius_angstrom=3.0,
                maximum_alpha_sphere_radius_angstrom=5.0,
            )
        ),
        cancel_check=cancel_after_tessellation,
    )

    assert result.spheres == ()
    assert [item.code for item in result.diagnostics] == ["pocket.detect.cancelled"]
    assert result.availability is Availability.NOT_APPLICABLE
    assert calls == 2


@pytest.mark.parametrize("failure", (ValueError("synthetic value failure"), MemoryError()))
def test_detect_alpha_spheres_contains_non_qhull_tessellation_failures(
    monkeypatch: pytest.MonkeyPatch,
    failure: Exception,
) -> None:
    def raising_delaunay(*args: object, **kwargs: object) -> object:
        raise failure

    monkeypatch.setattr("structlens.core.pockets.delaunay.Delaunay", raising_delaunay)
    result = detect_alpha_spheres(_regular_tetrahedron(), PocketDetectionSettings())

    assert result.spheres == ()
    assert [item.code for item in result.diagnostics] == ["pocket.detect.qhull_failure"]


def test_alpha_sphere_rejects_a_large_nonvertex_surface_intrusion() -> None:
    simplex = _regular_tetrahedron()
    intruder = _atom("large-intruder", (4.1, 0.0, 0.0), "99", radius=2.0)
    settings = PocketGeometrySettings(
        minimum_alpha_sphere_radius_angstrom=3.0,
        maximum_alpha_sphere_radius_angstrom=5.0,
    )

    assert alpha_sphere_from_simplex(simplex, simplex + (intruder,), settings) is None


def test_detect_alpha_spheres_is_repeatable_under_reordered_atoms() -> None:
    atoms = _dense_shell()
    settings = PocketDetectionSettings(
        geometry=PocketGeometrySettings(
            minimum_alpha_sphere_radius_angstrom=2.0,
            maximum_alpha_sphere_radius_angstrom=6.2,
        ),
    )

    first = detect_alpha_spheres(atoms, settings)
    second = detect_alpha_spheres(tuple(reversed(atoms)), settings)

    assert first.diagnostics == ()
    assert second.diagnostics == ()
    assert first.spheres
    assert tuple(item.sphere_id for item in first.spheres) == tuple(item.sphere_id for item in second.spheres)


def test_detect_alpha_spheres_is_stable_under_rigid_transform_and_small_jitter() -> None:
    atoms = _regular_tetrahedron()
    transformed = tuple(
        _atom(
            atom.atom_id,
            (-atom.coordinate[1] + 8.0, atom.coordinate[0] - 3.0, atom.coordinate[2] + 1.5),
            atom.residue_id.auth_seq_id,
        )
        for atom in atoms
    )
    jittered = tuple(
        _atom(
            atom.atom_id,
            (atom.coordinate[0] + 1.0e-7 * index, atom.coordinate[1], atom.coordinate[2]),
            atom.residue_id.auth_seq_id,
        )
        for index, atom in enumerate(atoms)
    )
    settings = PocketDetectionSettings(
        geometry=PocketGeometrySettings(
            minimum_alpha_sphere_radius_angstrom=2.0,
            maximum_alpha_sphere_radius_angstrom=3.0,
        )
    )

    baseline = detect_alpha_spheres(atoms, settings)
    rigid = detect_alpha_spheres(transformed, settings)
    jitter = detect_alpha_spheres(jittered, settings)

    assert len(baseline.spheres) == len(rigid.spheres) == len(jitter.spheres) == 1
    assert rigid.spheres[0].radius_angstrom == pytest.approx(baseline.spheres[0].radius_angstrom)
    assert jitter.spheres[0].radius_angstrom == pytest.approx(
        baseline.spheres[0].radius_angstrom,
        abs=1.0e-6,
    )


def test_detect_alpha_spheres_reports_insufficient_atoms_without_placeholder_spheres() -> None:
    result = detect_alpha_spheres(_regular_tetrahedron()[:3], PocketDetectionSettings())

    assert result.spheres == ()
    assert [item.code for item in result.diagnostics] == ["pocket.detect.insufficient_atoms"]
    assert result.availability is Availability.INVALID_INPUT


def test_detect_alpha_spheres_enforces_resource_limit_before_tessellation() -> None:
    result = detect_alpha_spheres(
        _dense_shell(),
        PocketDetectionSettings(max_estimated_simplices=10),
    )

    assert result.spheres == ()
    assert [item.code for item in result.diagnostics] == ["pocket.detect.resource_limit"]


def test_detect_alpha_spheres_enforces_indexed_clearance_work_limit() -> None:
    result = detect_alpha_spheres(
        _regular_tetrahedron(),
        PocketDetectionSettings(
            geometry=PocketGeometrySettings(
                minimum_alpha_sphere_radius_angstrom=2.0,
                maximum_alpha_sphere_radius_angstrom=3.0,
            ),
            max_clearance_atom_checks=3,
        ),
    )

    assert result.spheres == ()
    assert [item.code for item in result.diagnostics] == ["pocket.detect.resource_limit"]


def test_simplex_budget_uses_the_three_dimensional_delaunay_upper_bound(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    atoms = tuple(
        _atom(
            f"atom-{index}",
            (
                10.0 * math.cos(index),
                10.0 * math.sin(index),
                float(index % 7),
            ),
            str(index),
        )
        for index in range(100)
    )
    called = False

    class _EmptyTessellation:
        simplices = np.empty((0, 4), dtype=int)

    def recording_delaunay(coordinates: np.ndarray) -> _EmptyTessellation:
        nonlocal called
        called = coordinates.shape == (100, 3)
        return _EmptyTessellation()

    monkeypatch.setattr("structlens.core.pockets.delaunay.Delaunay", recording_delaunay)

    result = detect_alpha_spheres(
        atoms,
        PocketDetectionSettings(max_estimated_simplices=5000),
    )

    assert called is True
    assert result.diagnostics == ()


def test_detect_alpha_spheres_returns_typed_qhull_failure(monkeypatch: pytest.MonkeyPatch) -> None:
    atoms = _dense_shell()

    def raising_delaunay(*args: object, **kwargs: object) -> object:
        raise QhullError("synthetic qhull failure")

    monkeypatch.setattr("structlens.core.pockets.delaunay.Delaunay", raising_delaunay)
    result = detect_alpha_spheres(atoms, PocketDetectionSettings())

    assert result.spheres == ()
    assert [item.code for item in result.diagnostics] == ["pocket.detect.qhull_failure"]
    assert all(isinstance(item, Diagnostic) for item in result.diagnostics)
