from __future__ import annotations

import math

import numpy as np
import pytest

from structlens.core.evidence import Diagnostic
from structlens.core.pockets import (
    CircumsphereResult,
    tetrahedron_circumsphere,
)


def _rotate_z(point: tuple[float, float, float]) -> tuple[float, float, float]:
    return (-point[1], point[0], point[2])


def _translate(point: tuple[float, float, float], offset: tuple[float, float, float]) -> tuple[float, float, float]:
    return tuple(value + delta for value, delta in zip(point, offset, strict=True))


def test_tetrahedron_circumsphere_matches_known_regular_solution() -> None:
    points = (
        (1.0, 1.0, 1.0),
        (-1.0, -1.0, 1.0),
        (-1.0, 1.0, -1.0),
        (1.0, -1.0, -1.0),
    )

    result = tetrahedron_circumsphere(points)

    assert isinstance(result, CircumsphereResult)
    assert result.diagnostics == ()
    assert result.sphere is not None
    assert result.sphere.center_xyz == pytest.approx((0.0, 0.0, 0.0))
    assert result.sphere.radius_angstrom == pytest.approx(math.sqrt(3.0))


def test_tetrahedron_circumsphere_is_rigid_body_invariant() -> None:
    points = (
        (4.0, 3.0, 2.0),
        (2.0, 1.0, 2.0),
        (2.0, 3.0, 0.0),
        (4.0, 1.0, 0.0),
    )
    offset = (11.5, -7.25, 3.0)
    rotated = tuple(_translate(_rotate_z(point), offset) for point in points)

    first = tetrahedron_circumsphere(points)
    second = tetrahedron_circumsphere(rotated)

    assert first.sphere is not None
    assert second.sphere is not None
    assert first.diagnostics == ()
    assert second.diagnostics == ()
    assert second.sphere.radius_angstrom == pytest.approx(first.sphere.radius_angstrom)
    assert second.sphere.center_xyz == pytest.approx(_translate(_rotate_z(first.sphere.center_xyz), offset))


def test_tetrahedron_circumsphere_is_invariant_to_vertex_permutation() -> None:
    points = (
        (1.0, 1.0, 1.0),
        (-1.0, -1.0, 1.0),
        (-1.0, 1.0, -1.0),
        (1.0, -1.0, -1.0),
    )

    first = tetrahedron_circumsphere(points)
    permuted = tetrahedron_circumsphere((points[2], points[0], points[3], points[1]))

    assert first.sphere is not None
    assert permuted.sphere is not None
    assert permuted.sphere.center_xyz == pytest.approx(first.sphere.center_xyz, abs=1e-12)
    assert permuted.sphere.radius_angstrom == pytest.approx(first.sphere.radius_angstrom, abs=1e-12)


@pytest.mark.parametrize(
    "points, code",
    (
        (
            (
                (0.0, 0.0, 0.0),
                (1.0, 0.0, 0.0),
                (0.0, 1.0, 0.0),
                (1.0, 1.0, 0.0),
            ),
            "pocket.geometry.degenerate_simplex",
        ),
        (
            (
                (0.0, 0.0, 0.0),
                (1.0, 0.0, 0.0),
                (0.0, 1.0, 0.0),
                (float("nan"), 1.0, 1.0),
            ),
            "pocket.geometry.nonfinite_coordinate",
        ),
    ),
)
def test_tetrahedron_circumsphere_returns_typed_diagnostics_for_invalid_input(
    points: tuple[tuple[float, float, float], ...],
    code: str,
) -> None:
    result = tetrahedron_circumsphere(points)

    assert isinstance(result, CircumsphereResult)
    assert result.sphere is None
    assert [diagnostic.code for diagnostic in result.diagnostics] == [code]
    assert all(isinstance(diagnostic, Diagnostic) for diagnostic in result.diagnostics)
    assert all(diagnostic.severity.value == "error" for diagnostic in result.diagnostics)


def test_nearly_degenerate_tetrahedron_respects_condition_tolerance() -> None:
    points = (
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1e-10),
    )

    rejected = tetrahedron_circumsphere(points, tolerance=1e-8)
    accepted = tetrahedron_circumsphere(points, tolerance=1e-14)

    assert rejected.sphere is None
    assert [diagnostic.code for diagnostic in rejected.diagnostics] == ["pocket.geometry.degenerate_simplex"]
    assert accepted.sphere is not None
    assert accepted.diagnostics == ()


def test_tetrahedron_circumsphere_never_returns_nan_output() -> None:
    result = tetrahedron_circumsphere(
        (
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (1.0, 1.0, 0.0),
        )
    )

    assert result.sphere is None
    assert result.diagnostics
    assert all(isinstance(diagnostic, Diagnostic) for diagnostic in result.diagnostics)
    assert result.sphere is None or (
        np.isfinite(result.sphere.center_xyz).all() and np.isfinite(result.sphere.radius_angstrom)
    )


@pytest.mark.parametrize(
    "points",
    (
        ((0.0, 0.0, 0.0),) * 4,
        ((0.0, 0.0),) * 4,
        ((0.0, 0.0, 0.0),) * 3,
    ),
)
def test_tetrahedron_circumsphere_rejects_malformed_shapes(
    points: tuple[tuple[float, ...], ...],
) -> None:
    with pytest.raises(ValueError, match="four|shape|three"):
        tetrahedron_circumsphere(points)


def test_tetrahedron_circumsphere_rejects_invalid_tolerance() -> None:
    with pytest.raises(ValueError, match="tolerance"):
        tetrahedron_circumsphere(
            (
                (1.0, 1.0, 1.0),
                (-1.0, -1.0, 1.0),
                (-1.0, 1.0, -1.0),
                (1.0, -1.0, -1.0),
            ),
            tolerance=0.0,
        )


def test_tetrahedron_circumsphere_rejects_nonnumeric_vertices() -> None:
    with pytest.raises(ValueError, match="numeric"):
        tetrahedron_circumsphere(
            (
                ("bad", 0.0, 0.0),  # type: ignore[arg-type]
                (1.0, 0.0, 0.0),
                (0.0, 1.0, 0.0),
                (0.0, 0.0, 1.0),
            )
        )


def test_circumsphere_result_validates_internal_state() -> None:
    with pytest.raises(TypeError, match="Circumsphere"):
        CircumsphereResult(sphere=object())  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="Diagnostic"):
        CircumsphereResult(diagnostics=(object(),))  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="must not carry diagnostics"):
        CircumsphereResult(
            sphere=tetrahedron_circumsphere(
                (
                    (1.0, 1.0, 1.0),
                    (-1.0, -1.0, 1.0),
                    (-1.0, 1.0, -1.0),
                    (1.0, -1.0, -1.0),
                )
            ).sphere,
            diagnostics=(
                Diagnostic(
                    code="pocket.geometry.synthetic",
                    severity="error",
                    message="synthetic",
                ),
            ),
        )
    with pytest.raises(ValueError, match="must carry a diagnostic"):
        CircumsphereResult()


def test_tetrahedron_circumsphere_handles_linear_algebra_failure_paths(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    points = (
        (1.0, 1.0, 1.0),
        (-1.0, -1.0, 1.0),
        (-1.0, 1.0, -1.0),
        (1.0, -1.0, -1.0),
    )

    def raising_svd(*args: object, **kwargs: object) -> object:
        raise np.linalg.LinAlgError("boom")

    monkeypatch.setattr(np.linalg, "svd", raising_svd)
    result = tetrahedron_circumsphere(points)
    assert result.sphere is None
    assert [diagnostic.code for diagnostic in result.diagnostics] == ["pocket.geometry.degenerate_simplex"]


def test_tetrahedron_circumsphere_handles_nonfinite_solve_output(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    points = (
        (1.0, 1.0, 1.0),
        (-1.0, -1.0, 1.0),
        (-1.0, 1.0, -1.0),
        (1.0, -1.0, -1.0),
    )

    monkeypatch.setattr(np.linalg, "solve", lambda *args, **kwargs: np.asarray((math.inf, 0.0, 0.0)))
    result = tetrahedron_circumsphere(points)
    assert result.sphere is None
    assert [diagnostic.code for diagnostic in result.diagnostics] == ["pocket.geometry.degenerate_simplex"]
