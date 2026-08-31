"""Analytic geometry helpers for pocket detection."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Final

import numpy as np

from structlens.core.evidence import Diagnostic, DiagnosticSeverity

_DEGENERATE_CONDITION_LIMIT: Final[float] = 1.0e12


def _diagnostic(code: str, message: str) -> Diagnostic:
    return Diagnostic(
        code=code,
        severity=DiagnosticSeverity.ERROR,
        message=message,
        remediation="Review the simplex coordinates and exclude degenerate or non-finite atom sets.",
    )


@dataclass(frozen=True, slots=True)
class Circumsphere:
    center_xyz: tuple[float, float, float]
    radius_angstrom: float

    def __post_init__(self) -> None:
        try:
            center = tuple(float(value) for value in self.center_xyz)
        except (TypeError, ValueError) as exc:
            raise ValueError("center_xyz must contain numeric values") from exc
        if len(center) != 3 or any(not math.isfinite(value) for value in center):
            raise ValueError("center_xyz must contain exactly three finite values")
        radius = float(self.radius_angstrom)
        if not math.isfinite(radius) or radius <= 0.0:
            raise ValueError("radius_angstrom must be finite and positive")
        object.__setattr__(self, "center_xyz", center)
        object.__setattr__(self, "radius_angstrom", radius)


@dataclass(frozen=True, slots=True)
class CircumsphereResult:
    sphere: Circumsphere | None = None
    diagnostics: tuple[Diagnostic, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        if self.sphere is not None and not isinstance(self.sphere, Circumsphere):
            raise TypeError("sphere must be a Circumsphere or None")
        diagnostics = tuple(self.diagnostics)
        if any(not isinstance(item, Diagnostic) for item in diagnostics):
            raise TypeError("diagnostics must contain Diagnostic values")
        object.__setattr__(self, "diagnostics", diagnostics)
        if self.sphere is not None and diagnostics:
            raise ValueError("successful circumsphere results must not carry diagnostics")
        if self.sphere is None and not diagnostics:
            raise ValueError("failed circumsphere results must carry a diagnostic")


def tetrahedron_circumsphere(
    points: tuple[tuple[float, float, float], ...] | list[tuple[float, float, float]],
    *,
    tolerance: float = 1.0e-12,
) -> CircumsphereResult:
    """Return the unique circumsphere for four non-coplanar points."""

    if not math.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("tolerance must be finite and positive")
    try:
        array = np.asarray(points, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ValueError("points must contain four 3D numeric vertices") from exc
    if array.shape != (4, 3):
        raise ValueError("points must contain four 3D vertices")
    if len({tuple(float(value) for value in row) for row in array.tolist()}) != 4:
        raise ValueError("points must contain four distinct 3D vertices")
    if not np.isfinite(array).all():
        return CircumsphereResult(
            diagnostics=(
                _diagnostic(
                    "pocket.geometry.nonfinite_coordinate",
                    "Pocket simplex coordinates must be numeric and finite.",
                ),
            )
        )
    origin = array[0]
    system = 2.0 * (array[1:] - origin)
    rhs = np.sum(array[1:] * array[1:], axis=1) - np.sum(origin * origin)
    try:
        singular_values = np.linalg.svd(system, compute_uv=False)
        condition = float(np.linalg.cond(system))
    except np.linalg.LinAlgError:
        singular_values = np.asarray((0.0, 0.0, 0.0), dtype=np.float64)
        condition = math.inf
    if (
        not math.isfinite(condition)
        or condition >= min(_DEGENERATE_CONDITION_LIMIT, 1.0 / tolerance)
        or float(singular_values[-1]) <= tolerance
    ):
        return CircumsphereResult(
            diagnostics=(
                _diagnostic(
                    "pocket.geometry.degenerate_simplex",
                    "Pocket simplex points are degenerate or numerically ill-conditioned.",
                ),
            )
        )
    try:
        center = np.linalg.solve(system, rhs)
    except np.linalg.LinAlgError:
        return CircumsphereResult(
            diagnostics=(
                _diagnostic(
                    "pocket.geometry.degenerate_simplex",
                    "Pocket simplex points are degenerate or numerically ill-conditioned.",
                ),
            )
        )
    radius = float(np.linalg.norm(center - origin))
    if not np.isfinite(center).all() or not math.isfinite(radius) or radius <= 0.0:
        return CircumsphereResult(
            diagnostics=(
                _diagnostic(
                    "pocket.geometry.degenerate_simplex",
                    "Pocket simplex circumsphere could not be solved as a finite sphere.",
                ),
            )
        )
    return CircumsphereResult(
        sphere=Circumsphere(
            center_xyz=(float(center[0]), float(center[1]), float(center[2])),
            radius_angstrom=radius,
        )
    )


__all__ = [
    "Circumsphere",
    "CircumsphereResult",
    "tetrahedron_circumsphere",
]
