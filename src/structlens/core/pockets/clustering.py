"""Deterministic clustering and solvent-connectivity checks for pockets.

The detector produces individual alpha spheres.  This module turns those
spheres into candidate pockets using a small, reproducible graph operation.
It deliberately does not attach a trained ``ligandability`` probability to a
candidate: geometric support and solvent exposure remain separate evidence.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import cast

import numpy as np

from structlens.core.evidence import Diagnostic, DiagnosticSeverity

from .delaunay import PocketAtom
from .models import AlphaSphere, PocketCandidate, PocketDetectionSettings


@dataclass(frozen=True, slots=True)
class PocketClusteringResult:
    """Immutable output of alpha-sphere graph clustering."""

    candidates: tuple[PocketCandidate, ...] = ()
    diagnostics: tuple[Diagnostic, ...] = ()

    def __post_init__(self) -> None:
        candidates = tuple(self.candidates)
        diagnostics = tuple(self.diagnostics)
        if any(not isinstance(item, PocketCandidate) for item in candidates):
            raise TypeError("candidates must contain PocketCandidate values")
        if any(not isinstance(item, Diagnostic) for item in diagnostics):
            raise TypeError("diagnostics must contain Diagnostic values")
        object.__setattr__(self, "candidates", candidates)
        object.__setattr__(self, "diagnostics", diagnostics)


# A descriptive alias makes the result convenient for callers that use the
# shorter terminology without duplicating the contract.
ClusterResult = PocketClusteringResult


def _diagnostic(code: str, message: str, remediation: str) -> Diagnostic:
    return Diagnostic(
        code=code,
        severity=DiagnosticSeverity.ERROR,
        message=message,
        remediation=remediation,
    )


def _union_find(size: int) -> tuple[list[int], list[int]]:
    return list(range(size)), [1] * size


def _find(parent: list[int], index: int) -> int:
    root = index
    while parent[root] != root:
        root = parent[root]
    while parent[index] != index:
        next_index = parent[index]
        parent[index] = root
        index = next_index
    return root


def _union(parent: list[int], rank: list[int], left: int, right: int) -> None:
    left_root = _find(parent, left)
    right_root = _find(parent, right)
    if left_root == right_root:
        return
    if rank[left_root] < rank[right_root]:
        left_root, right_root = right_root, left_root
    parent[right_root] = left_root
    if rank[left_root] == rank[right_root]:
        rank[left_root] += 1


def _cluster_components(
    spheres: tuple[AlphaSphere, ...],
    padding_angstrom: float,
) -> tuple[tuple[AlphaSphere, ...], ...]:
    """Build components using distance-expanded sphere overlap.

    Iteration is over sphere IDs rather than input order.  Thus atom/sphere
    reordering cannot change component membership or the resulting candidate
    identity.
    """

    ordered = tuple(sorted(spheres, key=lambda item: item.sphere_id))
    parent, rank = _union_find(len(ordered))
    for left_index, left in enumerate(ordered):
        left_center = np.asarray(left.center_xyz, dtype=np.float64)
        for right_index in range(left_index + 1, len(ordered)):
            right = ordered[right_index]
            right_center = np.asarray(right.center_xyz, dtype=np.float64)
            distance = float(np.linalg.norm(left_center - right_center))
            threshold = left.radius_angstrom + right.radius_angstrom + padding_angstrom
            if distance <= threshold:
                _union(parent, rank, left_index, right_index)
    components: dict[int, list[AlphaSphere]] = {}
    for index, sphere in enumerate(ordered):
        components.setdefault(_find(parent, index), []).append(sphere)
    return tuple(
        tuple(sorted(component, key=lambda item: item.sphere_id))
        for component in sorted(components.values(), key=lambda item: item[0].sphere_id)
    )


def _validated_atoms(atoms: Sequence[PocketAtom]) -> tuple[PocketAtom, ...]:
    ordered = tuple(sorted(tuple(atoms), key=lambda item: item.atom_id))
    for atom in ordered:
        coordinate = tuple(float(value) for value in atom.coordinate)
        if len(coordinate) != 3 or any(not math.isfinite(value) for value in coordinate):
            raise ValueError("atom coordinates must contain exactly three finite values")
        radius = float(atom.radius_angstrom)
        if not math.isfinite(radius) or radius <= 0.0:
            raise ValueError("atom radii must be finite and positive")
    return ordered


def _setting(
    settings: PocketDetectionSettings,
    name: str,
    default: float | int,
    *,
    non_negative: bool = False,
) -> float | int:
    value = getattr(settings, name, default)
    if isinstance(default, float):
        value = float(value)
        if not math.isfinite(value) or (value < 0.0 if non_negative else value <= 0.0):
            bound = "non-negative" if non_negative else "positive"
            raise ValueError(f"{name} must be finite and {bound}")
    else:
        value = int(value)
        if value <= 0:
            raise ValueError(f"{name} must be positive")
    return value


@dataclass(frozen=True, slots=True)
class _SolventGrid:
    blocked: np.ndarray
    outside: np.ndarray
    lower: np.ndarray
    spacing: float
    shape: tuple[int, int, int]
    raster_work: int


def _cell_bounds(
    center: Sequence[float],
    radius: float,
    lower: np.ndarray,
    spacing: float,
    shape: tuple[int, int, int],
) -> tuple[tuple[int, int], tuple[int, int], tuple[int, int]] | None:
    center_array = np.asarray(center, dtype=np.float64)
    minimum = np.ceil((center_array - radius - lower) / spacing - 0.5).astype(int)
    maximum = np.floor((center_array + radius - lower) / spacing - 0.5).astype(int)
    bounds: list[tuple[int, int]] = []
    for axis in range(3):
        first = max(0, int(minimum[axis]))
        last = min(shape[axis] - 1, int(maximum[axis]))
        if first > last:
            return None
        bounds.append((first, last))
    return (bounds[0], bounds[1], bounds[2])


def _bounds_work(bounds: tuple[tuple[int, int], tuple[int, int], tuple[int, int]]) -> int:
    return math.prod(last - first + 1 for first, last in bounds)


def _local_squared_distances(
    bounds: tuple[tuple[int, int], tuple[int, int], tuple[int, int]],
    center: Sequence[float],
    lower: np.ndarray,
    spacing: float,
) -> np.ndarray:
    local_shape = tuple(last - first + 1 for first, last in bounds)
    indices = np.indices(local_shape, dtype=np.float64)
    squared = np.zeros(local_shape, dtype=np.float64)
    for axis, (first, _) in enumerate(bounds):
        coordinates = lower[axis] + (indices[axis] + first + 0.5) * spacing
        squared += (coordinates - float(center[axis])) ** 2
    return squared


def _slices(
    bounds: tuple[tuple[int, int], tuple[int, int], tuple[int, int]],
) -> tuple[slice, slice, slice]:
    return tuple(slice(first, last + 1) for first, last in bounds)  # type: ignore[return-value]


def _build_solvent_grid(
    atoms: tuple[PocketAtom, ...],
    settings: PocketDetectionSettings,
) -> _SolventGrid:
    """Rasterize inflated atoms once and flood-fill outside solvent."""

    spacing = float(_setting(settings, "solvent_grid_spacing_angstrom", 0.75))
    margin = float(
        _setting(settings, "solvent_boundary_margin_angstrom", 2.0, non_negative=True)
    )
    geometry = getattr(settings, "geometry", None)
    probe_default = float(getattr(geometry, "probe_radius_angstrom", 1.4))
    probe = float(_setting(settings, "probe_radius_angstrom", probe_default))
    grid_cap = int(_setting(settings, "max_solvent_grid_cells", 2_000_000))
    raster_cap = int(_setting(settings, "max_solvent_raster_cells", 20_000_000))

    coordinates = np.asarray([atom.coordinate for atom in atoms], dtype=np.float64)
    radii = np.asarray([atom.radius_angstrom + probe for atom in atoms], dtype=np.float64)
    outside_padding = float(np.max(radii)) + margin
    lower = coordinates.min(axis=0) - outside_padding
    upper = coordinates.max(axis=0) + outside_padding
    shape = cast(
        tuple[int, int, int],
        tuple(max(1, int(math.ceil((upper[index] - lower[index]) / spacing))) for index in range(3)),
    )
    cell_count = math.prod(shape)
    if cell_count > grid_cap:
        raise MemoryError("solvent grid exceeds the configured cell budget")

    atom_bounds = tuple(
        _cell_bounds(atom.coordinate, radius, lower, spacing, shape)
        for atom, radius in zip(atoms, radii, strict=True)
    )
    raster_work = sum(_bounds_work(bounds) for bounds in atom_bounds if bounds is not None)
    if raster_work > raster_cap:
        raise MemoryError("solvent rasterization exceeds the configured work budget")

    blocked = np.zeros(shape, dtype=bool)
    for atom, radius, bounds in zip(atoms, radii, atom_bounds, strict=True):
        if bounds is None:
            continue
        view = blocked[_slices(bounds)]
        view |= _local_squared_distances(bounds, atom.coordinate, lower, spacing) <= radius * radius

    outside = np.zeros(shape, dtype=bool)
    stack: list[tuple[int, int, int]] = []
    for index in np.ndindex(shape):
        if any(index[axis] in (0, shape[axis] - 1) for axis in range(3)) and not blocked[index]:
            outside[index] = True
            stack.append(index)
    while stack:
        current = stack.pop()
        for axis in range(3):
            for step in (-1, 1):
                neighbour = list(current)
                neighbour[axis] += step
                if not 0 <= neighbour[axis] < shape[axis]:
                    continue
                item = (neighbour[0], neighbour[1], neighbour[2])
                if not blocked[item] and not outside[item]:
                    outside[item] = True
                    stack.append(item)

    return _SolventGrid(blocked, outside, lower, spacing, shape, raster_work)


def _cluster_is_exposed(
    component: tuple[AlphaSphere, ...],
    grid: _SolventGrid,
    remaining_work: int,
) -> tuple[bool, int]:
    work = 0
    for sphere in component:
        bounds = _cell_bounds(
            sphere.center_xyz,
            sphere.radius_angstrom,
            grid.lower,
            grid.spacing,
            grid.shape,
        )
        if bounds is None:
            continue
        work += _bounds_work(bounds)
        if work > remaining_work:
            raise MemoryError("candidate rasterization exceeds the configured work budget")
        region = _slices(bounds)
        inside = (
            _local_squared_distances(
                bounds,
                sphere.center_xyz,
                grid.lower,
                grid.spacing,
            )
            <= sphere.radius_angstrom * sphere.radius_angstrom
        )
        if bool(np.any(inside & ~grid.blocked[region] & grid.outside[region])):
            return True, work
    return False, work


def cluster_alpha_spheres(
    spheres: Sequence[AlphaSphere],
    atoms: Sequence[PocketAtom],
    settings: PocketDetectionSettings,
) -> PocketClusteringResult:
    """Cluster spheres and retain minimum-size, buried components.

    The input is never mutated.  Candidates and diagnostics are sorted by
    stable content-derived IDs before being returned.
    """

    if not isinstance(settings, PocketDetectionSettings):
        raise TypeError("settings must be PocketDetectionSettings")
    sphere_values = tuple(spheres)
    if any(not isinstance(item, AlphaSphere) for item in sphere_values):
        raise TypeError("spheres must contain AlphaSphere values")
    if not sphere_values:
        return PocketClusteringResult()
    atom_values = _validated_atoms(atoms)
    padding = float(
        _setting(settings, "cluster_distance_padding_angstrom", 0.0, non_negative=True)
    )
    minimum_size = int(_setting(settings, "minimum_cluster_size", 1))
    components = _cluster_components(sphere_values, padding)
    eligible_components = tuple(component for component in components if len(component) >= minimum_size)
    candidates: list[PocketCandidate] = []
    diagnostics: list[Diagnostic] = []
    solvent_grid: _SolventGrid | None = None
    raster_work = 0
    if atom_values and eligible_components:
        try:
            solvent_grid = _build_solvent_grid(atom_values, settings)
            raster_work = solvent_grid.raster_work
        except MemoryError:
            return PocketClusteringResult(
                diagnostics=(
                    _diagnostic(
                        "pocket.detect.resource_limit",
                        "The solvent-connectivity grid exceeds the configured resource limit.",
                        "Increase the grid/work budget or use a coarser solvent-grid spacing.",
                    ),
                )
            )
    raster_cap = int(_setting(settings, "max_solvent_raster_cells", 20_000_000))
    for component in eligible_components:
        exposed = False
        if solvent_grid is not None:
            try:
                exposed, component_work = _cluster_is_exposed(
                    component,
                    solvent_grid,
                    raster_cap - raster_work,
                )
                raster_work += component_work
            except MemoryError:
                return PocketClusteringResult(
                    diagnostics=(
                        _diagnostic(
                            "pocket.detect.resource_limit",
                            "Pocket-region rasterization exceeds the configured resource limit.",
                            "Increase the grid/work budget or use a coarser solvent-grid spacing.",
                        ),
                    )
                )
        if exposed:
            diagnostics.append(
                _diagnostic(
                    "pocket.detect.solvent_exposed_cluster",
                    "The alpha-sphere cluster is connected to outside solvent.",
                    "Inspect the cavity opening or tighten the declared detector settings.",
                )
            )
            continue
        candidates.append(PocketCandidate(alpha_spheres=component))
    # Spatial order is useful to callers presenting an unranked detector
    # result; the content-derived ID remains the final tie-breaker.
    candidates.sort(key=lambda item: (*item.centroid_xyz, item.candidate_id))
    diagnostics.sort(key=lambda item: (item.code, item.message))
    return PocketClusteringResult(tuple(candidates), tuple(diagnostics))


__all__ = ["ClusterResult", "PocketClusteringResult", "cluster_alpha_spheres"]
