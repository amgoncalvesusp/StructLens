"""Bounded, deterministic free-volume measurements for pocket candidates.

The volume reported here is the union of accepted alpha-sphere interiors after
optional atom exclusions.  It is deliberately a grid measurement, not a
continuous cavity-volume estimate.  All mutable numerical work is kept local
to a chunk; the implementation never allocates a full three-dimensional grid.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.spatial import cKDTree  # type: ignore[import-untyped]

from structlens.core.evidence import Availability, Diagnostic, DiagnosticSeverity
from structlens.core.models import AtomRecord, ComponentKind, StructureComponent

from .models import PocketCandidate
from .radii import vdw_radius_angstrom
from .volume_models import (
    _GRID_PHASE,
    _GRID_SAMPLING,
    _VOLUME_UNITS,
    PocketVolumeComparison,
    PocketVolumeResult,
    PocketVolumeSensitivity,
    PocketVolumeSensitivityResult,
    PocketVolumeSettings,
    _diagnostic,
    _finite,
)


@dataclass(frozen=True, slots=True)
class _ExclusionAtom:
    coordinate: tuple[float, float, float]
    radius_angstrom: float


@dataclass(frozen=True, slots=True)
class _GridMeasurement:
    count: int
    shape: tuple[int, int, int]
    envelope_volume_angstrom3: float
    exclusion_neighbor_checks: int = 0


class _VolumeResourceLimit(RuntimeError):
    """Internal bounded-work stop converted to a typed public result."""


def sensitivity_from_volumes(
    coarse_volume_angstrom3: float | None,
    fine_volume_angstrom3: float | None,
) -> PocketVolumeSensitivity:
    """Calculate absolute and fine-relative grid sensitivity without NaN."""

    if coarse_volume_angstrom3 is None or fine_volume_angstrom3 is None:
        return PocketVolumeSensitivity()
    coarse = _finite(coarse_volume_angstrom3, "coarse_volume_angstrom3", non_negative=True)
    fine = _finite(fine_volume_angstrom3, "fine_volume_angstrom3", non_negative=True)
    absolute = abs(fine - coarse)
    relative = absolute / fine if fine > 0.0 else None
    return PocketVolumeSensitivity(absolute_angstrom3=absolute, relative_fraction=relative)


def measure_pocket_volume(
    candidate: PocketCandidate | None,
    *,
    protein_atoms: Sequence[AtomRecord] = (),
    retained_components: Sequence[StructureComponent] = (),
    settings: PocketVolumeSettings | None = None,
) -> PocketVolumeResult:
    """Measure the bounded union of candidate spheres on coarse and fine grids."""

    run_settings = settings if settings is not None else PocketVolumeSettings()
    if not isinstance(run_settings, PocketVolumeSettings):
        raise TypeError("settings must be PocketVolumeSettings or None")
    protein = _normalize_atoms(protein_atoms, "protein_atoms")
    components = _normalize_components(retained_components)
    if candidate is None:
        return PocketVolumeResult(
            availability=Availability.NOT_APPLICABLE,
            diagnostics=(
                _diagnostic(
                    "pocket.volume.no_candidate",
                    "No pocket candidate was supplied for free-volume measurement.",
                    remediation="Select a detected pocket or define a focused pocket region before measuring volume.",
                ),
            ),
            units=_VOLUME_UNITS,
            provenance_parameters={
                "radii_version": run_settings.radii_version,
                "component_exclusion_policy": run_settings.component_exclusion_policy,
                "grid_phase": _GRID_PHASE,
                "grid_sampling": _GRID_SAMPLING,
            },
            settings=run_settings,
        )
    if not isinstance(candidate, PocketCandidate):
        raise TypeError("candidate must be PocketCandidate or None")

    spheres = tuple(sorted(candidate.alpha_spheres, key=lambda item: item.sphere_id))
    diagnostics: list[Diagnostic] = []
    if len(spheres) > run_settings.max_sphere_count:
        diagnostics.append(
            _resource_diagnostic(
                "The pocket candidate exceeds the configured alpha-sphere limit.",
                "Use a bounded detected candidate or raise max_sphere_count deliberately.",
            )
        )
        return _resource_result_without_grid(candidate, run_settings, diagnostics)
    lower, upper = _bounds(spheres, run_settings.boundary_margin_angstrom)
    origin = (float(lower[0]), float(lower[1]), float(lower[2]))
    selected_components = tuple(
        component
        for component in components
        if component.kind in {ComponentKind.LIGAND, ComponentKind.ION, ComponentKind.OTHER}
        and component.metadata.get("selected_for_analysis") is True
    )
    exclusion_input_count = len(protein)
    if run_settings.component_exclusion_policy == "unoccupied":
        exclusion_input_count += sum(len(component.atoms) for component in selected_components)
    if exclusion_input_count > run_settings.max_exclusion_atom_count:
        diagnostics.append(
            _resource_diagnostic(
                "The selected structure exceeds the configured exclusion-atom limit.",
                "Narrow the selected structure or raise max_exclusion_atom_count deliberately.",
            )
        )
        return _resource_result_without_grid(candidate, run_settings, diagnostics)
    exclusion_atoms = list(_known_exclusion_atoms(protein, diagnostics, scope="protein"))
    if run_settings.component_exclusion_policy == "unoccupied":
        for component in selected_components:
            exclusion_atoms.extend(
                _known_exclusion_atoms(
                    component.atoms,
                    diagnostics,
                    scope=component.component_id,
                )
            )

    coarse_shape = _shape_for_bounds(lower, upper, run_settings.coarse_grid_spacing_angstrom)
    fine_shape = _shape_for_bounds(lower, upper, run_settings.fine_grid_spacing_angstrom)
    if any(item.code == "pocket.volume.unknown_radius" for item in diagnostics):
        return _result(
            candidate,
            run_settings,
            Availability.INVALID_INPUT,
            origin,
            coarse_shape,
            fine_shape,
            diagnostics,
            None,
            None,
            None,
            None,
            None,
            (),
        )
    active_exclusions = _active_exclusion_atoms(exclusion_atoms, lower, upper)
    tree, radii = _build_exclusion_index(active_exclusions)
    total_coarse = math.prod(coarse_shape)
    total_fine = math.prod(fine_shape)
    sphere_voxel_checks = (total_coarse + total_fine) * len(spheres)
    if (
        total_coarse > run_settings.max_voxel_count
        or total_fine > run_settings.max_voxel_count
        or sphere_voxel_checks > run_settings.max_sphere_voxel_checks
    ):
        diagnostics.append(
            _resource_diagnostic(
                "The requested pocket volume grid exceeds a configured voxel/work limit.",
                "Increase grid spacing, reduce the boundary margin, or raise the relevant work limit deliberately.",
            )
        )
        return _result(
            candidate,
            run_settings,
            Availability.NUMERICAL_FAILURE,
            origin,
            coarse_shape,
            fine_shape,
            diagnostics,
            None,
            None,
            None,
            None,
            None,
            active_exclusions,
        )

    try:
        coarse = _measure_grid(
            origin,
            coarse_shape,
            run_settings.coarse_grid_spacing_angstrom,
            spheres,
            tree,
            radii,
            run_settings.voxel_chunk_size,
            run_settings.max_exclusion_neighbor_checks,
        )
        fine = _measure_grid(
            origin,
            fine_shape,
            run_settings.fine_grid_spacing_angstrom,
            spheres,
            tree,
            radii,
            run_settings.voxel_chunk_size,
            run_settings.max_exclusion_neighbor_checks - coarse.exclusion_neighbor_checks,
        )
    except _VolumeResourceLimit:
        diagnostics.append(
            _resource_diagnostic(
                "Pocket volume exclusion queries exceeded the configured neighbor-check limit.",
                "Narrow the structure or raise max_exclusion_neighbor_checks deliberately.",
            )
        )
        return _result(
            candidate,
            run_settings,
            Availability.NUMERICAL_FAILURE,
            origin,
            coarse_shape,
            fine_shape,
            diagnostics,
            None,
            None,
            None,
            None,
            None,
            active_exclusions,
        )
    coarse_volume = coarse.count * run_settings.coarse_grid_spacing_angstrom**3
    fine_volume = fine.count * run_settings.fine_grid_spacing_angstrom**3
    sensitivity = sensitivity_from_volumes(coarse_volume, fine_volume)
    result = _result(
        candidate,
        run_settings,
        Availability.AVAILABLE,
        origin,
        coarse_shape,
        fine_shape,
        diagnostics,
        coarse.count,
        fine.count,
        coarse_volume,
        fine_volume,
        sensitivity,
        active_exclusions,
    )
    return result


def _normalize_atoms(atoms: Sequence[AtomRecord], name: str) -> tuple[AtomRecord, ...]:
    values = tuple(atoms)
    if any(not isinstance(atom, AtomRecord) for atom in values):
        raise TypeError(f"{name} must contain AtomRecord values")
    return values


def _normalize_components(components: Sequence[StructureComponent]) -> tuple[StructureComponent, ...]:
    values = tuple(components)
    if any(not isinstance(component, StructureComponent) for component in values):
        raise TypeError("retained_components must contain StructureComponent values")
    return values


def _known_exclusion_atoms(
    atoms: Sequence[AtomRecord],
    diagnostics: list[Diagnostic],
    *,
    scope: str,
) -> tuple[_ExclusionAtom, ...]:
    known: list[_ExclusionAtom] = []
    for index, atom in enumerate(atoms):
        radius = vdw_radius_angstrom(atom.element)
        if radius is None:
            atom_id = atom.source_atom_id or f"{scope}:{atom.name}:{index}"
            diagnostics.append(
                _diagnostic(
                    "pocket.volume.unknown_radius",
                    f"No validated pocket radius is available for element {atom.element!r}; volume was not measured.",
                    severity=DiagnosticSeverity.ERROR,
                    source_id=scope,
                    atom_id=atom_id,
                    remediation="Correct the element annotation or provide a validated radius before measuring volume.",
                )
            )
            continue
        coordinate = (float(atom.coordinate[0]), float(atom.coordinate[1]), float(atom.coordinate[2]))
        known.append(_ExclusionAtom(coordinate, radius))
    return tuple(known)


def _active_exclusion_atoms(
    atoms: Sequence[_ExclusionAtom],
    lower: np.ndarray,
    upper: np.ndarray,
) -> tuple[_ExclusionAtom, ...]:
    """Keep atom spheres whose bounds can intersect the pocket envelope."""

    return tuple(
        atom
        for atom in atoms
        if all(
            coordinate + atom.radius_angstrom >= float(low)
            and coordinate - atom.radius_angstrom <= float(high)
            for coordinate, low, high in zip(atom.coordinate, lower, upper, strict=True)
        )
    )


def _build_exclusion_index(
    atoms: Sequence[_ExclusionAtom],
) -> tuple[cKDTree | None, tuple[_ExclusionAtom, ...]]:
    ordered = tuple(atoms)
    if not ordered:
        return None, ()
    coordinates = np.asarray([atom.coordinate for atom in ordered], dtype=np.float64)
    return cKDTree(coordinates), ordered


def _bounds(
    spheres: Sequence[Any],
    margin: float,
) -> tuple[np.ndarray, np.ndarray]:
    centers = np.asarray([sphere.center_xyz for sphere in spheres], dtype=np.float64)
    radii = np.asarray([sphere.radius_angstrom for sphere in spheres], dtype=np.float64)
    lower = np.min(centers - radii[:, None], axis=0) - margin
    upper = np.max(centers + radii[:, None], axis=0) + margin
    return lower, upper


def _shape_for_bounds(lower: np.ndarray, upper: np.ndarray, spacing: float) -> tuple[int, int, int]:
    shape: list[int] = []
    for low, high in zip(lower, upper, strict=True):
        quotient = (float(high) - float(low)) / spacing
        # Coordinates produced by subtraction can sit a few ulps above an
        # integral boundary.  Treat that representation noise as the boundary
        # itself, while retaining genuinely larger extents.
        tolerance = 1.0e-12 * max(1.0, abs(quotient))
        shape.append(max(1, int(math.ceil(quotient - tolerance))))
    return (shape[0], shape[1], shape[2])


def _measure_grid(
    origin: tuple[float, float, float],
    shape: tuple[int, int, int],
    spacing: float,
    spheres: Sequence[Any],
    tree: cKDTree | None,
    exclusion_atoms: Sequence[_ExclusionAtom],
    chunk_size: int,
    max_exclusion_neighbor_checks: int,
) -> _GridMeasurement:
    total = math.prod(shape)
    count = 0
    nx, ny, _ = shape
    plane = nx * ny
    sphere_data = tuple((np.asarray(sphere.center_xyz, dtype=np.float64), sphere.radius_angstrom**2) for sphere in spheres)
    exclusion_neighbor_checks = 0
    for start in range(0, total, chunk_size):
        stop = min(total, start + chunk_size)
        linear = np.arange(start, stop, dtype=np.int64)
        ix = linear % nx
        iy = (linear // nx) % ny
        iz = linear // plane
        points = np.column_stack(
            (
                origin[0] + (ix + 0.5) * spacing,
                origin[1] + (iy + 0.5) * spacing,
                origin[2] + (iz + 0.5) * spacing,
            )
        )
        inside = np.zeros(len(points), dtype=bool)
        for center, radius_squared in sphere_data:
            delta = points - center
            inside |= np.einsum("ij,ij->i", delta, delta) <= radius_squared
        if tree is not None and np.any(inside):
            inside_indices = np.flatnonzero(inside)
            query_radius = max(atom.radius_angstrom for atom in exclusion_atoms) + 1.0e-12
            for batch_start in range(0, len(inside_indices), 256):
                batch_indices = inside_indices[batch_start : batch_start + 256]
                batch_points = points[batch_indices]
                neighbor_counts = tree.query_ball_point(
                    batch_points,
                    r=query_radius,
                    return_length=True,
                )
                batch_checks = int(np.sum(neighbor_counts))
                if exclusion_neighbor_checks + batch_checks > max_exclusion_neighbor_checks:
                    raise _VolumeResourceLimit
                exclusion_neighbor_checks += batch_checks
                nearby = tree.query_ball_point(batch_points, r=query_radius)
                for local_index, candidates in zip(batch_indices, nearby, strict=True):
                    point = points[int(local_index)]
                    for atom_index in candidates:
                        atom = exclusion_atoms[int(atom_index)]
                        delta = point - np.asarray(atom.coordinate, dtype=np.float64)
                        if float(np.dot(delta, delta)) <= atom.radius_angstrom**2:
                            inside[int(local_index)] = False
                            break
        count += int(np.count_nonzero(inside))
    envelope_volume = float(total) * spacing**3
    return _GridMeasurement(count, shape, envelope_volume, exclusion_neighbor_checks)


def _resource_diagnostic(message: str, remediation: str) -> Diagnostic:
    return _diagnostic(
        "pocket.volume.resource_limit",
        message,
        severity=DiagnosticSeverity.ERROR,
        remediation=remediation,
    )


def _resource_result_without_grid(
    candidate: PocketCandidate,
    settings: PocketVolumeSettings,
    diagnostics: Sequence[Diagnostic],
) -> PocketVolumeResult:
    return PocketVolumeResult(
        availability=Availability.NUMERICAL_FAILURE,
        units=_VOLUME_UNITS,
        provenance_parameters=_provenance_parameters(candidate, settings),
        diagnostics=tuple(diagnostics),
        candidate_id=candidate.candidate_id,
        settings=settings,
    )


def _provenance_parameters(
    candidate: PocketCandidate,
    settings: PocketVolumeSettings,
) -> dict[str, object]:
    return {
        "radii_version": settings.radii_version,
        "candidate_id": candidate.candidate_id,
        "sphere_ids": tuple(sphere.sphere_id for sphere in candidate.alpha_spheres),
        "sphere_count": len(candidate.alpha_spheres),
        "coarse_grid_spacing_angstrom": settings.coarse_grid_spacing_angstrom,
        "fine_grid_spacing_angstrom": settings.fine_grid_spacing_angstrom,
        "boundary_margin_angstrom": settings.boundary_margin_angstrom,
        "component_exclusion_policy": settings.component_exclusion_policy,
        "grid_phase": _GRID_PHASE,
        "grid_sampling": _GRID_SAMPLING,
    }


def _result(
    candidate: PocketCandidate,
    settings: PocketVolumeSettings,
    availability: Availability,
    origin: tuple[float, float, float],
    coarse_shape: tuple[int, int, int],
    fine_shape: tuple[int, int, int],
    diagnostics: Sequence[Diagnostic],
    coarse_count: int | None,
    fine_count: int | None,
    coarse_volume: float | None,
    fine_volume: float | None,
    sensitivity: PocketVolumeSensitivity | None,
    exclusion_atoms: Sequence[_ExclusionAtom],
) -> PocketVolumeResult:
    parameters = _provenance_parameters(candidate, settings)
    parameters["active_exclusion_atom_count"] = len(exclusion_atoms)
    parameters["rotation_error_bound"] = "two_times_half_voxel_diagonal_boundary_shells"
    signature = (
        settings.radii_version,
        candidate.candidate_id,
        tuple(sphere.sphere_id for sphere in candidate.alpha_spheres),
        settings.coarse_grid_spacing_angstrom,
        settings.fine_grid_spacing_angstrom,
        settings.boundary_margin_angstrom,
        settings.component_exclusion_policy,
        _GRID_PHASE,
        _GRID_SAMPLING,
    )
    half_diagonal = math.sqrt(3.0) * settings.fine_grid_spacing_angstrom / 2.0
    boundary_radii = tuple(sphere.radius_angstrom for sphere in candidate.alpha_spheres) + tuple(
        atom.radius_angstrom for atom in exclusion_atoms
    )
    error_bound = 2.0 * (4.0 / 3.0) * math.pi * sum(
        (radius + half_diagonal) ** 3 - max(0.0, radius - half_diagonal) ** 3
        for radius in boundary_radii
    )
    return PocketVolumeResult(
        availability=availability,
        coarse_voxel_count=coarse_count,
        fine_voxel_count=fine_count,
        coarse_volume_angstrom3=coarse_volume,
        fine_volume_angstrom3=fine_volume,
        grid_shape=coarse_shape,
        grid_origin_xyz=origin,
        grid_phase=_GRID_PHASE,
        coarse_grid_shape=coarse_shape,
        fine_grid_shape=fine_shape,
        rotation_error_bound_angstrom3=error_bound if availability is Availability.AVAILABLE else None,
        sensitivity=sensitivity,
        units=_VOLUME_UNITS,
        provenance_parameters=parameters,
        diagnostics=tuple(diagnostics),
        candidate_id=candidate.candidate_id,
        compatibility_signature=signature,
        settings=settings,
    )


def compare_pocket_volumes(
    reference: PocketVolumeResult,
    target: PocketVolumeResult,
) -> PocketVolumeComparison:
    """Compare compatible measurements as target minus reference."""

    if not isinstance(reference, PocketVolumeResult) or not isinstance(target, PocketVolumeResult):
        raise TypeError("reference and target must be PocketVolumeResult values")
    if reference.compatibility_signature != target.compatibility_signature:
        diagnostic = _diagnostic(
            "pocket.volume.incompatible_settings",
            "Pocket volume deltas are disabled because method-defining inputs differ.",
            remediation="Use the same radii table, sphere selection, grid settings, and component policy.",
        )
        return PocketVolumeComparison(
            availability=Availability.NOT_APPLICABLE,
            diagnostics=(diagnostic,),
            units=_VOLUME_UNITS,
        )
    if reference.availability is not Availability.AVAILABLE or target.availability is not Availability.AVAILABLE:
        availability = (
            Availability.NUMERICAL_FAILURE
            if Availability.NUMERICAL_FAILURE in {reference.availability, target.availability}
            else Availability.NOT_APPLICABLE
        )
        return PocketVolumeComparison(
            availability=availability,
            reference_volume_angstrom3=reference.fine_volume_angstrom3,
            target_volume_angstrom3=target.fine_volume_angstrom3,
            diagnostics=reference.diagnostics + target.diagnostics,
            units=_VOLUME_UNITS,
            compatibility_signature=reference.compatibility_signature,
        )
    assert reference.fine_volume_angstrom3 is not None
    assert target.fine_volume_angstrom3 is not None
    delta = target.fine_volume_angstrom3 - reference.fine_volume_angstrom3
    relative = delta / reference.fine_volume_angstrom3 if reference.fine_volume_angstrom3 > 0.0 else None
    return PocketVolumeComparison(
        availability=Availability.AVAILABLE,
        delta_angstrom3=delta,
        relative_delta_fraction=relative,
        reference_volume_angstrom3=reference.fine_volume_angstrom3,
        target_volume_angstrom3=target.fine_volume_angstrom3,
        units=_VOLUME_UNITS,
        compatibility_signature=reference.compatibility_signature,
    )


__all__ = [
    "PocketVolumeComparison",
    "PocketVolumeResult",
    "PocketVolumeSensitivity",
    "PocketVolumeSensitivityResult",
    "PocketVolumeSettings",
    "compare_pocket_volumes",
    "measure_pocket_volume",
    "sensitivity_from_volumes",
]
