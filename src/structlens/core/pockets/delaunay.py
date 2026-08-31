"""Deterministic Delaunay alpha-sphere generation.

This is a small clean-room implementation of the geometric alpha-sphere idea.
It only creates evidence from finite, known-radius atoms and refuses to run an
unbounded tessellation.  It does not assign a ligandability probability.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass

import numpy as np
from scipy.spatial import Delaunay, QhullError, cKDTree  # type: ignore[import-untyped]

from structlens.core.evidence import Diagnostic, DiagnosticSeverity
from structlens.core.models import ResidueId

from .geometry import tetrahedron_circumsphere
from .models import AlphaSphere, AlphaSphereDetectionResult, PocketDetectionSettings, PocketGeometrySettings
from .radii import vdw_radius_angstrom

_GEOMETRY_TOLERANCE = 1.0e-9


def _diagnostic(code: str, message: str, severity: DiagnosticSeverity) -> Diagnostic:
    return Diagnostic(
        code=code,
        severity=severity,
        message=message,
        remediation="Review the selected polymer atoms and detector settings before retrying.",
    )


@dataclass(frozen=True, slots=True)
class PocketAtom:
    """Immutable atom input for pocket geometry.

    An element absent from the selected radius table is retained as evidence
    but excluded from tessellation; no arbitrary radius is substituted.
    """

    atom_id: str
    residue_id: ResidueId
    coordinate: tuple[float, float, float]
    radius_angstrom: float
    element: str | None = None

    def __post_init__(self) -> None:
        atom_id = str(self.atom_id).strip()
        if not atom_id:
            raise ValueError("atom_id must not be empty")
        if not isinstance(self.residue_id, ResidueId):
            raise TypeError("residue_id must be a ResidueId")
        try:
            coordinate = tuple(float(value) for value in self.coordinate)
        except (TypeError, ValueError) as exc:
            raise ValueError("atom coordinate must contain numeric values") from exc
        if len(coordinate) != 3 or any(not math.isfinite(value) for value in coordinate):
            raise ValueError("atom coordinate must contain three finite values")
        radius = float(self.radius_angstrom)
        if not math.isfinite(radius) or radius <= 0.0:
            raise ValueError("radius_angstrom must be finite and positive")
        element = None if self.element is None else str(self.element).strip().upper() or None
        object.__setattr__(self, "atom_id", atom_id)
        object.__setattr__(self, "coordinate", coordinate)
        object.__setattr__(self, "radius_angstrom", radius)
        object.__setattr__(self, "element", element)

    @property
    def has_known_radius(self) -> bool:
        return self.element is None or vdw_radius_angstrom(self.element) is not None

    def to_json(self) -> dict[str, object]:
        return {
            "atom_id": self.atom_id,
            "residue_id": {
                "structure_id": self.residue_id.structure_id,
                "model_id": self.residue_id.model_id,
                "chain_id": self.residue_id.chain_id,
                "auth_seq_id": self.residue_id.auth_seq_id,
                "insertion_code": self.residue_id.insertion_code,
                "residue_name": self.residue_id.residue_name,
            },
            "coordinate": list(self.coordinate),
            "radius_angstrom": self.radius_angstrom,
            "element": self.element,
        }


def _normalize_atoms(atoms: Iterable[PocketAtom]) -> tuple[PocketAtom, ...]:
    values = tuple(atoms)
    if any(not isinstance(item, PocketAtom) for item in values):
        raise TypeError("atoms must contain PocketAtom values")
    return tuple(sorted(values, key=lambda item: item.atom_id))


@dataclass(frozen=True, slots=True)
class _AtomSpatialIndex:
    atoms: tuple[PocketAtom, ...]
    tree: cKDTree
    maximum_radius_angstrom: float

    @classmethod
    def build(cls, atoms: Sequence[PocketAtom]) -> _AtomSpatialIndex:
        values = tuple(atoms)
        if not values:
            raise ValueError("all_atoms must contain the simplex atoms")
        coordinates = np.asarray([atom.coordinate for atom in values], dtype=np.float64)
        return cls(values, cKDTree(coordinates), max(atom.radius_angstrom for atom in values))

    def surface_neighbors(
        self,
        center: tuple[float, float, float],
        maximum_surface_clearance: float,
    ) -> tuple[PocketAtom, ...]:
        center_array = np.asarray(center, dtype=np.float64)
        center_cutoff = max(
            0.0,
            maximum_surface_clearance + self.maximum_radius_angstrom + _GEOMETRY_TOLERANCE,
        )
        indices = self.tree.query_ball_point(center_array, center_cutoff)
        return tuple(self.atoms[int(index)] for index in sorted(int(item) for item in indices))


def _residue_sort_key(item: ResidueId) -> tuple[str, str, str, str]:
    return (item.chain_id, item.auth_seq_id, item.insertion_code or "", item.residue_name)


def _alpha_sphere_from_index(
    simplex_values: tuple[PocketAtom, ...],
    spatial_index: _AtomSpatialIndex,
    settings: PocketGeometrySettings,
) -> tuple[AlphaSphere | None, int]:
    simplex_ids = tuple(item.atom_id for item in simplex_values)
    circumsphere = tetrahedron_circumsphere(tuple(item.coordinate for item in simplex_values))
    if circumsphere.sphere is None:
        return None, 0
    center = circumsphere.sphere.center_xyz
    vertex_clearance = min(
        _distance_from_point(center, atom.coordinate) - atom.radius_angstrom
        for atom in simplex_values
    )
    if not (
        settings.minimum_alpha_sphere_radius_angstrom
        <= vertex_clearance
        <= settings.maximum_alpha_sphere_radius_angstrom
    ):
        return None, 0

    relevant = spatial_index.surface_neighbors(
        center,
        vertex_clearance + settings.lining_contact_slack_angstrom,
    )
    simplex_id_set = set(simplex_ids)
    surface_clearances = tuple(
        (
            atom,
            _distance_from_point(center, atom.coordinate) - atom.radius_angstrom,
        )
        for atom in relevant
    )
    if any(
        clearance < vertex_clearance - _GEOMETRY_TOLERANCE
        for atom, clearance in surface_clearances
        if atom.atom_id not in simplex_id_set
    ):
        return None, len(relevant)
    lining_residues = tuple(
        sorted(
            {
                atom.residue_id
                for atom, clearance in surface_clearances
                if clearance
                <= vertex_clearance + settings.lining_contact_slack_angstrom + _GEOMETRY_TOLERANCE
            },
            key=_residue_sort_key,
        )
    )
    if not lining_residues:
        lining_residues = tuple(sorted({item.residue_id for item in simplex_values}, key=_residue_sort_key))
    return (
        AlphaSphere(
            center_xyz=center,
            radius_angstrom=vertex_clearance,
            touching_atom_ids=simplex_ids,
            lining_residues=lining_residues,
            source_simplex_atom_ids=simplex_ids,
        ),
        len(relevant),
    )


def alpha_sphere_from_simplex(
    simplex: Sequence[PocketAtom],
    all_atoms: Sequence[PocketAtom],
    settings: PocketGeometrySettings,
) -> AlphaSphere | None:
    """Build one alpha sphere from four vertices, or return ``None``.

    The Delaunay circumsphere supplies the center, while the stored radius is
    the minimum atom-surface clearance ``distance - vdw_radius``.  The declared
    alpha range is applied to that physical clearance.  An indexed query also
    rejects a non-vertex surface that is closer than the defining vertices.
    """

    if not isinstance(settings, PocketGeometrySettings):
        raise TypeError("settings must be PocketGeometrySettings")
    simplex_values = tuple(simplex)
    all_values = tuple(all_atoms)
    if len(simplex_values) != 4:
        raise ValueError("simplex must contain exactly four PocketAtom values")
    if any(not isinstance(item, PocketAtom) for item in simplex_values + all_values):
        raise TypeError("simplex and all_atoms must contain PocketAtom values")
    simplex_ids = tuple(item.atom_id for item in simplex_values)
    if len(set(simplex_ids)) != 4:
        raise ValueError("simplex must contain four distinct atom IDs")
    if any(not item.has_known_radius for item in simplex_values):
        return None

    all_ids = {item.atom_id for item in all_values}
    if not set(simplex_ids).issubset(all_ids):
        raise ValueError("all_atoms must contain every simplex atom")
    known_atoms = tuple(item for item in all_values if item.has_known_radius)
    if not known_atoms:
        return None
    candidate, _ = _alpha_sphere_from_index(
        simplex_values,
        _AtomSpatialIndex.build(known_atoms),
        settings,
    )
    return candidate


def _distance_from_point(first: Sequence[float], second: Sequence[float]) -> float:
    return float(np.linalg.norm(np.asarray(first, dtype=np.float64) - np.asarray(second, dtype=np.float64)))


def _estimated_simplex_count(atom_count: int) -> int:
    if atom_count < 4:
        return 0
    # A 3-D Delaunay tessellation is the lower hull of points lifted into 4-D.
    # The upper-bound theorem therefore gives at most n(n-3)/2 tetrahedral
    # facets.  Using n-choose-4 would reject realistic proteins even though
    # those combinations are never materialized by Qhull.
    return atom_count * (atom_count - 3) // 2


def detect_alpha_spheres(
    atoms: Iterable[PocketAtom],
    settings: PocketDetectionSettings | None = None,
    *,
    cancel_check: Callable[[], bool] | None = None,
) -> AlphaSphereDetectionResult:
    """Generate deterministic alpha spheres from a bounded atom selection."""

    detector_settings = settings or PocketDetectionSettings()
    if not isinstance(detector_settings, PocketDetectionSettings):
        raise TypeError("settings must be PocketDetectionSettings")
    normalized = _normalize_atoms(atoms)
    if cancel_check is not None and cancel_check():
        return AlphaSphereDetectionResult(
            diagnostics=(
                _diagnostic(
                    "pocket.detect.cancelled",
                    "Pocket detection was cancelled before tessellation.",
                    DiagnosticSeverity.INFO,
                ),
            )
        )
    diagnostics: list[Diagnostic] = []
    seen: set[str] = set()
    unique: list[PocketAtom] = []
    for atom in normalized:
        if atom.atom_id in seen:
            diagnostics.append(
                _diagnostic(
                    "pocket.detect.duplicate_atom_id",
                    "Duplicate atom identities cannot be used for deterministic tessellation.",
                    DiagnosticSeverity.ERROR,
                )
            )
            continue
        seen.add(atom.atom_id)
        unique.append(atom)
        if not atom.has_known_radius:
            diagnostics.append(
                _diagnostic(
                    "pocket.detect.unknown_radius",
                    "An atom without a validated van der Waals radius was excluded from tessellation.",
                    DiagnosticSeverity.WARNING,
                )
            )

    if any(item.code == "pocket.detect.duplicate_atom_id" for item in diagnostics):
        return AlphaSphereDetectionResult(diagnostics=tuple(_sort_diagnostics(diagnostics)))

    if len(unique) > detector_settings.max_atom_count:
        diagnostics.append(
            _diagnostic(
                "pocket.detect.resource_limit",
                "The selected atom count exceeds the detector resource limit.",
                DiagnosticSeverity.ERROR,
            )
        )
        return AlphaSphereDetectionResult(diagnostics=tuple(_sort_diagnostics(diagnostics)))

    known = tuple(item for item in unique if item.has_known_radius)
    if len(known) < 4:
        diagnostics.append(
            _diagnostic(
                "pocket.detect.insufficient_atoms",
                "At least four atoms with validated radii are required for a 3D tessellation.",
                DiagnosticSeverity.WARNING,
            )
        )
        return AlphaSphereDetectionResult(diagnostics=tuple(_sort_diagnostics(diagnostics)))

    estimated = _estimated_simplex_count(len(known))
    if estimated > detector_settings.max_estimated_simplices:
        diagnostics.append(
            _diagnostic(
                "pocket.detect.resource_limit",
                "The estimated Delaunay simplex count exceeds the detector resource limit.",
                DiagnosticSeverity.ERROR,
            )
        )
        return AlphaSphereDetectionResult(diagnostics=tuple(_sort_diagnostics(diagnostics)))

    coordinates = np.asarray([item.coordinate for item in known], dtype=np.float64)
    try:
        tessellation = Delaunay(coordinates)
    except QhullError:
        diagnostics.append(
            _diagnostic(
                "pocket.detect.qhull_failure",
                "The selected atom coordinates could not be tessellated as a stable 3D hull.",
                DiagnosticSeverity.ERROR,
            )
        )
        return AlphaSphereDetectionResult(diagnostics=tuple(_sort_diagnostics(diagnostics)))
    except (MemoryError, ValueError):
        diagnostics.append(
            _diagnostic(
                "pocket.detect.qhull_failure",
                "The selected atom coordinates could not be tessellated within numerical limits.",
                DiagnosticSeverity.ERROR,
            )
        )
        return AlphaSphereDetectionResult(diagnostics=tuple(_sort_diagnostics(diagnostics)))

    spheres: dict[str, AlphaSphere] = {}
    spatial_index = _AtomSpatialIndex.build(known)
    clearance_atom_checks = 0
    simplex_indices = sorted(tuple(sorted(int(index) for index in simplex)) for simplex in tessellation.simplices)
    for simplex_index in simplex_indices:
        if cancel_check is not None and cancel_check():
            diagnostics.append(
                _diagnostic(
                    "pocket.detect.cancelled",
                    "Pocket detection was cancelled during alpha-sphere generation.",
                    DiagnosticSeverity.INFO,
                )
            )
            return AlphaSphereDetectionResult(diagnostics=tuple(_sort_diagnostics(diagnostics)))
        simplex = tuple(known[index] for index in simplex_index)
        candidate, atom_checks = _alpha_sphere_from_index(
            simplex,
            spatial_index,
            detector_settings.geometry,
        )
        clearance_atom_checks += atom_checks
        if clearance_atom_checks > detector_settings.max_clearance_atom_checks:
            diagnostics.append(
                _diagnostic(
                    "pocket.detect.resource_limit",
                    "The indexed atom-clearance work exceeds the detector resource limit.",
                    DiagnosticSeverity.ERROR,
                )
            )
            return AlphaSphereDetectionResult(diagnostics=tuple(_sort_diagnostics(diagnostics)))
        if candidate is not None:
            spheres[candidate.sphere_id] = candidate

    return AlphaSphereDetectionResult(
        spheres=tuple(spheres.values()),
        diagnostics=tuple(_sort_diagnostics(diagnostics)),
    )


def _sort_diagnostics(values: Iterable[Diagnostic]) -> list[Diagnostic]:
    return sorted(values, key=lambda item: (item.code, item.source_id or "", item.atom_id or "", item.residue_id or ""))


# Descriptive aliases keep the public vocabulary discoverable without creating
# separate result contracts.
generate_alpha_spheres = detect_alpha_spheres
detect_pocket_spheres = detect_alpha_spheres


__all__ = [
    "PocketAtom",
    "alpha_sphere_from_simplex",
    "detect_alpha_spheres",
    "detect_pocket_spheres",
    "generate_alpha_spheres",
]
