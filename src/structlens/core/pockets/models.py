"""Immutable pocket geometry primitives.

Task 6 only establishes validated, serialization-ready geometry contracts.
Detection, measurement, and comparison logic arrive in later tasks.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from dataclasses import dataclass, field

from structlens.core.evidence import Availability, Diagnostic
from structlens.core.models import ResidueId

JSONScalar = str | int | float | bool | None
JSONValue = JSONScalar | list["JSONValue"] | dict[str, "JSONValue"]
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def _canonical_bytes(payload: dict[str, JSONValue]) -> bytes:
    return json.dumps(
        payload,
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")


def _finite(value: float, name: str, *, positive: bool = False, non_negative: bool = False) -> float:
    numeric = float(value)
    if not math.isfinite(numeric):
        raise ValueError(f"{name} must be finite")
    if positive and numeric <= 0.0:
        raise ValueError(f"{name} must be positive")
    if non_negative and numeric < 0.0:
        raise ValueError(f"{name} must be non-negative")
    return numeric


def _positive_int(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _atom_ids(values: tuple[str, ...], name: str, *, expected_count: int | None = None) -> tuple[str, ...]:
    normalized = tuple(sorted(str(item).strip() for item in values))
    if not normalized or any(not item for item in normalized):
        raise ValueError(f"{name} must contain non-empty atom IDs")
    if len(set(normalized)) != len(normalized):
        raise ValueError(f"{name} must not contain duplicates")
    if expected_count is not None and len(normalized) != expected_count:
        raise ValueError(f"{name} must contain exactly {expected_count} atom IDs")
    return normalized


def _residue_sort_key(residue: ResidueId) -> tuple[str, int, str, str, str, str, str]:
    auth_seq_id = residue.auth_seq_id.strip()
    try:
        return (
            residue.chain_id,
            0,
            f"{int(auth_seq_id):012d}",
            auth_seq_id,
            residue.insertion_code or "",
            residue.residue_name,
            residue.model_id,
        )
    except ValueError:
        return (
            residue.chain_id,
            1,
            auth_seq_id,
            auth_seq_id,
            residue.insertion_code or "",
            residue.residue_name,
            residue.model_id,
        )


def _residue_json(residue: ResidueId) -> dict[str, JSONValue]:
    return {
        "structure_id": residue.structure_id,
        "model_id": residue.model_id,
        "chain_id": residue.chain_id,
        "auth_seq_id": residue.auth_seq_id,
        "insertion_code": residue.insertion_code,
        "residue_name": residue.residue_name,
    }


@dataclass(frozen=True, slots=True)
class PocketGeometrySettings:
    """Physical pocket settings; alpha radii are free VDW clearances."""

    minimum_alpha_sphere_radius_angstrom: float = 2.8
    maximum_alpha_sphere_radius_angstrom: float = 6.2
    probe_radius_angstrom: float = 1.4
    lining_contact_slack_angstrom: float = 0.5

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "minimum_alpha_sphere_radius_angstrom",
            _finite(
                self.minimum_alpha_sphere_radius_angstrom,
                "minimum_alpha_sphere_radius_angstrom",
                positive=True,
            ),
        )
        object.__setattr__(
            self,
            "maximum_alpha_sphere_radius_angstrom",
            _finite(
                self.maximum_alpha_sphere_radius_angstrom,
                "maximum_alpha_sphere_radius_angstrom",
                positive=True,
            ),
        )
        object.__setattr__(
            self,
            "probe_radius_angstrom",
            _finite(self.probe_radius_angstrom, "probe_radius_angstrom", positive=True),
        )
        object.__setattr__(
            self,
            "lining_contact_slack_angstrom",
            _finite(
                self.lining_contact_slack_angstrom,
                "lining_contact_slack_angstrom",
                non_negative=True,
            ),
        )
        if self.minimum_alpha_sphere_radius_angstrom > self.maximum_alpha_sphere_radius_angstrom:
            raise ValueError(
                "minimum_alpha_sphere_radius_angstrom must not exceed maximum_alpha_sphere_radius_angstrom"
            )

    def to_json(self) -> dict[str, JSONValue]:
        return {
            "minimum_alpha_sphere_radius_angstrom": self.minimum_alpha_sphere_radius_angstrom,
            "maximum_alpha_sphere_radius_angstrom": self.maximum_alpha_sphere_radius_angstrom,
            "probe_radius_angstrom": self.probe_radius_angstrom,
            "lining_contact_slack_angstrom": self.lining_contact_slack_angstrom,
        }


@dataclass(frozen=True, slots=True)
class PocketDetectionSettings:
    """Bounded, reproducible settings for blind pocket detection.

    The detector is deliberately parameterized at the domain boundary.  In
    particular, Delaunay tessellation is never allowed to allocate based only
    on untrusted input size; both atom and estimated-simplex caps are applied
    before SciPy is called.
    """

    geometry: PocketGeometrySettings = field(default_factory=PocketGeometrySettings)
    cluster_distance_padding_angstrom: float = 1.0
    solvent_grid_spacing_angstrom: float = 0.5
    solvent_boundary_margin_angstrom: float = 4.0
    minimum_cluster_size: int = 2
    maximum_candidates: int = 20
    max_atom_count: int = 50_000
    max_estimated_simplices: int = 2_000_000
    max_solvent_grid_cells: int = 2_000_000
    max_clearance_atom_checks: int = 20_000_000
    max_solvent_raster_cells: int = 20_000_000

    def __post_init__(self) -> None:
        if not isinstance(self.geometry, PocketGeometrySettings):
            raise TypeError("geometry must be PocketGeometrySettings")
        for name in (
            "cluster_distance_padding_angstrom",
            "solvent_boundary_margin_angstrom",
        ):
            value = _finite(getattr(self, name), name, non_negative=True)
            object.__setattr__(self, name, value)
        spacing = _finite(self.solvent_grid_spacing_angstrom, "solvent_grid_spacing_angstrom", positive=True)
        object.__setattr__(self, "solvent_grid_spacing_angstrom", spacing)
        for name in (
            "minimum_cluster_size",
            "maximum_candidates",
            "max_atom_count",
            "max_estimated_simplices",
            "max_solvent_grid_cells",
            "max_clearance_atom_checks",
            "max_solvent_raster_cells",
        ):
            object.__setattr__(self, name, _positive_int(getattr(self, name), name))

    def to_json(self) -> dict[str, JSONValue]:
        """Return a fresh JSON-compatible copy of the detector policy."""

        return {
            "geometry": self.geometry.to_json(),
            "cluster_distance_padding_angstrom": self.cluster_distance_padding_angstrom,
            "solvent_grid_spacing_angstrom": self.solvent_grid_spacing_angstrom,
            "solvent_boundary_margin_angstrom": self.solvent_boundary_margin_angstrom,
            "minimum_cluster_size": self.minimum_cluster_size,
            "maximum_candidates": self.maximum_candidates,
            "max_atom_count": self.max_atom_count,
            "max_estimated_simplices": self.max_estimated_simplices,
            "max_solvent_grid_cells": self.max_solvent_grid_cells,
            "max_clearance_atom_checks": self.max_clearance_atom_checks,
            "max_solvent_raster_cells": self.max_solvent_raster_cells,
        }


@dataclass(frozen=True, slots=True)
class AlphaSphere:
    """A Delaunay center with ``radius_angstrom`` equal to free VDW clearance."""

    center_xyz: tuple[float, float, float]
    radius_angstrom: float
    touching_atom_ids: tuple[str, ...]
    lining_residues: tuple[ResidueId, ...]
    source_simplex_atom_ids: tuple[str, ...]
    sphere_id: str = field(init=False)

    def __post_init__(self) -> None:
        try:
            center = tuple(float(value) for value in self.center_xyz)
        except (TypeError, ValueError) as exc:
            raise ValueError("center_xyz must contain numeric values") from exc
        if len(center) != 3 or any(not math.isfinite(value) for value in center):
            raise ValueError("center_xyz must contain exactly three finite values")
        object.__setattr__(self, "center_xyz", center)
        object.__setattr__(self, "radius_angstrom", _finite(self.radius_angstrom, "radius_angstrom", positive=True))
        object.__setattr__(
            self,
            "touching_atom_ids",
            _atom_ids(tuple(self.touching_atom_ids), "touching_atom_ids", expected_count=4),
        )
        simplex_ids = _atom_ids(tuple(self.source_simplex_atom_ids), "source_simplex_atom_ids", expected_count=4)
        object.__setattr__(self, "source_simplex_atom_ids", simplex_ids)
        residues = tuple(self.lining_residues)
        if not residues or any(not isinstance(item, ResidueId) for item in residues):
            raise ValueError("lining_residues must contain ResidueId values")
        residues = tuple(sorted(set(residues), key=_residue_sort_key))
        object.__setattr__(self, "lining_residues", residues)
        payload = self._payload(include_id=False)
        object.__setattr__(self, "sphere_id", hashlib.sha256(_canonical_bytes(payload)).hexdigest())

    def _payload(self, *, include_id: bool) -> dict[str, JSONValue]:
        payload: dict[str, JSONValue] = {
            "center_xyz": list(self.center_xyz),
            "radius_angstrom": self.radius_angstrom,
            "touching_atom_ids": list(self.touching_atom_ids),
            "lining_residues": [_residue_json(item) for item in self.lining_residues],
            "source_simplex_atom_ids": list(self.source_simplex_atom_ids),
        }
        if include_id:
            payload["sphere_id"] = self.sphere_id
        return payload

    def to_json(self) -> dict[str, JSONValue]:
        return self._payload(include_id=True)

    def canonical_json_bytes(self) -> bytes:
        return _canonical_bytes(self.to_json())


@dataclass(frozen=True, slots=True)
class PocketCandidate:
    alpha_spheres: tuple[AlphaSphere, ...]
    source_content_id: str | None = None
    selection_id: str | None = None
    centroid_xyz: tuple[float, float, float] = field(init=False)
    lining_residues: tuple[ResidueId, ...] = field(init=False)
    touching_atom_ids: tuple[str, ...] = field(init=False)
    candidate_id: str = field(init=False)

    def __post_init__(self) -> None:
        spheres = tuple(self.alpha_spheres)
        if not spheres or any(not isinstance(item, AlphaSphere) for item in spheres):
            raise ValueError("alpha_spheres must contain AlphaSphere values")
        spheres = tuple(sorted(spheres, key=lambda item: item.sphere_id))
        object.__setattr__(self, "alpha_spheres", spheres)
        centroid = tuple(
            sum(sphere.center_xyz[index] for sphere in spheres) / len(spheres) for index in range(3)
        )
        object.__setattr__(self, "centroid_xyz", centroid)
        residues = tuple(
            sorted(
                {residue for sphere in spheres for residue in sphere.lining_residues},
                key=_residue_sort_key,
            )
        )
        object.__setattr__(self, "lining_residues", residues)
        atom_ids = tuple(sorted({atom_id for sphere in spheres for atom_id in sphere.touching_atom_ids}))
        object.__setattr__(self, "touching_atom_ids", atom_ids)
        if (self.source_content_id is None) != (self.selection_id is None):
            raise ValueError("source_content_id and selection_id must be provided together")
        for name in ("source_content_id", "selection_id"):
            value = getattr(self, name)
            if value is not None:
                token = str(value).strip().lower()
                if _SHA256_RE.fullmatch(token) is None:
                    raise ValueError(f"{name} must be a lowercase SHA-256 digest")
                object.__setattr__(self, name, token)
        payload = self._payload(include_id=False)
        object.__setattr__(self, "candidate_id", hashlib.sha256(_canonical_bytes(payload)).hexdigest())

    def _payload(self, *, include_id: bool) -> dict[str, JSONValue]:
        payload: dict[str, JSONValue] = {
            "alpha_spheres": [sphere.to_json() for sphere in self.alpha_spheres],
            "centroid_xyz": list(self.centroid_xyz),
            "lining_residues": [_residue_json(item) for item in self.lining_residues],
            "touching_atom_ids": list(self.touching_atom_ids),
        }
        if include_id:
            payload["candidate_id"] = self.candidate_id
            payload["lineage"] = {
                "source_content_id": self.source_content_id,
                "selection_id": self.selection_id,
            }
        return payload

    def to_json(self) -> dict[str, JSONValue]:
        return self._payload(include_id=True)

    def canonical_json_bytes(self) -> bytes:
        return _canonical_bytes(self.to_json())


@dataclass(frozen=True, slots=True)
class AlphaSphereDetectionResult:
    """Immutable alpha-sphere output with explicit no-pocket semantics."""

    spheres: tuple[AlphaSphere, ...] = field(default_factory=tuple)
    diagnostics: tuple[Diagnostic, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        spheres = tuple(self.spheres)
        if any(not isinstance(item, AlphaSphere) for item in spheres):
            raise TypeError("spheres must contain AlphaSphere values")
        diagnostics = tuple(self.diagnostics)
        if any(not isinstance(item, Diagnostic) for item in diagnostics):
            raise TypeError("diagnostics must contain Diagnostic values")
        object.__setattr__(self, "spheres", tuple(sorted(spheres, key=lambda item: item.sphere_id)))
        object.__setattr__(self, "diagnostics", diagnostics)

    @property
    def availability(self) -> Availability:
        codes = {item.code for item in self.diagnostics}
        if "pocket.detect.cancelled" in codes:
            return Availability.NOT_APPLICABLE
        if codes & {"pocket.detect.resource_limit", "pocket.detect.qhull_failure"}:
            return Availability.NUMERICAL_FAILURE
        if codes & {"pocket.detect.insufficient_atoms", "pocket.detect.duplicate_atom_id"}:
            return Availability.INVALID_INPUT
        if self.spheres:
            return Availability.AVAILABLE
        return Availability.NOT_DETECTED

    @property
    def status(self) -> Availability:
        """Compatibility alias for consumers that call state ``status``."""

        return self.availability

    @property
    def has_spheres(self) -> bool:
        return bool(self.spheres)

    def to_json(self) -> dict[str, object]:
        return {
            "availability": self.availability.value,
            "spheres": [item.to_json() for item in self.spheres],
            "diagnostics": [item.to_json() for item in self.diagnostics],
        }


__all__ = [
    "AlphaSphere",
    "AlphaSphereDetectionResult",
    "PocketCandidate",
    "PocketDetectionSettings",
    "PocketGeometrySettings",
]
