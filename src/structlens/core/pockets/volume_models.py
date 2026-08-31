"""Immutable contracts and validation helpers for pocket free-volume runs.

The measurement implementation lives in :mod:`volume`; this module keeps the
public data contracts focused and dependency-light so callers can construct,
serialize, and compare results without importing numerical machinery.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any, cast

from structlens.core.evidence import Availability, Diagnostic, DiagnosticSeverity
from structlens.core.provenance import MethodProvenance

from .radii import POCKET_RADII_VERSION

_GRID_PHASE = "cell_center"
_GRID_SAMPLING = "origin + (index + 0.5) * spacing"
_VOLUME_UNITS = MappingProxyType({"volume": "angstrom^3", "length": "angstrom"})
_COMPONENT_POLICIES = frozenset({"protein_only", "unoccupied"})


def _finite(value: object, name: str, *, positive: bool = False, non_negative: bool = False) -> float:
    try:
        converted = float(cast(Any, value))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(converted):
        raise ValueError(f"{name} must be finite")
    if positive and converted <= 0.0:
        raise ValueError(f"{name} must be positive")
    if non_negative and converted < 0.0:
        raise ValueError(f"{name} must be non-negative")
    return converted


def _positive_int(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _freeze_json(value: Any) -> Any:
    """Recursively freeze a JSON-like value for immutable result contracts."""

    if isinstance(value, Mapping):
        frozen: dict[str, Any] = {}
        for key, item in value.items():
            if not isinstance(key, str) or not key.strip():
                raise ValueError("mapping keys must be non-empty strings")
            frozen[key] = _freeze_json(item)
        return MappingProxyType(frozen)
    if isinstance(value, (list, tuple)):
        return tuple(_freeze_json(item) for item in value)
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError("JSON values must be finite")
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    raise TypeError(f"unsupported JSON value type: {type(value).__name__}")


def _thaw_json(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _thaw_json(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_thaw_json(item) for item in value]
    return value


def _diagnostic(code: str, message: str, *, remediation: str | None = None) -> Diagnostic:
    return Diagnostic(
        code=code,
        severity=DiagnosticSeverity.ERROR if code.endswith("resource_limit") else DiagnosticSeverity.WARNING,
        message=message,
        remediation=remediation,
    )


@dataclass(frozen=True, slots=True)
class PocketVolumeSettings:
    """Validated settings for a dual-resolution pocket free-volume run."""

    coarse_grid_spacing_angstrom: float = 1.0
    fine_grid_spacing_angstrom: float = 0.5
    boundary_margin_angstrom: float = 0.0
    component_exclusion_policy: str = "protein_only"
    max_voxel_count: int = 2_000_000
    voxel_chunk_size: int = 16_384
    radii_version: str = POCKET_RADII_VERSION

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "coarse_grid_spacing_angstrom",
            _finite(self.coarse_grid_spacing_angstrom, "coarse_grid_spacing_angstrom", positive=True),
        )
        object.__setattr__(
            self,
            "fine_grid_spacing_angstrom",
            _finite(self.fine_grid_spacing_angstrom, "fine_grid_spacing_angstrom", positive=True),
        )
        object.__setattr__(
            self,
            "boundary_margin_angstrom",
            _finite(self.boundary_margin_angstrom, "boundary_margin_angstrom", non_negative=True),
        )
        policy = str(self.component_exclusion_policy).strip()
        if policy not in _COMPONENT_POLICIES:
            raise ValueError(
                "component_exclusion_policy must be one of: protein_only, unoccupied"
            )
        object.__setattr__(self, "component_exclusion_policy", policy)
        object.__setattr__(self, "max_voxel_count", _positive_int(self.max_voxel_count, "max_voxel_count"))
        object.__setattr__(self, "voxel_chunk_size", _positive_int(self.voxel_chunk_size, "voxel_chunk_size"))
        radii_version = str(self.radii_version).strip()
        if radii_version != POCKET_RADII_VERSION:
            raise ValueError("radii_version must identify the active pocket radii table")
        object.__setattr__(self, "radii_version", radii_version)

    @property
    def compatibility_parameters(self) -> Mapping[str, object]:
        """Method-defining settings, excluding resource/execution controls."""

        return MappingProxyType(
            {
                "coarse_grid_spacing_angstrom": self.coarse_grid_spacing_angstrom,
                "fine_grid_spacing_angstrom": self.fine_grid_spacing_angstrom,
                "boundary_margin_angstrom": self.boundary_margin_angstrom,
                "component_exclusion_policy": self.component_exclusion_policy,
                "radii_version": self.radii_version,
                "grid_phase": _GRID_PHASE,
                "grid_sampling": _GRID_SAMPLING,
            }
        )

    @property
    def compatibility_signature(self) -> tuple[tuple[str, object], ...]:
        return tuple(self.compatibility_parameters.items())

    def to_json(self) -> dict[str, object]:
        """Return all settings, including bounded execution controls."""

        return {
            "coarse_grid_spacing_angstrom": self.coarse_grid_spacing_angstrom,
            "fine_grid_spacing_angstrom": self.fine_grid_spacing_angstrom,
            "boundary_margin_angstrom": self.boundary_margin_angstrom,
            "component_exclusion_policy": self.component_exclusion_policy,
            "max_voxel_count": self.max_voxel_count,
            "voxel_chunk_size": self.voxel_chunk_size,
            "radii_version": self.radii_version,
            "grid_phase": _GRID_PHASE,
            "grid_sampling": _GRID_SAMPLING,
        }


@dataclass(frozen=True, slots=True)
class PocketVolumeSensitivity:
    """Absolute and fine-grid-relative dual-resolution sensitivity."""

    absolute_angstrom3: float | None = None
    relative_fraction: float | None = None

    def __post_init__(self) -> None:
        for name in ("absolute_angstrom3", "relative_fraction"):
            value = getattr(self, name)
            if value is not None:
                value = _finite(value, name, non_negative=True)
                object.__setattr__(self, name, value)
        if self.relative_fraction is not None and self.absolute_angstrom3 is None:
            raise ValueError("relative_fraction requires an absolute sensitivity")

    @property
    def absolute_sensitivity_angstrom3(self) -> float | None:
        return self.absolute_angstrom3

    @property
    def relative_sensitivity(self) -> float | None:
        return self.relative_fraction

    def to_json(self) -> dict[str, float | None]:
        return {
            "absolute_angstrom3": self.absolute_angstrom3,
            "relative_fraction": self.relative_fraction,
        }


# A descriptive alias retained for integrations that use the longer name.
PocketVolumeSensitivityResult = PocketVolumeSensitivity


@dataclass(frozen=True, slots=True)
class PocketVolumeResult:
    """Immutable result of one bounded dual-grid volume measurement."""

    availability: Availability
    coarse_voxel_count: int | None = None
    fine_voxel_count: int | None = None
    coarse_volume_angstrom3: float | None = None
    fine_volume_angstrom3: float | None = None
    grid_shape: tuple[int, int, int] | None = None
    grid_origin_xyz: tuple[float, float, float] | None = None
    grid_phase: str = _GRID_PHASE
    coarse_grid_shape: tuple[int, int, int] | None = None
    fine_grid_shape: tuple[int, int, int] | None = None
    rotation_error_bound_angstrom3: float | None = None
    sensitivity: PocketVolumeSensitivity | None = None
    units: Mapping[str, object] = field(default_factory=lambda: _VOLUME_UNITS)
    provenance_parameters: Mapping[str, object] = field(default_factory=dict)
    diagnostics: tuple[Diagnostic, ...] = ()
    candidate_id: str | None = None
    compatibility_signature: tuple[object, ...] = ()
    settings: PocketVolumeSettings | None = field(default=None, compare=False, repr=False)
    provenance: MethodProvenance | None = None

    def __post_init__(self) -> None:
        availability = self.availability
        if not isinstance(availability, Availability):
            try:
                availability = Availability(availability)
            except (TypeError, ValueError) as exc:
                raise ValueError(f"unknown pocket volume availability: {self.availability!r}") from exc
            object.__setattr__(self, "availability", availability)
        for name in ("coarse_voxel_count", "fine_voxel_count"):
            value = getattr(self, name)
            if value is not None:
                object.__setattr__(self, name, _positive_or_zero_int(value, name))
        for name in ("coarse_volume_angstrom3", "fine_volume_angstrom3", "rotation_error_bound_angstrom3"):
            value = getattr(self, name)
            if value is not None:
                object.__setattr__(self, name, _finite(value, name, non_negative=True))
        if self.grid_shape is not None:
            object.__setattr__(self, "grid_shape", _grid_shape(self.grid_shape))
        if self.coarse_grid_shape is not None:
            object.__setattr__(self, "coarse_grid_shape", _grid_shape(self.coarse_grid_shape))
        if self.fine_grid_shape is not None:
            object.__setattr__(self, "fine_grid_shape", _grid_shape(self.fine_grid_shape))
        if self.grid_origin_xyz is not None:
            raw_origin = tuple(self.grid_origin_xyz)
            if len(raw_origin) != 3:
                raise ValueError("grid_origin_xyz must contain three finite values")
            origin = (float(raw_origin[0]), float(raw_origin[1]), float(raw_origin[2]))
            if any(not math.isfinite(value) for value in origin):
                raise ValueError("grid_origin_xyz must contain three finite values")
            object.__setattr__(self, "grid_origin_xyz", origin)
        if self.grid_phase != _GRID_PHASE:
            raise ValueError("grid_phase must be 'cell_center'")
        sensitivity = self.sensitivity
        if sensitivity is not None and not isinstance(sensitivity, PocketVolumeSensitivity):
            raise TypeError("sensitivity must be PocketVolumeSensitivity or None")
        settings = self.settings
        if settings is not None and not isinstance(settings, PocketVolumeSettings):
            raise TypeError("settings must be PocketVolumeSettings or None")
        object.__setattr__(self, "units", _freeze_json(self.units))
        object.__setattr__(self, "provenance_parameters", _freeze_json(self.provenance_parameters))
        if self.provenance is not None and not isinstance(self.provenance, MethodProvenance):
            raise TypeError("provenance must be a MethodProvenance or None")
        diagnostics = tuple(self.diagnostics)
        if any(not isinstance(item, Diagnostic) for item in diagnostics):
            raise TypeError("diagnostics must contain Diagnostic values")
        object.__setattr__(self, "diagnostics", diagnostics)
        object.__setattr__(self, "compatibility_signature", tuple(self.compatibility_signature))
        if self.candidate_id is not None:
            candidate_id = str(self.candidate_id).strip()
            if not candidate_id:
                raise ValueError("candidate_id must not be empty when provided")
            object.__setattr__(self, "candidate_id", candidate_id)

    @property
    def status(self) -> Availability:
        return self.availability

    @property
    def absolute_sensitivity_angstrom3(self) -> float | None:
        return self.sensitivity.absolute_angstrom3 if self.sensitivity is not None else None

    @property
    def relative_sensitivity(self) -> float | None:
        return self.sensitivity.relative_fraction if self.sensitivity is not None else None

    def to_json(self) -> dict[str, object]:
        settings = self.settings.to_json() if self.settings is not None else None
        return {
            "availability": self.availability.value,
            "coarse_voxel_count": self.coarse_voxel_count,
            "fine_voxel_count": self.fine_voxel_count,
            "coarse_volume_angstrom3": self.coarse_volume_angstrom3,
            "fine_volume_angstrom3": self.fine_volume_angstrom3,
            "coarse": {
                "voxel_count": self.coarse_voxel_count,
                "volume_angstrom3": self.coarse_volume_angstrom3,
                "grid_shape": list(self.coarse_grid_shape or self.grid_shape or ()),
            },
            "fine": {
                "voxel_count": self.fine_voxel_count,
                "volume_angstrom3": self.fine_volume_angstrom3,
                "grid_shape": list(self.fine_grid_shape or ()),
            },
            "sensitivity": self.sensitivity.to_json() if self.sensitivity is not None else None,
            "absolute_sensitivity_angstrom3": self.absolute_sensitivity_angstrom3,
            "relative_sensitivity": self.relative_sensitivity,
            "rotation_error_bound_angstrom3": self.rotation_error_bound_angstrom3,
            "grid": {
                "shape": list(self.grid_shape) if self.grid_shape is not None else None,
                "coarse_shape": list(self.coarse_grid_shape) if self.coarse_grid_shape is not None else None,
                "fine_shape": list(self.fine_grid_shape) if self.fine_grid_shape is not None else None,
                "origin_xyz": list(self.grid_origin_xyz) if self.grid_origin_xyz is not None else None,
                "phase": self.grid_phase,
                "sampling": _GRID_SAMPLING,
            },
            "units": _thaw_json(self.units),
            "provenance_parameters": _thaw_json(self.provenance_parameters),
            "provenance": self.provenance.to_json() if self.provenance is not None else None,
            "settings": settings,
            "candidate_id": self.candidate_id,
            "compatibility_signature": _thaw_json(self.compatibility_signature),
            "diagnostics": [item.to_json() for item in self.diagnostics],
        }


@dataclass(frozen=True, slots=True)
class PocketVolumeComparison:
    """Typed target-minus-reference comparison of compatible measurements."""

    availability: Availability
    delta_angstrom3: float | None = None
    relative_delta_fraction: float | None = None
    reference_volume_angstrom3: float | None = None
    target_volume_angstrom3: float | None = None
    diagnostics: tuple[Diagnostic, ...] = ()
    units: Mapping[str, object] = field(default_factory=lambda: _VOLUME_UNITS)
    compatibility_signature: tuple[object, ...] = ()

    def __post_init__(self) -> None:
        availability = self.availability
        if not isinstance(availability, Availability):
            try:
                availability = Availability(availability)
            except (TypeError, ValueError) as exc:
                raise ValueError(f"unknown pocket volume comparison availability: {self.availability!r}") from exc
            object.__setattr__(self, "availability", availability)
        for name in (
            "delta_angstrom3",
            "relative_delta_fraction",
            "reference_volume_angstrom3",
            "target_volume_angstrom3",
        ):
            value = getattr(self, name)
            if value is not None:
                value = _finite(value, name)
                object.__setattr__(self, name, value)
        object.__setattr__(self, "units", _freeze_json(self.units))
        diagnostics = tuple(self.diagnostics)
        if any(not isinstance(item, Diagnostic) for item in diagnostics):
            raise TypeError("diagnostics must contain Diagnostic values")
        object.__setattr__(self, "diagnostics", diagnostics)
        object.__setattr__(self, "compatibility_signature", tuple(self.compatibility_signature))

    @property
    def status(self) -> Availability:
        return self.availability

    @property
    def relative_delta(self) -> float | None:
        return self.relative_delta_fraction

    def to_json(self) -> dict[str, object]:
        return {
            "availability": self.availability.value,
            "delta_angstrom3": self.delta_angstrom3,
            "relative_delta_fraction": self.relative_delta_fraction,
            "relative_delta": self.relative_delta_fraction,
            "reference_volume_angstrom3": self.reference_volume_angstrom3,
            "target_volume_angstrom3": self.target_volume_angstrom3,
            "units": _thaw_json(self.units),
            "compatibility_signature": _thaw_json(self.compatibility_signature),
            "diagnostics": [item.to_json() for item in self.diagnostics],
        }


def _positive_or_zero_int(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a non-negative integer")
    return value


def _grid_shape(shape: Sequence[int]) -> tuple[int, int, int]:
    values = tuple(shape)
    if len(values) != 3 or any(isinstance(value, bool) or not isinstance(value, int) or value <= 0 for value in values):
        raise ValueError("grid shape must contain three positive integers")
    return (values[0], values[1], values[2])


__all__ = [
    "PocketVolumeComparison",
    "PocketVolumeResult",
    "PocketVolumeSensitivity",
    "PocketVolumeSensitivityResult",
    "PocketVolumeSettings",
]
