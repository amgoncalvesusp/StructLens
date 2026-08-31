"""Immutable contracts and validation helpers for pocket free-volume runs."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any, cast

from structlens.core.evidence import Availability, Diagnostic, DiagnosticSeverity
from structlens.core.provenance import MethodProvenance, freeze_bounded_string_map

from .radii import POCKET_RADII_VERSION

_GRID_PHASE = "cell_center"
_GRID_SAMPLING = "origin + (index + 0.5) * spacing"
_VOLUME_UNITS = MappingProxyType({"volume": "angstrom^3", "length": "angstrom"})
_COMPONENT_POLICIES = frozenset({"protein_only", "unoccupied"})
_MAX_VOLUME_DIAGNOSTICS = 100_000
_CORE_VOLUME_PRODUCER = ("structlens.pocket_free_volume", "0.4")
_APPLICATION_VOLUME_PRODUCER = ("structlens.pocket.volume", "0.4.0")
_APPLICATION_SCOPE_KEYS = ("altloc_policy", "hydrogen_policy")
_APPLICATION_RESULT_SOURCE_KEYS = ("raw_source_hash", "logical_content")
_ALLOWED_ALTLOC_POLICIES = frozenset({"highest_occupancy", "first", "all"})
_ALLOWED_ASSEMBLY_SCOPES = frozenset({"asymmetric_unit", "biological_assembly"})
_ALLOWED_STRUCTURE_FORMATS = frozenset({"pdb", "mmcif"})
_VOLUME_HYDROGEN_POLICY = "deposited_heavy_atoms_only"
_VOLUME_ATOM_SCOPE = "selected_primary_polymer_heavy_atoms"
_VOLUME_COMPONENT_SCOPE = "selected_retained_ligand_ion_other_components"
_VOLUME_COMPONENT_RULES_VERSION = "structlens-pocket-known-ligands-1"
_METHOD_PARAMETER_KEYS = (
    "coarse_grid_spacing_angstrom",
    "fine_grid_spacing_angstrom",
    "boundary_margin_angstrom",
    "component_exclusion_policy",
    "radii_version",
    "grid_phase",
    "grid_sampling",
)
_CANDIDATE_PARAMETER_KEYS = (
    "candidate_id",
    "source_content_id",
    "selection_id",
    "candidate_lineage",
    "sphere_count",
    "sphere_ids",
    "sphere_radii_angstrom",
)
_RESULT_OPTIONAL_PARAMETER_KEYS = frozenset({"active_exclusion_atom_count", "rotation_error_bound"})
_APPLICATION_PARAMETER_KEYS = frozenset(_CANDIDATE_PARAMETER_KEYS) | frozenset(
    {
        "selection",
        "model_id",
        "author_chain_ids",
        "label_chain_ids",
        "altloc_policy",
        "assembly_scope",
        "atom_scope",
        "component_scope",
        "component_rules_version",
        "hydrogen_policy",
        "polymer_atom_count",
        "component_ids",
    }
)
_APPLICATION_SELECTION_KEYS = frozenset(
    "selection_id format model_id author_chain_ids label_chain_ids chain_locators altloc_policy assembly_scope".split()
)
_METHOD_UNITS = {
    "volume": "angstrom^3",
    "length": "angstrom",
    "coarse_grid_spacing_angstrom": "angstrom",
    "fine_grid_spacing_angstrom": "angstrom",
    "boundary_margin_angstrom": "angstrom",
    "sphere_radii_angstrom": "angstrom",
}


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


def _diagnostic(
    code: str,
    message: str,
    *,
    severity: DiagnosticSeverity | None = None,
    source_id: str | None = None,
    atom_id: str | None = None,
    remediation: str | None = None,
) -> Diagnostic:
    return Diagnostic(
        code=code,
        severity=(
            severity
            if severity is not None
            else DiagnosticSeverity.ERROR
            if code.endswith("resource_limit")
            else DiagnosticSeverity.WARNING
        ),
        message=message,
        source_id=source_id,
        atom_id=atom_id,
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
    max_sphere_count: int = 10_000
    max_exclusion_atom_count: int = 50_000
    max_sphere_voxel_checks: int = 100_000_000
    max_exclusion_neighbor_checks: int = 5_000_000
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
            raise ValueError("component_exclusion_policy must be one of: protein_only, unoccupied")
        object.__setattr__(self, "component_exclusion_policy", policy)
        object.__setattr__(self, "max_voxel_count", _positive_int(self.max_voxel_count, "max_voxel_count"))
        object.__setattr__(self, "voxel_chunk_size", _positive_int(self.voxel_chunk_size, "voxel_chunk_size"))
        for name in (
            "max_sphere_count",
            "max_exclusion_atom_count",
            "max_sphere_voxel_checks",
            "max_exclusion_neighbor_checks",
        ):
            object.__setattr__(self, name, _positive_int(getattr(self, name), name))
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
            "max_sphere_count": self.max_sphere_count,
            "max_exclusion_atom_count": self.max_exclusion_atom_count,
            "max_sphere_voxel_checks": self.max_sphere_voxel_checks,
            "max_exclusion_neighbor_checks": self.max_exclusion_neighbor_checks,
            "radii_version": self.radii_version,
            "grid_phase": _GRID_PHASE,
            "grid_sampling": _GRID_SAMPLING,
        }


def _validate_schema_keys(
    parameters: Mapping[str, object],
    required: set[str] | frozenset[str],
    allowed: set[str] | frozenset[str],
    label: str,
) -> None:
    keys = set(parameters)
    missing = required - keys
    unknown = keys - allowed
    if missing or unknown:
        raise ValueError(
            f"{label} schema is incomplete or contains unknown fields "
            f"(missing={sorted(missing)}, unknown={sorted(unknown)})"
        )


def _bounded_sequence(value: object, name: str, maximum: int) -> tuple[object, ...]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{name} must be a bounded sequence")
    if len(value) > maximum:
        raise ValueError(f"{name} exceeds the declared sphere limit")
    return tuple(value)


def _validate_candidate_parameters(
    value: PocketVolumeResult,
    result_parameters: Mapping[str, object],
    method_parameters: Mapping[str, object],
    settings: PocketVolumeSettings,
) -> tuple[str | None, str | None]:
    for name in _CANDIDATE_PARAMETER_KEYS:
        if result_parameters[name] != method_parameters[name]:
            raise ValueError(f"pocket volume candidate provenance field {name} is incoherent")
    candidate_id = result_parameters["candidate_id"]
    if not isinstance(candidate_id, str) or candidate_id != value.candidate_id:
        raise ValueError("pocket volume candidate_id provenance is incoherent")
    source_content_id = result_parameters["source_content_id"]
    selection_id = result_parameters["selection_id"]
    if source_content_id is not None and not isinstance(source_content_id, str):
        raise ValueError("pocket volume source_content_id provenance is invalid")
    if selection_id is not None and not isinstance(selection_id, str):
        raise ValueError("pocket volume selection_id provenance is invalid")
    if (source_content_id is None) != (selection_id is None):
        raise ValueError("pocket volume source and selection provenance must be declared together")
    expected_lineage = {
        "candidate_id": candidate_id,
        "source_content_id": source_content_id,
        "selection_id": selection_id,
    }
    lineage = result_parameters["candidate_lineage"]
    if not isinstance(lineage, Mapping) or dict(lineage) != expected_lineage:
        raise ValueError("pocket volume candidate_lineage provenance is incoherent")
    sphere_ids = _bounded_sequence(result_parameters["sphere_ids"], "sphere_ids", settings.max_sphere_count)
    sphere_radii = _bounded_sequence(
        result_parameters["sphere_radii_angstrom"],
        "sphere_radii_angstrom",
        settings.max_sphere_count,
    )
    sphere_count = result_parameters["sphere_count"]
    if isinstance(sphere_count, bool) or not isinstance(sphere_count, int) or sphere_count <= 0:
        raise ValueError("sphere_count must be a positive integer")
    if sphere_count != len(sphere_ids) or sphere_count != len(sphere_radii):
        raise ValueError("sphere_count is incoherent with sphere IDs and radii")
    if any(not isinstance(item, str) or not item for item in sphere_ids) or len(set(sphere_ids)) != len(sphere_ids):
        raise ValueError("sphere_ids must contain unique non-empty strings")
    if any(_finite(item, "sphere_radii_angstrom", positive=True) <= 0.0 for item in sphere_radii):
        raise ValueError("sphere_radii_angstrom must contain positive finite values")
    return source_content_id, selection_id


def _validate_core_schema(
    value: PocketVolumeResult,
    provenance: MethodProvenance,
    result_parameters: Mapping[str, object],
    method_parameters: Mapping[str, object],
    source_content_id: str | None,
) -> None:
    required = set(_METHOD_PARAMETER_KEYS) | set(_CANDIDATE_PARAMETER_KEYS)
    _validate_schema_keys(method_parameters, required, required, "core pocket volume provenance")
    _validate_schema_keys(
        result_parameters,
        required,
        required | _RESULT_OPTIONAL_PARAMETER_KEYS,
        "core pocket volume result provenance",
    )
    expected_hashes = {
        "candidate": value.candidate_id,
        "source_content": source_content_id or value.candidate_id,
    }
    if dict(provenance.input_hashes) != expected_hashes:
        raise ValueError("core pocket volume input hashes are incoherent with candidate provenance")
    if provenance.backend_versions or provenance.analyzed_representation != "as_file":
        raise ValueError("core pocket volume producer schema is incoherent")
    expected_units = {name: unit for name, unit in _METHOD_UNITS.items() if name != "length"}
    if dict(provenance.units) != expected_units:
        raise ValueError("core pocket volume provenance units are incoherent")


def _validate_text_sequence(value: object, name: str) -> tuple[object, ...]:
    values = _bounded_sequence(value, name, _MAX_VOLUME_DIAGNOSTICS)
    if any(not isinstance(item, str) or not item.strip() for item in values):
        raise ValueError(f"{name} must contain non-empty strings")
    return values


def _validate_application_selection(parameters: Mapping[str, object], provenance: MethodProvenance) -> None:
    selection = parameters["selection"]
    if not isinstance(selection, Mapping):
        raise ValueError("application pocket volume selection provenance must be a mapping")
    _validate_schema_keys(
        selection,
        _APPLICATION_SELECTION_KEYS,
        _APPLICATION_SELECTION_KEYS,
        "application pocket volume selection",
    )
    for name in ("selection_id", "model_id", "author_chain_ids", "label_chain_ids", "altloc_policy", "assembly_scope"):
        if selection[name] != parameters[name]:
            raise ValueError(f"application pocket volume selection field {name} is incoherent")
    if selection["format"] not in _ALLOWED_STRUCTURE_FORMATS:
        raise ValueError("application pocket volume selection format is not recognized")
    if parameters["altloc_policy"] not in _ALLOWED_ALTLOC_POLICIES:
        raise ValueError("application pocket volume altloc policy is not recognized")
    if parameters["assembly_scope"] not in _ALLOWED_ASSEMBLY_SCOPES:
        raise ValueError("application pocket volume assembly policy is not recognized")
    if provenance.analyzed_representation != parameters["assembly_scope"]:
        raise ValueError("application pocket volume analyzed representation is incoherent")
    _validate_text_sequence(parameters["author_chain_ids"], "author_chain_ids")
    _validate_text_sequence(parameters["label_chain_ids"], "label_chain_ids")
    chain_locators = _bounded_sequence(selection["chain_locators"], "chain_locators", _MAX_VOLUME_DIAGNOSTICS)
    expected_locator_keys = {"author_chain_id", "label_chain_id", "entity_id"}
    if any(not isinstance(item, Mapping) or set(item) != expected_locator_keys for item in chain_locators):
        raise ValueError("application pocket volume chain locator schema is incoherent")


def _validate_application_schema(
    value: PocketVolumeResult,
    settings: PocketVolumeSettings,
    provenance: MethodProvenance,
    result_parameters: Mapping[str, object],
    method_parameters: Mapping[str, object],
    source_content_id: str | None,
    selection_id: str | None,
) -> None:
    if source_content_id is None or selection_id is None:
        raise ValueError("application pocket volume requires complete source/selection provenance")
    settings_payload = settings.to_json()
    method_required = set(settings_payload) | set(_APPLICATION_PARAMETER_KEYS)
    result_required = (
        set(_METHOD_PARAMETER_KEYS)
        | set(_CANDIDATE_PARAMETER_KEYS)
        | set(_APPLICATION_SCOPE_KEYS)
        | set(_APPLICATION_RESULT_SOURCE_KEYS)
    )
    _validate_schema_keys(method_parameters, method_required, method_required, "application pocket volume provenance")
    _validate_schema_keys(
        result_parameters,
        result_required,
        result_required | _RESULT_OPTIONAL_PARAMETER_KEYS,
        "application pocket volume result provenance",
    )
    for name, expected in settings_payload.items():
        if method_parameters[name] != expected:
            raise ValueError(f"application pocket volume setting {name} is incoherent")
    _validate_application_selection(method_parameters, provenance)
    fixed_policies = {
        "atom_scope": _VOLUME_ATOM_SCOPE,
        "component_scope": _VOLUME_COMPONENT_SCOPE,
        "component_rules_version": _VOLUME_COMPONENT_RULES_VERSION,
        "hydrogen_policy": _VOLUME_HYDROGEN_POLICY,
    }
    if any(method_parameters[name] != expected for name, expected in fixed_policies.items()):
        raise ValueError("application pocket volume policy provenance is incoherent")
    polymer_atom_count = method_parameters["polymer_atom_count"]
    if isinstance(polymer_atom_count, bool) or not isinstance(polymer_atom_count, int) or polymer_atom_count < 0:
        raise ValueError("application pocket volume polymer_atom_count is invalid")
    _validate_text_sequence(method_parameters["component_ids"], "component_ids")
    expected_hashes = {
        "raw_source",
        "logical_content",
        "candidate",
        "source_content",
    }
    if set(provenance.input_hashes) != expected_hashes:
        raise ValueError("application pocket volume input hash schema is incomplete")
    if (
        provenance.input_hashes["candidate"] != value.candidate_id
        or provenance.input_hashes["source_content"] != source_content_id
        or provenance.input_hashes["logical_content"] != source_content_id
        or result_parameters["raw_source_hash"] != provenance.input_hashes["raw_source"]
        or result_parameters["logical_content"] != provenance.input_hashes["logical_content"]
    ):
        raise ValueError("application pocket volume input hashes are incoherent with candidate provenance")
    if set(provenance.backend_versions) != {"numpy", "scipy", "pocket_radii"}:
        raise ValueError("application pocket volume backend schema is incomplete")
    if dict(provenance.units) != _METHOD_UNITS:
        raise ValueError("application pocket volume provenance units are incoherent")


def validate_pocket_volume_result(value: PocketVolumeResult) -> None:
    """Reject measurements that violate the strict schema of their producer."""

    if not isinstance(value, PocketVolumeResult):
        raise TypeError("value must be a PocketVolumeResult")
    settings = value.settings
    provenance = value.provenance
    if not isinstance(settings, PocketVolumeSettings) or not isinstance(provenance, MethodProvenance):
        raise ValueError("pocket volume result requires validated settings and method provenance")
    result_parameters = value.provenance_parameters
    method_parameters = provenance.parameters
    if not isinstance(result_parameters, Mapping) or not isinstance(method_parameters, Mapping):
        raise ValueError("pocket volume provenance parameters must be mappings")
    producer = (provenance.method_id, provenance.method_version)
    if producer not in {_CORE_VOLUME_PRODUCER, _APPLICATION_VOLUME_PRODUCER}:
        raise ValueError("pocket volume method producer is not recognized")
    if producer == _APPLICATION_VOLUME_PRODUCER and any(
        result_parameters.get(name) != method_parameters.get(name) for name in _APPLICATION_SCOPE_KEYS
    ):
        raise ValueError("application pocket volume scope provenance is incoherent")
    scope = tuple((name, result_parameters.get(name)) for name in _APPLICATION_SCOPE_KEYS)
    expected_signature = (("method", "pocket_free_volume"),) + settings.compatibility_signature
    if producer == _APPLICATION_VOLUME_PRODUCER:
        expected_signature += scope
    if value.compatibility_signature != expected_signature:
        raise ValueError("pocket volume compatibility signature is incoherent with producer settings")
    for label, parameters in (("result", result_parameters), ("method", method_parameters)):
        for name, expected in settings.compatibility_parameters.items():
            if parameters.get(name) != expected:
                raise ValueError(f"pocket volume {label} parameter {name} is incoherent with settings")
    if value.grid_phase != settings.compatibility_parameters["grid_phase"]:
        raise ValueError("pocket volume grid phase is incoherent with settings")
    if dict(value.units) != dict(_VOLUME_UNITS):
        raise ValueError("pocket volume result units are incoherent")
    source_content_id, selection_id = _validate_candidate_parameters(
        value,
        result_parameters,
        method_parameters,
        settings,
    )
    if producer == _CORE_VOLUME_PRODUCER:
        _validate_core_schema(value, provenance, result_parameters, method_parameters, source_content_id)
    else:
        _validate_application_schema(
            value,
            settings,
            provenance,
            result_parameters,
            method_parameters,
            source_content_id,
            selection_id,
        )


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
    units: Mapping[str, str] = field(default_factory=lambda: _VOLUME_UNITS)
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
        object.__setattr__(self, "units", freeze_bounded_string_map(self.units, field_name="units"))
        object.__setattr__(self, "provenance_parameters", _freeze_json(self.provenance_parameters))
        if self.provenance is not None and not isinstance(self.provenance, MethodProvenance):
            raise TypeError("provenance must be a MethodProvenance or None")
        if not isinstance(self.diagnostics, Sequence) or isinstance(self.diagnostics, (str, bytes)):
            raise TypeError("diagnostics must be a bounded Sequence")
        if len(self.diagnostics) > _MAX_VOLUME_DIAGNOSTICS:
            raise ValueError("diagnostics exceed the supported limit")
        diagnostics = tuple(self.diagnostics)
        if any(not isinstance(item, Diagnostic) for item in diagnostics):
            raise TypeError("diagnostics must contain Diagnostic values")
        object.__setattr__(self, "diagnostics", diagnostics)
        signature = tuple(self.compatibility_signature)
        if not signature and self.settings is not None:
            signature = (("method", "pocket_free_volume"),) + self.settings.compatibility_signature
        identity_tokens = {
            "candidate_id",
            "sphere_id",
            "sphere_ids",
            "sphere_count",
            "sphere_radii_angstrom",
            "grid_origin_xyz",
            "grid_shape",
        }
        for item in signature:
            if isinstance(item, (tuple, list)) and item and item[0] in identity_tokens:
                raise ValueError("compatibility_signature must not contain candidate or sphere identity")
        object.__setattr__(self, "compatibility_signature", signature)
        if self.candidate_id is not None:
            candidate_id = str(self.candidate_id).strip()
            if not candidate_id:
                raise ValueError("candidate_id must not be empty when provided")
            object.__setattr__(self, "candidate_id", candidate_id)
        quantitative = (
            self.coarse_voxel_count,
            self.fine_voxel_count,
            self.coarse_volume_angstrom3,
            self.fine_volume_angstrom3,
        )
        if availability is Availability.AVAILABLE:
            if any(value is None for value in quantitative) or self.sensitivity is None:
                raise ValueError("available pocket volume results require both grid counts, volumes, and sensitivity")
        elif (
            any(value is not None for value in quantitative)
            or self.sensitivity is not None
            or self.rotation_error_bound_angstrom3 is not None
        ):
            raise ValueError("unavailable pocket volume results must not carry quantitative measurements")

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
    units: Mapping[str, str] = field(default_factory=lambda: _VOLUME_UNITS)
    compatibility_signature: tuple[object, ...] = ()
    reference_sensitivity: PocketVolumeSensitivity | None = None
    target_sensitivity: PocketVolumeSensitivity | None = None
    provenance: tuple[MethodProvenance, ...] = ()

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
        object.__setattr__(self, "units", freeze_bounded_string_map(self.units, field_name="units"))
        if not isinstance(self.diagnostics, Sequence) or isinstance(self.diagnostics, (str, bytes)):
            raise TypeError("diagnostics must be a bounded Sequence")
        if len(self.diagnostics) > _MAX_VOLUME_DIAGNOSTICS:
            raise ValueError("diagnostics exceed the supported limit")
        diagnostics = tuple(self.diagnostics)
        if any(not isinstance(item, Diagnostic) for item in diagnostics):
            raise TypeError("diagnostics must contain Diagnostic values")
        object.__setattr__(self, "diagnostics", diagnostics)
        signature = tuple(self.compatibility_signature)
        for item in signature:
            if (
                isinstance(item, (tuple, list))
                and item
                and item[0]
                in {
                    "candidate_id",
                    "sphere_id",
                    "sphere_ids",
                    "sphere_count",
                    "sphere_radii_angstrom",
                    "grid_origin_xyz",
                    "grid_shape",
                }
            ):
                raise ValueError("compatibility_signature must not contain candidate or sphere identity")
        object.__setattr__(self, "compatibility_signature", signature)
        for name in ("reference_sensitivity", "target_sensitivity"):
            value = getattr(self, name)
            if value is not None and not isinstance(value, PocketVolumeSensitivity):
                raise TypeError(f"{name} must be PocketVolumeSensitivity or None")
        provenance = tuple(self.provenance)
        if any(not isinstance(item, MethodProvenance) for item in provenance):
            raise TypeError("provenance must contain MethodProvenance values")
        object.__setattr__(self, "provenance", provenance)

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
            "reference_sensitivity": self.reference_sensitivity.to_json()
            if self.reference_sensitivity is not None
            else None,
            "target_sensitivity": self.target_sensitivity.to_json() if self.target_sensitivity is not None else None,
            "provenance": [item.to_json() for item in self.provenance],
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
