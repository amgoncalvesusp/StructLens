"""Immutable, deterministic provenance for scientific artifacts.

Wall-clock audit information is intentionally represented by :class:`AuditEvent`
and kept outside the canonical scientific payload.  This makes rerunning the
same method on the same inputs byte-identical while still allowing a workflow
to record when and where an action occurred.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from collections.abc import Mapping
from dataclasses import dataclass, field
from datetime import UTC, datetime
from types import MappingProxyType
from typing import TypeAlias, cast

JSONScalar: TypeAlias = str | int | float | bool | None
FrozenJSON: TypeAlias = JSONScalar | tuple["FrozenJSON", ...] | Mapping[str, "FrozenJSON"]

_SHA256 = re.compile(r"^[0-9a-fA-F]{64}$")
_PHYSICAL_TOKENS = (
    "angstrom",
    "nanometer",
    "distance",
    "radius",
    "spacing",
    "volume",
    "area",
    "length",
    "coordinate",
    "centroid",
    "rmsd",
    "sasa",
    "probe",
)
_SUPPORTED_UNITS = {
    "angstrom",
    "angstroms",
    "å",
    "a",
    "angstrom^2",
    "angstrom2",
    "angstrom²",
    "a2",
    "angstrom^3",
    "angstrom3",
    "angstrom³",
    "a3",
    "degree",
    "degrees",
    "deg",
    "radian",
    "radians",
    "rad",
    "count",
    "counts",
    "dimensionless",
    "fraction",
    "percent",
    "%",
}
MAX_UNIT_MAP_ITEMS = 256
MAX_UNIT_STRING_LENGTH = 1_024


def _freeze(value: object, *, path: str) -> FrozenJSON:
    """Recursively normalize JSON values into immutable values.

    Lists are copied to tuples and mappings to read-only mapping proxies.  Sets,
    bytes, and arbitrary objects are rejected because they do not have a stable
    JSON representation suitable for a scientific artifact ID.
    """

    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError(f"{path} must contain only finite numbers")
        return value
    if isinstance(value, Mapping):
        frozen: dict[str, FrozenJSON] = {}
        for key, item in value.items():
            if not isinstance(key, str) or not key:
                raise TypeError(f"{path} mapping keys must be non-empty strings")
            frozen[key] = _freeze(item, path=f"{path}.{key}")
        return cast(FrozenJSON, MappingProxyType(frozen))
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(item, path=f"{path}[{index}]") for index, item in enumerate(value))
    raise TypeError(f"{path} contains a mutable or non-JSON value: {type(value).__name__}")


def _json_ready(value: FrozenJSON) -> JSONScalar | list[object] | dict[str, object]:
    if isinstance(value, Mapping):
        return {key: _json_ready(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_json_ready(item) for item in value]
    return value


def freeze_json(value: object, *, path: str = "value") -> FrozenJSON:
    """Public boundary helper returning recursively immutable JSON data."""

    return _freeze(value, path=path)


def json_ready(value: FrozenJSON) -> JSONScalar | list[object] | dict[str, object]:
    """Return a fresh mutable JSON-compatible copy of frozen data."""

    return _json_ready(value)


def _canonical_bytes(value: Mapping[str, object]) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")


def _validate_text(value: str, field_name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{field_name} must be a non-empty string")
    if any(char.isspace() for char in value):
        raise ValueError(f"{field_name} must not contain whitespace")
    return value


def freeze_bounded_string_map(
    values: Mapping[str, str],
    *,
    field_name: str,
    maximum_items: int = MAX_UNIT_MAP_ITEMS,
    maximum_string_length: int = MAX_UNIT_STRING_LENGTH,
) -> Mapping[str, str]:
    """Validate and snapshot a small public string map before materialization."""

    if not isinstance(values, Mapping):
        raise TypeError(f"{field_name} must be a mapping")
    if len(values) > maximum_items:
        raise ValueError(f"{field_name} exceed the maximum item limit")
    frozen: dict[str, str] = {}
    for key, value in values.items():
        if not isinstance(key, str) or not key.strip() or not isinstance(value, str) or not value.strip():
            raise ValueError(f"{field_name} must map non-empty names to non-empty strings")
        if len(key) > maximum_string_length or len(value) > maximum_string_length:
            raise ValueError(f"{field_name} exceed the maximum string length")
        frozen[key] = value
    return MappingProxyType(frozen)


def _freeze_string_map(values: Mapping[str, str], field_name: str) -> Mapping[str, str]:
    return freeze_bounded_string_map(values, field_name=field_name)


def _freeze_units(values: Mapping[str, str]) -> Mapping[str, str]:
    frozen = _freeze_string_map(values, "units")
    unsupported = [f"{name}={unit!r}" for name, unit in frozen.items() if unit.casefold() not in _SUPPORTED_UNITS]
    if unsupported:
        raise ValueError(f"unsupported unit(s): {', '.join(sorted(unsupported))}")
    return frozen


def _parameter_paths(value: FrozenJSON, path: str = "") -> tuple[tuple[str, float | int], ...]:
    if isinstance(value, Mapping):
        entries: list[tuple[str, float | int]] = []
        for key, child in value.items():
            child_path = f"{path}.{key}" if path else key
            entries.extend(_parameter_paths(child, child_path))
        return tuple(entries)
    if isinstance(value, tuple):
        entries = []
        for index, child in enumerate(value):
            entries.extend(_parameter_paths(child, f"{path}[{index}]"))
        return tuple(entries)
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return ((path, value),)
    return ()


def _requires_unit(path: str) -> bool:
    leaf = re.sub(r"\[\d+\]", "", path).rsplit(".", 1)[-1].casefold()
    if leaf == "count" or leaf.endswith("_count") or leaf.endswith("counts"):
        return False
    return _dimension_for(path) is not None


def _dimension_for(path: str) -> str | None:
    leaf = re.sub(r"\[\d+\]", "", path).rsplit(".", 1)[-1].casefold()
    if leaf == "count" or leaf.endswith("_count") or leaf.endswith("counts"):
        return "dimensionless"
    if "volume" in leaf:
        return "volume"
    if "area" in leaf or "sasa" in leaf:
        return "area"
    if "angle" in leaf:
        return "angle"
    if any(token in leaf for token in _PHYSICAL_TOKENS):
        return "length"
    return None


def _unit_dimension(unit: str) -> str | None:
    normalized = unit.casefold()
    if normalized in {"angstrom", "angstroms", "å", "a"}:
        return "length"
    if normalized in {"angstrom^2", "angstrom2", "angstrom²", "a2"}:
        return "area"
    if normalized in {"angstrom^3", "angstrom3", "angstrom³", "a3"}:
        return "volume"
    if normalized in {"degree", "degrees", "deg", "radian", "radians", "rad"}:
        return "angle"
    if normalized in {"count", "counts", "dimensionless", "fraction", "percent", "%"}:
        return "dimensionless"
    return None


def _unit_key(path: str) -> str:
    return re.sub(r"\[\d+\]", "", path)


@dataclass(frozen=True, slots=True)
class AuditEvent:
    """Workflow metadata kept outside canonical scientific provenance."""

    event_type: str
    occurred_at: datetime
    user: str | None = None
    workstation: str | None = None
    context: Mapping[str, FrozenJSON] = field(default_factory=dict)

    def __post_init__(self) -> None:
        _validate_text(self.event_type, "event_type")
        if not isinstance(self.occurred_at, datetime):
            raise ValueError("occurred_at must be a datetime")
        if self.occurred_at.tzinfo is None or self.occurred_at.utcoffset() is None:
            raise ValueError("occurred_at must be timezone-aware")
        for field_name in ("user", "workstation"):
            value = getattr(self, field_name)
            if value is not None:
                _validate_text(value, field_name)
        context = _freeze(self.context, path="context")
        if not isinstance(context, Mapping):
            raise TypeError("context must be a mapping")
        object.__setattr__(self, "context", context)

    def to_json(self) -> dict[str, object]:
        timestamp = self.occurred_at.astimezone(UTC).isoformat().replace("+00:00", "Z")
        return {
            "event_type": self.event_type,
            "occurred_at": timestamp,
            "user": self.user,
            "workstation": self.workstation,
            "context": _json_ready(self.context),
        }


@dataclass(frozen=True, slots=True)
class MethodProvenance:
    """Validated method identity, parameters, inputs, and deterministic ID."""

    method_id: str
    method_version: str
    parameters: Mapping[str, FrozenJSON] = field(default_factory=dict)
    units: Mapping[str, str] = field(default_factory=dict)
    backend_versions: Mapping[str, str] = field(default_factory=dict)
    input_hashes: Mapping[str, str] = field(default_factory=dict)
    analyzed_representation: str = "as_file"
    artifact_id: str = field(init=False)

    def __post_init__(self) -> None:
        _validate_text(self.method_id, "method_id")
        _validate_text(self.method_version, "method_version")
        _validate_text(self.analyzed_representation, "analyzed_representation")

        parameters = _freeze(self.parameters, path="parameters")
        if not isinstance(parameters, Mapping):
            raise TypeError("parameters must be a mapping")
        units = _freeze_units(self.units)
        backend_versions = _freeze_string_map(self.backend_versions, "backend_versions")
        input_hashes = _freeze_string_map(self.input_hashes, "input_hashes")
        if not input_hashes:
            raise ValueError("at least one input hash is required")
        normalized_sources: dict[str, str] = {}
        for source, digest in input_hashes.items():
            if not _SHA256.fullmatch(digest):
                raise ValueError(f"input hash for {source!r} must be a SHA-256 hexadecimal digest")
            logical_source = source.casefold().replace("-", "_")
            if logical_source in normalized_sources:
                raise ValueError(f"inconsistent duplicate source hashes for {source!r}")
            normalized_sources[logical_source] = digest.lower()
        input_hashes = MappingProxyType({source: digest.lower() for source, digest in input_hashes.items()})

        missing_units = [
            path for path, _ in _parameter_paths(parameters) if _requires_unit(path) and _unit_key(path) not in units
        ]
        if missing_units:
            joined = ", ".join(sorted(missing_units))
            raise ValueError(f"physical parameter(s) {joined} require a unit")
        incompatible_units = [
            (path, units[_unit_key(path)], _dimension_for(path))
            for path, _ in _parameter_paths(parameters)
            if _dimension_for(path) is not None
            and _unit_key(path) in units
            and _unit_dimension(units[_unit_key(path)]) != _dimension_for(path)
        ]
        if incompatible_units:
            details = ", ".join(
                f"{path}={unit!r} (expected {dimension})" for path, unit, dimension in incompatible_units
            )
            raise ValueError(f"incompatible unit(s): {details}")

        object.__setattr__(self, "parameters", parameters)
        object.__setattr__(self, "units", units)
        object.__setattr__(self, "backend_versions", backend_versions)
        object.__setattr__(self, "input_hashes", input_hashes)
        scientific = self._scientific_payload(include_artifact_id=False)
        artifact_id = hashlib.sha256(_canonical_bytes(scientific)).hexdigest()
        object.__setattr__(self, "artifact_id", artifact_id)

    @classmethod
    def from_json(cls, payload: Mapping[str, object] | str) -> MethodProvenance:
        """Restore and verify a canonical provenance payload.

        The stored artifact ID is never trusted: it is recomputed by the normal
        constructor and compared with the serialized value.
        """

        decoded: object = json.loads(payload) if isinstance(payload, str) else payload
        if not isinstance(decoded, Mapping):
            raise TypeError("provenance JSON must contain an object")
        stored_id = decoded.get("artifact_id")
        if not isinstance(stored_id, str) or not _SHA256.fullmatch(stored_id):
            raise ValueError("artifact_id must be a SHA-256 hexadecimal digest")
        provenance = cls(
            method_id=cast(str, decoded.get("method_id")),
            method_version=cast(str, decoded.get("method_version")),
            parameters=cast(Mapping[str, FrozenJSON], decoded.get("parameters", {})),
            units=cast(Mapping[str, str], decoded.get("units", {})),
            backend_versions=cast(Mapping[str, str], decoded.get("backend_versions", {})),
            input_hashes=cast(Mapping[str, str], decoded.get("input_hashes", {})),
            analyzed_representation=cast(str, decoded.get("analyzed_representation")),
        )
        if provenance.artifact_id != stored_id.lower():
            raise ValueError("artifact_id does not match the canonical provenance content")
        return provenance

    def _scientific_payload(self, *, include_artifact_id: bool) -> dict[str, object]:
        payload: dict[str, object] = {
            "method_id": self.method_id,
            "method_version": self.method_version,
            "parameters": _json_ready(self.parameters),
            "units": dict(self.units),
            "backend_versions": dict(self.backend_versions),
            "input_hashes": dict(self.input_hashes),
            "analyzed_representation": self.analyzed_representation,
        }
        if include_artifact_id:
            payload["artifact_id"] = self.artifact_id
        return payload

    def to_json(self) -> dict[str, object]:
        """Return the canonical scientific payload as fresh JSON-ready data."""

        return self._scientific_payload(include_artifact_id=True)

    def canonical_bytes(self) -> bytes:
        """Serialize scientific content deterministically, excluding audit time."""

        return _canonical_bytes(self.to_json())

    def with_audit_event(self, event: AuditEvent) -> dict[str, object]:
        """Wrap canonical content with non-canonical workflow metadata."""

        if not isinstance(event, AuditEvent):
            raise TypeError("event must be an AuditEvent")
        return {"scientific": self.to_json(), "audit_event": event.to_json()}

    def compatibility_view(self) -> Mapping[str, str]:
        """Expose a flat ``Mapping[str, str]`` for v0.3 consumers."""

        values = {
            "method_id": self.method_id,
            "method_version": self.method_version,
            "analyzed_representation": self.analyzed_representation,
            "artifact_id": self.artifact_id,
        }
        values.update({f"input_hash:{source}": digest for source, digest in self.input_hashes.items()})
        return MappingProxyType(values)

    # A descriptive alias for callers migrating from legacy provenance maps.
    as_legacy_mapping = compatibility_view

    def __hash__(self) -> int:
        return hash(self.artifact_id)


__all__ = [
    "AuditEvent",
    "FrozenJSON",
    "JSONScalar",
    "MAX_UNIT_MAP_ITEMS",
    "MAX_UNIT_STRING_LENGTH",
    "MethodProvenance",
    "freeze_bounded_string_map",
    "freeze_json",
    "json_ready",
]
