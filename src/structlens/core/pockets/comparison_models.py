"""Immutable contracts used by pocket comparison.

The calculation orchestration stays in :mod:`comparison`; this module keeps
the serializable result models small and dependency-light.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any, cast

from structlens.core.evidence import Availability, Diagnostic
from structlens.core.interactions import InteractionDifference, InteractionRecord
from structlens.core.models import MutationEvent, ResidueId
from structlens.core.provenance import (
    FrozenJSON,
    MethodProvenance,
    freeze_bounded_string_map,
    freeze_json,
    json_ready,
)

from .matching import PocketMatch
from .volume_models import PocketVolumeComparison

MAX_COMPARISON_ITEMS = 100_000
MAX_SURFACE_PROVENANCE_ITEMS = 100_000
MAX_SURFACE_PROVENANCE_DEPTH = 64
MAX_SURFACE_PROVENANCE_STRING_LENGTH = 1_000_000


def _finite(value: object, name: str, *, non_negative: bool = False) -> float:
    try:
        numeric = float(cast(Any, value))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(numeric) or (non_negative and numeric < 0.0):
        suffix = " and non-negative" if non_negative else ""
        raise ValueError(f"{name} must be finite{suffix}")
    return numeric


def residue_key(residue: ResidueId) -> tuple[str, str, str, str, str, str]:
    return (
        residue.structure_id,
        residue.model_id,
        residue.chain_id,
        residue.auth_seq_id,
        residue.insertion_code or "",
        residue.residue_name,
    )


def _freeze_surface_json(
    value: object,
    *,
    path: str = "surface.provenance",
    depth: int = 0,
    node_count: list[int] | None = None,
) -> FrozenJSON:
    if depth > MAX_SURFACE_PROVENANCE_DEPTH:
        raise ValueError(f"{path} exceeds the maximum JSON depth")
    counter = node_count if node_count is not None else [0]
    counter[0] += 1
    if counter[0] > MAX_SURFACE_PROVENANCE_ITEMS:
        raise ValueError(f"{path} exceeds the maximum JSON item limit")
    if isinstance(value, Mapping):
        if len(value) > MAX_SURFACE_PROVENANCE_ITEMS:
            raise ValueError(f"{path} exceeds the maximum JSON item limit")
        frozen: dict[str, FrozenJSON] = {}
        for key, item in value.items():
            if not isinstance(key, str) or not key:
                raise TypeError(f"{path} mapping keys must be non-empty strings")
            if len(key) > MAX_SURFACE_PROVENANCE_STRING_LENGTH:
                raise ValueError(f"{path} exceeds the maximum JSON string length")
            frozen[key] = _freeze_surface_json(
                item,
                path=f"{path}.{key}",
                depth=depth + 1,
                node_count=counter,
            )
        return cast(FrozenJSON, MappingProxyType(frozen))
    if isinstance(value, (tuple, list)):
        if len(value) > MAX_SURFACE_PROVENANCE_ITEMS:
            raise ValueError(f"{path} exceeds the maximum JSON item limit")
        return tuple(
            _freeze_surface_json(
                item,
                path=f"{path}[{index}]",
                depth=depth + 1,
                node_count=counter,
            )
            for index, item in enumerate(value)
        )
    if isinstance(value, str) and len(value) > MAX_SURFACE_PROVENANCE_STRING_LENGTH:
        raise ValueError(f"{path} exceeds the maximum JSON string length")
    return freeze_json(value, path=path)


def _surface_provenance(value: Mapping[str, object] | MethodProvenance) -> FrozenJSON | MethodProvenance:
    if hasattr(value, "to_json"):
        return _freeze_surface_json(cast(Any, value).to_json())
    frozen = _freeze_surface_json(value)
    if not isinstance(frozen, Mapping):
        raise TypeError("surface provenance must be a mapping")
    return frozen


def _surface_units(value: Mapping[str, str]) -> Mapping[str, str]:
    units = freeze_bounded_string_map(value, field_name="surface units")
    if units.get("surface_area") != "angstrom^2":
        raise ValueError("surface measurements must declare surface_area in angstrom^2")
    return units


@dataclass(frozen=True, slots=True)
class PocketSurfaceMeasurement:
    """Source/selection-bound surface measurement for application inputs."""

    candidate_id: str
    source_content_id: str
    selection_id: str
    surface_area_angstrom2: float
    method: str
    provenance: Mapping[str, object] | MethodProvenance
    units: Mapping[str, str] = field(default_factory=lambda: MappingProxyType({"surface_area": "angstrom^2"}))

    def __post_init__(self) -> None:
        for name in ("candidate_id", "source_content_id", "selection_id", "method"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"{name} must be a non-empty string")
        object.__setattr__(
            self,
            "surface_area_angstrom2",
            _finite(self.surface_area_angstrom2, "surface_area_angstrom2", non_negative=True),
        )
        object.__setattr__(self, "provenance", _surface_provenance(self.provenance))
        object.__setattr__(self, "units", _surface_units(self.units))

    def to_json(self) -> dict[str, object]:
        provenance = (
            self.provenance.to_json()
            if isinstance(self.provenance, MethodProvenance)
            else cast(dict[str, object], json_ready(cast(FrozenJSON, self.provenance)))
        )
        return {
            "candidate_id": self.candidate_id,
            "source_content_id": self.source_content_id,
            "selection_id": self.selection_id,
            "surface_area_angstrom2": self.surface_area_angstrom2,
            "method": self.method,
            "provenance": provenance,
            "units": dict(self.units),
        }


def _interaction_record_sort_key(record: InteractionRecord | None) -> tuple[object, ...]:
    if record is None:
        return (0,)
    return (
        1,
        record.structure_id,
        record.interaction_type.value,
        residue_key(record.residue_a),
        residue_key(record.residue_b) if record.residue_b is not None else (),
        record.atom_a or "",
        record.atom_b or "",
        record.distance_angstrom,
        record.angle_degrees if record.angle_degrees is not None else float("-inf"),
        record.ligand_or_metal_id or "",
        record.evidence_mode,
    )


def interaction_sort_key(item: InteractionDifference) -> tuple[object, ...]:
    return (
        item.key.interaction_type.value,
        item.key.reference_position_a,
        item.key.reference_position_b or "",
        item.key.external_partner_id or "",
        item.change.value,
        _interaction_record_sort_key(item.reference_record),
        _interaction_record_sort_key(item.target_record),
    )


def residue_to_json(residue: ResidueId | None) -> dict[str, object] | None:
    if residue is None:
        return None
    return {
        "structure_id": residue.structure_id,
        "model_id": residue.model_id,
        "chain_id": residue.chain_id,
        "auth_seq_id": residue.auth_seq_id,
        "insertion_code": residue.insertion_code,
        "residue_name": residue.residue_name,
    }


def _record_to_json(record: InteractionRecord | None) -> dict[str, object] | None:
    if record is None:
        return None
    return {
        "structure_id": record.structure_id,
        "interaction_type": record.interaction_type.value,
        "residue_a": residue_to_json(record.residue_a),
        "residue_b": residue_to_json(record.residue_b),
        "atom_a": record.atom_a,
        "atom_b": record.atom_b,
        "distance_angstrom": record.distance_angstrom,
        "angle_degrees": record.angle_degrees,
        "ligand_or_metal_id": record.ligand_or_metal_id,
        "evidence_mode": record.evidence_mode,
    }


def interaction_to_json(item: InteractionDifference) -> dict[str, object]:
    return {
        "key": {
            "interaction_type": item.key.interaction_type.value,
            "reference_position_a": item.key.reference_position_a,
            "reference_position_b": item.key.reference_position_b,
            "external_partner_id": item.key.external_partner_id,
        },
        "change": item.change.value,
        "reference_record": _record_to_json(item.reference_record),
        "target_record": _record_to_json(item.target_record),
    }


def _normalise_displacements(values: Sequence[tuple[ResidueId, float]]) -> tuple[tuple[ResidueId, float], ...]:
    output = tuple(
        (residue, _finite(magnitude, "magnitude_angstrom", non_negative=True)) for residue, magnitude in values
    )
    if len({item[0] for item in output}) != len(output):
        raise ValueError("local displacements must not contain duplicate residues")
    return tuple(sorted(output, key=lambda item: residue_key(item[0])))


def _bounded(values: object, name: str) -> tuple[Any, ...]:
    if not isinstance(values, Sequence) or isinstance(values, (str, bytes)):
        raise TypeError(f"{name} must be a bounded Sequence")
    if len(values) > MAX_COMPARISON_ITEMS:
        raise ValueError(f"{name} exceed the bounded comparison limit")
    return tuple(values)


@dataclass(frozen=True, slots=True)
class PocketComparison:
    """Independent evidence channels for one candidate match."""

    availability: Availability
    match: PocketMatch
    volume: PocketVolumeComparison | None = None
    surface_delta_angstrom2: float | None = None
    relative_surface_delta_fraction: float | None = None
    reference_surface_area_angstrom2: float | None = None
    target_surface_area_angstrom2: float | None = None
    reference_surface_method: str | None = None
    target_surface_method: str | None = None
    reference_surface_provenance: FrozenJSON | MethodProvenance | None = None
    target_surface_provenance: FrozenJSON | MethodProvenance | None = None
    reference_surface_units: Mapping[str, str] | None = None
    target_surface_units: Mapping[str, str] | None = None
    lining_residue_conserved: tuple[ResidueId, ...] = ()
    lining_residue_gains: tuple[ResidueId, ...] = ()
    lining_residue_losses: tuple[ResidueId, ...] = ()
    associated_mutations: tuple[MutationEvent, ...] = ()
    interaction_changes: tuple[InteractionDifference, ...] = ()
    local_displacement_angstrom: float | None = None
    ca_displacement_angstrom: float | None = None
    sidechain_displacement_angstrom: float | None = None
    local_displacements: tuple[tuple[ResidueId, float], ...] = ()
    qc_availability: Availability | None = None
    qc_diagnostics: tuple[Diagnostic, ...] = ()
    diagnostics: tuple[Diagnostic, ...] = ()
    units: Mapping[str, str] = field(
        default_factory=lambda: MappingProxyType(
            {
                "volume_delta_angstrom3": "angstrom^3",
                "relative_volume_delta_fraction": "fraction",
                "surface_delta_angstrom2": "angstrom^2",
                "relative_surface_delta_fraction": "fraction",
                "local_displacement_angstrom": "angstrom",
                "ca_displacement_angstrom": "angstrom",
                "sidechain_displacement_angstrom": "angstrom",
            }
        )
    )

    def __post_init__(self) -> None:
        if not isinstance(self.match, PocketMatch):
            raise TypeError("match must be a PocketMatch")
        availability = self.availability
        if not isinstance(availability, Availability):
            availability = Availability(availability)
            object.__setattr__(self, "availability", availability)
        if self.volume is not None and not isinstance(self.volume, PocketVolumeComparison):
            raise TypeError("volume must be a PocketVolumeComparison or None")
        for name in (
            "surface_delta_angstrom2",
            "relative_surface_delta_fraction",
            "reference_surface_area_angstrom2",
            "target_surface_area_angstrom2",
            "local_displacement_angstrom",
            "ca_displacement_angstrom",
            "sidechain_displacement_angstrom",
        ):
            value = getattr(self, name)
            if value is not None:
                object.__setattr__(self, name, _finite(value, name))
        for name in ("reference_surface_method", "target_surface_method"):
            value = getattr(self, name)
            if value is not None and (not isinstance(value, str) or not value.strip()):
                raise ValueError(f"{name} must be a non-empty string when provided")
        for name in ("reference_surface_provenance", "target_surface_provenance"):
            value = getattr(self, name)
            if value is not None and not isinstance(value, MethodProvenance):
                frozen = _freeze_surface_json(value, path=name)
                object.__setattr__(self, name, frozen)
        for name in ("reference_surface_units", "target_surface_units"):
            value = getattr(self, name)
            if value is not None:
                object.__setattr__(self, name, _surface_units(value))
        for name in ("lining_residue_conserved", "lining_residue_gains", "lining_residue_losses"):
            values = _bounded(getattr(self, name), name)
            if any(not isinstance(item, ResidueId) for item in values):
                raise TypeError(f"{name} must contain ResidueId values")
            object.__setattr__(self, name, tuple(sorted(set(values), key=residue_key)))
        mutations = _bounded(self.associated_mutations, "associated_mutations")
        if any(not isinstance(item, MutationEvent) for item in mutations):
            raise TypeError("associated_mutations must contain MutationEvent values")
        object.__setattr__(
            self,
            "associated_mutations",
            tuple(
                sorted(
                    mutations,
                    key=lambda item: (
                        item.alignment_index,
                        item.canonical_notation,
                        item.reference_label,
                        item.target_label,
                    ),
                )
            ),
        )
        changes = _bounded(self.interaction_changes, "interaction_changes")
        if any(not isinstance(item, InteractionDifference) for item in changes):
            raise TypeError("interaction_changes must contain InteractionDifference values")
        object.__setattr__(self, "interaction_changes", tuple(sorted(changes, key=interaction_sort_key)))
        displacements = _bounded(self.local_displacements, "local_displacements")
        if any(
            not isinstance(item, tuple) or len(item) != 2 or not isinstance(item[0], ResidueId)
            for item in displacements
        ):
            raise TypeError("local_displacements must contain (ResidueId, magnitude) pairs")
        object.__setattr__(
            self,
            "local_displacements",
            _normalise_displacements(cast(Sequence[tuple[ResidueId, float]], displacements)),
        )
        if self.qc_availability is not None and not isinstance(self.qc_availability, Availability):
            object.__setattr__(self, "qc_availability", Availability(self.qc_availability))
        for name in ("qc_diagnostics", "diagnostics"):
            values = _bounded(getattr(self, name), name)
            if any(not isinstance(item, Diagnostic) for item in values):
                raise TypeError(f"{name} must contain Diagnostic values")
            object.__setattr__(self, name, values)
        object.__setattr__(self, "units", freeze_bounded_string_map(self.units, field_name="units"))

    @property
    def status(self) -> Availability:
        return self.availability

    @property
    def volume_delta_angstrom3(self) -> float | None:
        return self.volume.delta_angstrom3 if self.volume is not None else None

    @property
    def reference_volume_angstrom3(self) -> float | None:
        return self.volume.reference_volume_angstrom3 if self.volume is not None else None

    @property
    def target_volume_angstrom3(self) -> float | None:
        return self.volume.target_volume_angstrom3 if self.volume is not None else None

    @property
    def relative_volume_delta_fraction(self) -> float | None:
        return self.volume.relative_delta_fraction if self.volume is not None else None

    @property
    def lining_gains(self) -> tuple[ResidueId, ...]:
        return self.lining_residue_gains

    @property
    def lining_losses(self) -> tuple[ResidueId, ...]:
        return self.lining_residue_losses

    @property
    def mutation_associations(self) -> tuple[MutationEvent, ...]:
        return self.associated_mutations

    @property
    def interaction_differences(self) -> tuple[InteractionDifference, ...]:
        return self.interaction_changes

    @property
    def interaction_change(self) -> tuple[InteractionDifference, ...]:
        return self.interaction_changes

    @property
    def ca_displacement(self) -> float | None:
        return self.ca_displacement_angstrom

    @property
    def sidechain_displacement(self) -> float | None:
        return self.sidechain_displacement_angstrom

    def to_json(self) -> dict[str, Any]:
        def surface_provenance(value: FrozenJSON | MethodProvenance | None) -> object:
            if value is None:
                return None
            return value.to_json() if isinstance(value, MethodProvenance) else json_ready(value)

        return {
            "availability": self.availability.value,
            "match": self.match.to_json(),
            "volume": self.volume.to_json() if self.volume is not None else None,
            "surface_delta_angstrom2": self.surface_delta_angstrom2,
            "relative_surface_delta_fraction": self.relative_surface_delta_fraction,
            "reference_surface_area_angstrom2": self.reference_surface_area_angstrom2,
            "target_surface_area_angstrom2": self.target_surface_area_angstrom2,
            "reference_surface_method": self.reference_surface_method,
            "target_surface_method": self.target_surface_method,
            "reference_surface_provenance": surface_provenance(self.reference_surface_provenance),
            "target_surface_provenance": surface_provenance(self.target_surface_provenance),
            "reference_surface_units": dict(self.reference_surface_units)
            if self.reference_surface_units is not None
            else None,
            "target_surface_units": dict(self.target_surface_units) if self.target_surface_units is not None else None,
            "lining_residue_conserved": [residue_to_json(item) for item in self.lining_residue_conserved],
            "lining_residue_gains": [residue_to_json(item) for item in self.lining_residue_gains],
            "lining_residue_losses": [residue_to_json(item) for item in self.lining_residue_losses],
            "associated_mutations": [item.canonical_notation for item in self.associated_mutations],
            "interaction_changes": [interaction_to_json(item) for item in self.interaction_changes],
            "local_displacement_angstrom": self.local_displacement_angstrom,
            "ca_displacement_angstrom": self.ca_displacement_angstrom,
            "sidechain_displacement_angstrom": self.sidechain_displacement_angstrom,
            "local_displacements": [
                {"residue": residue_to_json(item), "magnitude_angstrom": magnitude}
                for item, magnitude in self.local_displacements
            ],
            "qc_availability": self.qc_availability.value if self.qc_availability is not None else None,
            "qc_diagnostics": [item.to_json() for item in self.qc_diagnostics],
            "diagnostics": [item.to_json() for item in self.diagnostics],
            "units": dict(self.units),
        }


__all__ = ["PocketComparison", "PocketSurfaceMeasurement", "interaction_sort_key", "residue_key"]
