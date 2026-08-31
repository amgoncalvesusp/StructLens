"""Evidence Card builder and typed pocket concordance channels."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any, cast

from structlens.core.models import ResidueId
from structlens.core.msa import SequenceResidueRef
from structlens.core.provenance import freeze_bounded_string_map, freeze_json

from . import EvidenceCard, EvidenceQuality, InteractionEvidence, SequenceEvidence, SiteEvidence, StructureEvidence
from .status import Availability, Diagnostic

POCKET_CONCORDANCE_CHANNELS = (
    "geometry",
    "ligand_support",
    "parameter_persistence",
    "volume_sensitivity",
    "match_ambiguity",
    "qc",
    "interactions",
)
MAX_EVIDENCE_ITEMS = 100_000
MAX_EVIDENCE_JSON_DEPTH = 64
MAX_EVIDENCE_STRING_LENGTH = 1_000_000


def _freeze_channel_value(
    value: object,
    *,
    path: str = "channel",
    depth: int = 0,
    node_count: list[int] | None = None,
) -> object:
    if depth > MAX_EVIDENCE_JSON_DEPTH:
        raise ValueError(f"{path} exceeds the maximum JSON depth")
    counter = node_count if node_count is not None else [0]
    counter[0] += 1
    if counter[0] > MAX_EVIDENCE_ITEMS:
        raise ValueError(f"{path} exceeds the maximum JSON item limit")
    if hasattr(value, "to_json"):
        payload = cast(Any, value).to_json()
        return _freeze_channel_value(payload, path=f"{path}.to_json", depth=depth + 1, node_count=counter)
    if isinstance(value, Mapping):
        if len(value) > MAX_EVIDENCE_ITEMS:
            raise ValueError(f"{path} exceeds the maximum JSON item limit")
        frozen: dict[str, object] = {}
        for key, item in value.items():
            if not isinstance(key, str) or not key:
                raise TypeError(f"{path} mapping keys must be non-empty strings")
            if len(key) > MAX_EVIDENCE_STRING_LENGTH:
                raise ValueError(f"{path} exceeds the maximum JSON string length")
            frozen[key] = _freeze_channel_value(
                item,
                path=f"{path}.{key}",
                depth=depth + 1,
                node_count=counter,
            )
        return MappingProxyType(frozen)
    if isinstance(value, (tuple, list)):
        if len(value) > MAX_EVIDENCE_ITEMS:
            raise ValueError(f"{path} exceeds the maximum JSON item limit")
        return tuple(
            _freeze_channel_value(
                item,
                path=f"{path}[{index}]",
                depth=depth + 1,
                node_count=counter,
            )
            for index, item in enumerate(value)
        )
    if isinstance(value, str) and len(value) > MAX_EVIDENCE_STRING_LENGTH:
        raise ValueError(f"{path} exceeds the maximum JSON string length")
    return freeze_json(value, path=path)


def _serialise_channel_value(value: object) -> object:
    """Convert typed native measures/collections into JSON-compatible values."""

    if isinstance(value, Mapping):
        return {key: _serialise_channel_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_serialise_channel_value(item) for item in value]
    return value


@dataclass(frozen=True, slots=True)
class PocketEvidenceChannel:
    """One concordance channel with its native measure and availability."""

    availability: Availability
    measure: object | None = None
    units: Mapping[str, str] = field(default_factory=dict)
    diagnostics: tuple[Diagnostic, ...] = ()
    provenance: object | None = None

    def __post_init__(self) -> None:
        availability = self.availability
        if not isinstance(availability, Availability):
            availability = Availability(availability)
            object.__setattr__(self, "availability", availability)
        units = freeze_bounded_string_map(self.units, field_name="channel units")
        if not isinstance(self.diagnostics, Sequence) or isinstance(self.diagnostics, (str, bytes)):
            raise TypeError("channel diagnostics must be a bounded Sequence")
        if len(self.diagnostics) > MAX_EVIDENCE_ITEMS:
            raise ValueError("channel diagnostics exceed the supported limit")
        diagnostics = tuple(self.diagnostics)
        if any(not isinstance(item, Diagnostic) for item in diagnostics):
            raise TypeError("channel diagnostics must contain Diagnostic values")
        if self.measure is not None:
            object.__setattr__(self, "measure", _freeze_channel_value(self.measure, path="channel.measure"))
        if self.provenance is not None:
            object.__setattr__(
                self,
                "provenance",
                _freeze_channel_value(self.provenance, path="channel.provenance"),
            )
        object.__setattr__(self, "units", units)
        object.__setattr__(self, "diagnostics", diagnostics)

    def to_json(self) -> dict[str, object]:
        measure = _serialise_channel_value(self.measure)
        provenance = _serialise_channel_value(self.provenance)
        return {
            "availability": self.availability.value,
            "measure": measure,
            "units": dict(self.units),
            "diagnostics": [item.to_json() for item in self.diagnostics],
            "provenance": provenance,
        }


@dataclass(frozen=True, slots=True)
class PocketEvidenceSections:
    """The seven independent pocket concordance channels."""

    geometry: PocketEvidenceChannel | Availability
    ligand_support: PocketEvidenceChannel | Availability
    parameter_persistence: PocketEvidenceChannel | Availability
    volume_sensitivity: PocketEvidenceChannel | Availability
    match_ambiguity: PocketEvidenceChannel | Availability
    qc: PocketEvidenceChannel | Availability
    interactions: PocketEvidenceChannel | Availability

    def __post_init__(self) -> None:
        for name in POCKET_CONCORDANCE_CHANNELS:
            value = getattr(self, name)
            if isinstance(value, Availability):
                object.__setattr__(self, name, PocketEvidenceChannel(value))
            elif not isinstance(value, PocketEvidenceChannel):
                raise TypeError(f"{name} must be a PocketEvidenceChannel or Availability")

    def to_mapping(self) -> Mapping[str, PocketEvidenceChannel]:
        return {name: getattr(self, name) for name in POCKET_CONCORDANCE_CHANNELS}

    def availability_mapping(self) -> Mapping[str, Availability]:
        return {name: getattr(self, name).availability for name in POCKET_CONCORDANCE_CHANNELS}

    def to_json(self) -> dict[str, object]:
        return {name: getattr(self, name).to_json() for name in POCKET_CONCORDANCE_CHANNELS}


def _normalise_sections(
    sections: PocketEvidenceSections | Mapping[str, Availability | str | PocketEvidenceChannel] | None,
) -> PocketEvidenceSections | None:
    if sections is None:
        return None
    if isinstance(sections, PocketEvidenceSections):
        return sections
    if not isinstance(sections, Mapping):
        raise TypeError("pocket_sections must be PocketEvidenceSections or a mapping")
    unknown = set(sections) - set(POCKET_CONCORDANCE_CHANNELS)
    if unknown:
        raise ValueError(f"unknown pocket evidence section: {sorted(unknown)}")
    values = {name: sections.get(name, Availability.NOT_APPLICABLE) for name in POCKET_CONCORDANCE_CHANNELS}
    return PocketEvidenceSections(**cast(dict[str, object], values))  # type: ignore[arg-type]


def _pocket_quality(sections: PocketEvidenceSections | None, quality: EvidenceQuality | None) -> EvidenceQuality | None:
    if sections is None:
        return quality
    states = sections.availability_mapping()
    available = tuple(name for name in POCKET_CONCORDANCE_CHANNELS if states[name] is Availability.AVAILABLE)
    unavailable = tuple(name for name in POCKET_CONCORDANCE_CHANNELS if states[name] is not Availability.AVAILABLE)
    if quality is None:
        status = "available" if available and not unavailable else "partial" if available else "unavailable"
        return EvidenceQuality(status, available, unavailable, source_count=len(POCKET_CONCORDANCE_CHANNELS))
    status = "available" if available and not unavailable else "partial" if available else "unavailable"
    return EvidenceQuality(
        status, available, unavailable, quality.warnings, quality.coverage_fraction, quality.source_count
    )


class EvidenceCardBuilder:
    def __init__(self, reference_residue: ResidueId | SequenceResidueRef, *, provenance: Sequence[str] = ()) -> None:
        self.reference_residue = reference_residue
        self.provenance = tuple(provenance)

    def build(
        self,
        *,
        target_id: str | None = None,
        sequence: SequenceEvidence | None = None,
        structure: StructureEvidence | None = None,
        interactions: InteractionEvidence | None = None,
        site: SiteEvidence | None = None,
        quality: EvidenceQuality | None = None,
        pocket_sections: PocketEvidenceSections
        | Mapping[str, Availability | str | PocketEvidenceChannel]
        | None = None,
    ) -> EvidenceCard:
        resolved_sections = _normalise_sections(pocket_sections)
        return EvidenceCard(
            self.reference_residue,
            target_id,
            sequence,
            structure,
            interactions,
            site,
            _pocket_quality(resolved_sections, quality),
            provenance=self.provenance,
            pocket_sections=resolved_sections,
        )


def build_evidence_card(
    reference_residue: ResidueId | SequenceResidueRef,
    *,
    target_id: str | None = None,
    sequence: SequenceEvidence | None = None,
    structure: StructureEvidence | None = None,
    interactions: InteractionEvidence | None = None,
    site: SiteEvidence | None = None,
    quality: EvidenceQuality | None = None,
    pocket_sections: PocketEvidenceSections | Mapping[str, Availability | str | PocketEvidenceChannel] | None = None,
    provenance: Sequence[str] = (),
) -> EvidenceCard:
    return EvidenceCardBuilder(reference_residue, provenance=provenance).build(
        target_id=target_id,
        sequence=sequence,
        structure=structure,
        interactions=interactions,
        site=site,
        quality=quality,
        pocket_sections=pocket_sections,
    )


__all__ = [
    "EvidenceCardBuilder",
    "PocketEvidenceChannel",
    "PocketEvidenceSections",
    "POCKET_CONCORDANCE_CHANNELS",
    "build_evidence_card",
]
