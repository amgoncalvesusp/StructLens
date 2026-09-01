"""Immutable pocket evidence snapshots for canonical analysis reports.

The application pocket services deliberately stay outside ``core``.  These
small report values retain their typed outputs at the serialization boundary,
without allowing arbitrary dictionaries to become scientific evidence.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import cast

from structlens.core.evidence import Availability, Diagnostic, PocketEvidenceSections
from structlens.core.models import ResidueId
from structlens.core.pockets.comparison_models import PocketComparison
from structlens.core.pockets.matching import PocketMatchingResult
from structlens.core.pockets.models import (
    AlphaSphereDetectionResult,
    PocketCandidate,
    PocketDetectionSettings,
)
from structlens.core.pockets.volume_models import PocketVolumeResult
from structlens.core.provenance import MethodProvenance, freeze_bounded_string_map

_ROLES = frozenset({"reference", "target"})
_MAX_ITEMS = 100_000


def _role(value: str) -> str:
    if not isinstance(value, str) or value not in _ROLES:
        raise ValueError("pocket evidence role must be 'reference' or 'target'")
    return value


def _sequence(values: object, name: str) -> tuple[object, ...]:
    if not isinstance(values, Sequence) or isinstance(values, (str, bytes)):
        raise TypeError(f"{name} must be a bounded sequence")
    if len(values) > _MAX_ITEMS:
        raise ValueError(f"{name} exceeds the supported item limit")
    return tuple(values)


@dataclass(frozen=True, slots=True)
class PocketDetectionSnapshot:
    """One typed detector output and its candidate lineage."""

    role: str
    availability: Availability
    candidates: tuple[PocketCandidate, ...] = ()
    detection: AlphaSphereDetectionResult | None = None
    settings: PocketDetectionSettings | None = None
    diagnostics: tuple[Diagnostic, ...] = ()
    counts: Mapping[str, int] = field(default_factory=dict)
    provenance: MethodProvenance | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "role", _role(self.role))
        availability = self.availability
        if not isinstance(availability, Availability):
            availability = Availability(availability)
            object.__setattr__(self, "availability", availability)
        candidates = _sequence(self.candidates, "candidates")
        if any(not isinstance(item, PocketCandidate) for item in candidates):
            raise TypeError("candidates must contain PocketCandidate values")
        object.__setattr__(self, "candidates", candidates)
        if self.detection is not None and not isinstance(self.detection, AlphaSphereDetectionResult):
            raise TypeError("detection must be AlphaSphereDetectionResult or None")
        if self.detection is not None and self.detection.availability is not availability:
            raise ValueError("detection availability does not match pocket availability")
        if self.settings is not None and not isinstance(self.settings, PocketDetectionSettings):
            raise TypeError("settings must be PocketDetectionSettings or None")
        diagnostics = _sequence(self.diagnostics, "diagnostics")
        if any(not isinstance(item, Diagnostic) for item in diagnostics):
            raise TypeError("diagnostics must contain Diagnostic values")
        object.__setattr__(self, "diagnostics", diagnostics)
        if not isinstance(self.counts, Mapping):
            raise TypeError("counts must be a mapping")
        normalized_counts: dict[str, int] = {}
        for key, value in self.counts.items():
            if not isinstance(key, str) or not key.strip():
                raise ValueError("count keys must be non-empty strings")
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ValueError("counts must contain non-negative integers")
            normalized_counts[key] = value
        object.__setattr__(self, "counts", MappingProxyType(normalized_counts))
        if self.provenance is not None and not isinstance(self.provenance, MethodProvenance):
            raise TypeError("provenance must be MethodProvenance or None")

    def to_json(self) -> dict[str, object]:
        return {
            "role": self.role,
            "availability": self.availability.value,
            "detection": self.detection.to_json() if self.detection is not None else None,
            "candidates": [item.to_json() for item in self.candidates],
            "settings": self.settings.to_json() if self.settings is not None else None,
            "diagnostics": [item.to_json() for item in self.diagnostics],
            "counts": dict(self.counts),
            "provenance": self.provenance.to_json() if self.provenance is not None else None,
        }


@dataclass(frozen=True, slots=True)
class PocketVolumeSnapshot:
    """One typed free-volume result bound to an input role."""

    role: str
    result: PocketVolumeResult | None = None
    availability: Availability | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "role", _role(self.role))
        if self.result is not None and not isinstance(self.result, PocketVolumeResult):
            raise TypeError("result must be PocketVolumeResult or None")
        availability = self.availability
        if availability is None:
            availability = self.result.availability if self.result is not None else Availability.NOT_APPLICABLE
            object.__setattr__(self, "availability", availability)
        if not isinstance(availability, Availability):
            availability = Availability(availability)
            object.__setattr__(self, "availability", availability)
        if self.result is not None:
            if self.result.availability is not availability:
                raise ValueError("volume result availability does not match pocket availability")

    def to_json(self) -> dict[str, object]:
        return {
            "role": self.role,
            "availability": cast(Availability, self.availability).value,
            "result": self.result.to_json() if self.result is not None else None,
        }


@dataclass(frozen=True, slots=True)
class PocketLiningSnapshot:
    """Lining residue evidence for one role/candidate pair."""

    role: str
    candidate_id: str
    residues: tuple[ResidueId, ...] = ()
    availability: Availability = Availability.AVAILABLE
    diagnostics: tuple[Diagnostic, ...] = ()
    units: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "role", _role(self.role))
        if not isinstance(self.candidate_id, str) or not self.candidate_id.strip():
            raise ValueError("candidate_id must be a non-empty string")
        availability = self.availability
        if not isinstance(availability, Availability):
            availability = Availability(availability)
            object.__setattr__(self, "availability", availability)
        residues = _sequence(self.residues, "residues")
        if any(not isinstance(item, ResidueId) for item in residues):
            raise TypeError("residues must contain ResidueId values")
        object.__setattr__(self, "residues", tuple(sorted(set(residues), key=repr)))
        diagnostics = _sequence(self.diagnostics, "diagnostics")
        if any(not isinstance(item, Diagnostic) for item in diagnostics):
            raise TypeError("diagnostics must contain Diagnostic values")
        object.__setattr__(self, "diagnostics", diagnostics)
        object.__setattr__(self, "units", freeze_bounded_string_map(self.units, field_name="lining units"))

    def to_json(self) -> dict[str, object]:
        return {
            "role": self.role,
            "candidate_id": self.candidate_id,
            "residues": [
                {
                    "structure_id": residue.structure_id,
                    "model_id": residue.model_id,
                    "chain_id": residue.chain_id,
                    "auth_seq_id": residue.auth_seq_id,
                    "insertion_code": residue.insertion_code,
                    "residue_name": residue.residue_name,
                }
                for residue in self.residues
            ],
            "availability": self.availability.value,
            "diagnostics": [item.to_json() for item in self.diagnostics],
            "units": dict(self.units),
        }


@dataclass(frozen=True, slots=True)
class PocketReportSnapshot:
    """Complete typed pocket evidence retained by one ``AnalysisReport``."""

    detections: tuple[PocketDetectionSnapshot, ...] = ()
    volumes: tuple[PocketVolumeSnapshot, ...] = ()
    matching: PocketMatchingResult | None = None
    comparisons: tuple[PocketComparison, ...] = ()
    lining_residues: tuple[PocketLiningSnapshot, ...] = ()
    concordance: PocketEvidenceSections | None = None
    diagnostics: tuple[Diagnostic, ...] = ()
    provenance: MethodProvenance | None = None
    units: Mapping[str, str] = field(default_factory=dict)
    availability: Availability | None = None

    def __post_init__(self) -> None:
        detections = cast(tuple[PocketDetectionSnapshot, ...], _sequence(self.detections, "detections"))
        if any(not isinstance(item, PocketDetectionSnapshot) for item in detections):
            raise TypeError("detections must contain PocketDetectionSnapshot values")
        volumes = cast(tuple[PocketVolumeSnapshot, ...], _sequence(self.volumes, "volumes"))
        if any(not isinstance(item, PocketVolumeSnapshot) for item in volumes):
            raise TypeError("volumes must contain PocketVolumeSnapshot values")
        comparisons = cast(tuple[PocketComparison, ...], _sequence(self.comparisons, "comparisons"))
        if any(not isinstance(item, PocketComparison) for item in comparisons):
            raise TypeError("comparisons must contain PocketComparison values")
        lining = cast(tuple[PocketLiningSnapshot, ...], _sequence(self.lining_residues, "lining_residues"))
        if any(not isinstance(item, PocketLiningSnapshot) for item in lining):
            raise TypeError("lining_residues must contain PocketLiningSnapshot values")
        diagnostics = cast(tuple[Diagnostic, ...], _sequence(self.diagnostics, "diagnostics"))
        if any(not isinstance(item, Diagnostic) for item in diagnostics):
            raise TypeError("diagnostics must contain Diagnostic values")
        if self.matching is not None and not isinstance(self.matching, PocketMatchingResult):
            raise TypeError("matching must be PocketMatchingResult or None")
        if self.concordance is not None and not isinstance(self.concordance, PocketEvidenceSections):
            raise TypeError("concordance must be PocketEvidenceSections or None")
        if self.provenance is not None and not isinstance(self.provenance, MethodProvenance):
            raise TypeError("provenance must be MethodProvenance or None")
        availability = self.availability
        if availability is None:
            states = tuple(item.availability for item in detections + volumes)
            if self.matching is not None:
                states += (self.matching.availability,)
            states += tuple(item.availability for item in comparisons)
            states += tuple(item.availability for item in lining)
            if self.concordance is not None:
                states += tuple(item.availability for item in self.concordance.to_mapping().values())
            availability = (
                Availability.AVAILABLE
                if any(state is Availability.AVAILABLE for state in states)
                else next(
                    (state for state in states if state is not Availability.NOT_APPLICABLE),
                    Availability.NOT_APPLICABLE,
                )
            )
            object.__setattr__(self, "availability", availability)
        elif not isinstance(availability, Availability):
            availability = Availability(availability)
            object.__setattr__(self, "availability", availability)
        has_evidence = bool(detections or volumes or comparisons or lining or diagnostics) or any(
            value is not None for value in (self.matching, self.concordance)
        )
        if availability is Availability.AVAILABLE and not has_evidence:
            raise ValueError("available pocket evidence requires typed evidence")
        child_states = tuple(item.availability for item in detections + volumes + comparisons + lining)
        if self.matching is not None:
            child_states += (self.matching.availability,)
        if self.concordance is not None:
            child_states += tuple(item.availability for item in self.concordance.to_mapping().values())
        if self.availability is not None and any(state is Availability.AVAILABLE for state in child_states):
            if availability is not Availability.AVAILABLE:
                raise ValueError("pocket aggregate availability is incoherent with available child evidence")
        elif self.availability is not None and availability is Availability.AVAILABLE:
            raise ValueError("available pocket aggregate requires an available child evidence channel")
        detection_roles = tuple(item.role for item in detections)
        volume_roles = tuple(item.role for item in volumes)
        if len(set(detection_roles)) != len(detection_roles) or len(set(volume_roles)) != len(volume_roles):
            raise ValueError("pocket detection and volume roles must be unique")
        candidates_by_role: dict[str, set[str]] = {role: set() for role in _ROLES}
        lineage_by_role: dict[str, dict[str, tuple[str | None, str | None]]] = {role: {} for role in _ROLES}
        for detection in detections:
            for candidate in detection.candidates:
                candidates_by_role[detection.role].add(candidate.candidate_id)
                lineage = (candidate.source_content_id, candidate.selection_id)
                previous = lineage_by_role[detection.role].get(candidate.candidate_id)
                if previous is not None and previous != lineage:
                    raise ValueError("detection candidates with one ID must have one lineage tuple")
                lineage_by_role[detection.role][candidate.candidate_id] = lineage

        def require_candidate(candidate: PocketCandidate, role: str, label: str) -> None:
            expected = lineage_by_role[role].get(candidate.candidate_id)
            actual = (candidate.source_content_id, candidate.selection_id)
            if expected is None or actual != expected:
                raise ValueError(f"pocket {label} candidate has no matching detection lineage tuple")

        for volume in volumes:
            if (
                volume.result is not None
                and volume.result.candidate_id is not None
                and volume.result.candidate_id not in candidates_by_role[volume.role]
            ):
                raise ValueError("pocket volume candidate_id has no role-matched detection candidate lineage")
        for lining_snapshot in lining:
            if lining_snapshot.candidate_id not in candidates_by_role[lining_snapshot.role]:
                raise ValueError("pocket lining candidate_id has no role-matched detection candidate lineage")
        if self.matching is not None:
            for match in self.matching.matches:
                if match.reference_candidate is not None:
                    require_candidate(match.reference_candidate, "reference", "match reference")
                if match.target_candidate is not None:
                    require_candidate(match.target_candidate, "target", "match target")
        for comparison in comparisons:
            if comparison.match.reference_candidate is not None:
                require_candidate(comparison.match.reference_candidate, "reference", "comparison reference")
            if comparison.match.target_candidate is not None:
                require_candidate(comparison.match.target_candidate, "target", "comparison target")
        object.__setattr__(self, "detections", detections)
        object.__setattr__(self, "volumes", volumes)
        object.__setattr__(self, "comparisons", comparisons)
        object.__setattr__(self, "lining_residues", lining)
        object.__setattr__(self, "diagnostics", diagnostics)
        object.__setattr__(self, "units", freeze_bounded_string_map(self.units, field_name="pocket units"))

    def to_json(self) -> dict[str, object]:
        if self.availability is None:
            raise RuntimeError("pocket availability was not initialized")
        return {
            "detections": [item.to_json() for item in self.detections],
            "volumes": [item.to_json() for item in self.volumes],
            "matching": self.matching.to_json() if self.matching is not None else None,
            "comparisons": [item.to_json() for item in self.comparisons],
            "lining_residues": [item.to_json() for item in self.lining_residues],
            "concordance": self.concordance.to_json() if self.concordance is not None else None,
            "diagnostics": [item.to_json() for item in self.diagnostics],
            "provenance": self.provenance.to_json() if self.provenance is not None else None,
            "units": dict(self.units),
            "availability": self.availability.value,
        }


__all__ = [
    "PocketDetectionSnapshot",
    "PocketLiningSnapshot",
    "PocketReportSnapshot",
    "PocketVolumeSnapshot",
]
