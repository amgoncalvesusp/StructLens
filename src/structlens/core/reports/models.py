"""Canonical, immutable report value objects.

The report layer is the serialization boundary for the application.  Mutable
working objects (most notably :class:`ResidueCorrespondence` and NumPy
matrices) are copied into small typed snapshots before they can enter a
report.  A report's identity is the SHA-256 digest of its canonical scientific
JSON, with the identity field itself excluded from the digest input.
"""

from __future__ import annotations

import hashlib
import math
from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import cast

from structlens.core.difference_maps import DistanceDifferenceMatrix, ResidueDisplacementVector
from structlens.core.evidence import (
    Availability,
    Diagnostic,
    EvidenceCard,
    InteractionEvidence,
    SiteEvidence,
)
from structlens.core.interactions import InteractionDifference, InteractionRecord
from structlens.core.models import (
    AnalysisResult,
    ComparisonMode,
    CorrespondenceStatus,
    MutationEvent,
    ResidueCorrespondence,
    ResidueId,
    StructuralTransform,
)
from structlens.core.msa import MultipleSequenceAlignment
from structlens.core.parsing import InputSelection
from structlens.core.provenance import MethodProvenance
from structlens.core.quality import StructureQualityReport
from structlens.core.sites import SiteDefinition, SiteMetrics

from .serialization import (
    JSONValue,
    canonical_bytes,
    copy_msa,
    evidence_card_json,
    interaction_json,
    interaction_record_json,
    matrix_rows,
    msa_json,
    mutation_json,
    residue_json,
    selection_json,
    site_definition_json,
    site_json,
    sites_present,
    transform_json,
    vector_json,
)


def _finite(value: float | None, name: str) -> float | None:
    if value is None:
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite or None")
    return result


def _text(value: str, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty string")
    return value.strip()


def _freeze_string_map(values: Mapping[str, str], name: str) -> Mapping[str, str]:
    if not isinstance(values, Mapping):
        raise TypeError(f"{name} must be a mapping")
    normalized: dict[str, str] = {}
    for key, value in values.items():
        normalized[_text(key, f"{name} key")] = _text(value, f"{name} value")
    return MappingProxyType(normalized)


@dataclass(frozen=True, slots=True)
class CorrespondenceSnapshot:
    """Deep immutable copy of one alignment correspondence."""

    alignment_index: int
    reference: ResidueId | None
    target: ResidueId | None
    reference_one_letter: str | None
    target_one_letter: str | None
    status: CorrespondenceStatus
    sequence_score: float | None = None
    ca_displacement_angstrom: float | None = None
    backbone_rmsd_angstrom: float | None = None
    sidechain_rmsd_angstrom: float | None = None
    all_heavy_atom_rmsd_angstrom: float | None = None
    is_outlier: bool = False
    is_key_residue: bool = False
    mapping_source: str = "unknown"
    mapping_locked: bool = False

    def __post_init__(self) -> None:
        if isinstance(self.status, str):
            object.__setattr__(self, "status", CorrespondenceStatus(self.status))
        if self.alignment_index < 0:
            raise ValueError("alignment_index must be non-negative")
        for name in (
            "sequence_score",
            "ca_displacement_angstrom",
            "backbone_rmsd_angstrom",
            "sidechain_rmsd_angstrom",
            "all_heavy_atom_rmsd_angstrom",
        ):
            object.__setattr__(self, name, _finite(getattr(self, name), name))
        if not isinstance(self.mapping_source, str) or not self.mapping_source.strip():
            raise ValueError("mapping_source must be non-empty")

    @classmethod
    def from_record(cls, record: ResidueCorrespondence) -> CorrespondenceSnapshot:
        if not isinstance(record, ResidueCorrespondence):
            raise TypeError("record must be a ResidueCorrespondence")
        return cls(
            alignment_index=record.alignment_index,
            reference=record.reference,
            target=record.target,
            reference_one_letter=record.reference_one_letter,
            target_one_letter=record.target_one_letter,
            status=record.status,
            sequence_score=record.sequence_score,
            ca_displacement_angstrom=record.ca_displacement_angstrom,
            backbone_rmsd_angstrom=record.backbone_rmsd_angstrom,
            sidechain_rmsd_angstrom=record.sidechain_rmsd_angstrom,
            all_heavy_atom_rmsd_angstrom=record.all_heavy_atom_rmsd_angstrom,
            is_outlier=record.is_outlier,
            is_key_residue=record.is_key_residue,
            mapping_source=record.mapping_source,
            mapping_locked=record.mapping_locked,
        )

    def to_json(self) -> dict[str, JSONValue]:
        return {
            "alignment_index": self.alignment_index,
            "reference": residue_json(self.reference),
            "target": residue_json(self.target),
            "reference_one_letter": self.reference_one_letter,
            "target_one_letter": self.target_one_letter,
            "status": self.status.value,
            "sequence_score": self.sequence_score,
            "ca_displacement_angstrom": self.ca_displacement_angstrom,
            "backbone_rmsd_angstrom": self.backbone_rmsd_angstrom,
            "sidechain_rmsd_angstrom": self.sidechain_rmsd_angstrom,
            "all_heavy_atom_rmsd_angstrom": self.all_heavy_atom_rmsd_angstrom,
            "is_outlier": self.is_outlier,
            "is_key_residue": self.is_key_residue,
            "mapping_source": self.mapping_source,
            "mapping_locked": self.mapping_locked,
        }


@dataclass(frozen=True, slots=True)
class AnalysisSnapshot:
    """Deep immutable, JSON-ready snapshot of an :class:`AnalysisResult`."""

    reference_id: str
    target_id: str
    correspondences: tuple[CorrespondenceSnapshot, ...]
    mutations: tuple[MutationEvent, ...]
    sequence_identity: float
    sequence_coverage: float
    alignment_decision: str
    sequence_similarity: float | None = None
    strict_rmsd_angstrom: float | None = None
    refined_rmsd_angstrom: float | None = None
    mapped_residue_count: int = 0
    refined_residue_count: int | None = None
    excluded_alignment_indices: tuple[int, ...] = ()
    tm_score: float | None = None
    legacy_provenance: Mapping[str, str] = field(default_factory=dict)
    transform: StructuralTransform | None = None
    method_provenance: MethodProvenance | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "reference_id", _text(self.reference_id, "reference_id"))
        object.__setattr__(self, "target_id", _text(self.target_id, "target_id"))
        correspondences = tuple(self.correspondences)
        if any(not isinstance(item, CorrespondenceSnapshot) for item in correspondences):
            raise TypeError("correspondences must contain CorrespondenceSnapshot values")
        object.__setattr__(self, "correspondences", correspondences)
        mutations = tuple(self.mutations)
        if any(not isinstance(item, MutationEvent) for item in mutations):
            raise TypeError("mutations must contain MutationEvent values")
        object.__setattr__(self, "mutations", mutations)
        for name in (
            "sequence_identity",
            "sequence_coverage",
            "sequence_similarity",
            "strict_rmsd_angstrom",
            "refined_rmsd_angstrom",
            "tm_score",
        ):
            object.__setattr__(self, name, _finite(getattr(self, name), name))
        for name in ("mapped_residue_count", "refined_residue_count"):
            value = getattr(self, name)
            if value is not None and (isinstance(value, bool) or value < 0):
                raise ValueError(f"{name} must be non-negative or None")
        object.__setattr__(self, "excluded_alignment_indices", tuple(self.excluded_alignment_indices))
        object.__setattr__(self, "legacy_provenance", _freeze_string_map(self.legacy_provenance, "legacy_provenance"))
        if self.method_provenance is not None and not isinstance(self.method_provenance, MethodProvenance):
            raise TypeError("method_provenance must be MethodProvenance or None")

    @classmethod
    def from_result(cls, result: AnalysisResult) -> AnalysisSnapshot:
        if not isinstance(result, AnalysisResult):
            raise TypeError("result must be an AnalysisResult")
        return cls(
            reference_id=result.reference_id,
            target_id=result.target_id,
            correspondences=tuple(CorrespondenceSnapshot.from_record(item) for item in result.correspondences),
            mutations=tuple(result.mutations),
            sequence_identity=result.sequence_identity,
            sequence_coverage=result.sequence_coverage,
            alignment_decision=result.alignment_decision,
            sequence_similarity=result.sequence_similarity,
            strict_rmsd_angstrom=result.strict_rmsd_angstrom,
            refined_rmsd_angstrom=result.refined_rmsd_angstrom,
            mapped_residue_count=result.mapped_residue_count,
            refined_residue_count=result.refined_residue_count,
            excluded_alignment_indices=tuple(result.excluded_alignment_indices),
            tm_score=result.tm_score,
            legacy_provenance=dict(result.provenance),
            transform=result.transform,
            method_provenance=result.method_provenance,
        )

    def to_json(self) -> dict[str, JSONValue]:
        return {
            "reference_id": self.reference_id,
            "target_id": self.target_id,
            "correspondences": [item.to_json() for item in self.correspondences],
            "mutations": [mutation_json(item) for item in self.mutations],
            "sequence_identity": self.sequence_identity,
            "sequence_coverage": self.sequence_coverage,
            "alignment_decision": self.alignment_decision,
            "sequence_similarity": self.sequence_similarity,
            "strict_rmsd_angstrom": self.strict_rmsd_angstrom,
            "refined_rmsd_angstrom": self.refined_rmsd_angstrom,
            "mapped_residue_count": self.mapped_residue_count,
            "refined_residue_count": self.refined_residue_count,
            "excluded_alignment_indices": list(self.excluded_alignment_indices),
            "tm_score": self.tm_score,
            "legacy_provenance": dict(self.legacy_provenance),
            "transform": transform_json(self.transform),
            "method_provenance": (
                cast(JSONValue, self.method_provenance.to_json()) if self.method_provenance is not None else None
            ),
        }


@dataclass(frozen=True, slots=True)
class DistanceMapSnapshot:
    """Deep immutable copy of a distance-difference matrix."""

    reference_positions: tuple[str, ...]
    reference_distances_angstrom: tuple[tuple[float, ...], ...]
    target_distances_angstrom: tuple[tuple[float, ...], ...]
    delta_angstrom: tuple[tuple[float, ...], ...]
    valid_mask: tuple[tuple[bool, ...], ...]

    @classmethod
    def from_matrix(cls, matrix: DistanceDifferenceMatrix) -> DistanceMapSnapshot:
        if not isinstance(matrix, DistanceDifferenceMatrix):
            raise TypeError("matrix must be a DistanceDifferenceMatrix")
        return cls(
            tuple(matrix.reference_positions),
            cast(tuple[tuple[float, ...], ...], matrix_rows(matrix.reference_distances_angstrom, "reference_distances_angstrom")),
            cast(tuple[tuple[float, ...], ...], matrix_rows(matrix.target_distances_angstrom, "target_distances_angstrom")),
            cast(tuple[tuple[float, ...], ...], matrix_rows(matrix.delta_angstrom, "delta_angstrom")),
            cast(tuple[tuple[bool, ...], ...], matrix_rows(matrix.valid_mask, "valid_mask", boolean=True)),
        )

    def __post_init__(self) -> None:
        positions = tuple(_text(item, "reference position") for item in self.reference_positions)
        size = len(positions)
        object.__setattr__(self, "reference_positions", positions)
        for name in ("reference_distances_angstrom", "target_distances_angstrom", "delta_angstrom"):
            rows = tuple(tuple(float(item) for item in row) for row in getattr(self, name))
            if len(rows) != size or any(len(row) != size for row in rows):
                raise ValueError(f"{name} shape must match reference_positions")
            if any(not math.isfinite(item) for row in rows for item in row):
                raise ValueError(f"{name} must contain finite values")
            object.__setattr__(self, name, rows)
        mask = tuple(tuple(bool(item) for item in row) for row in self.valid_mask)
        if len(mask) != size or any(len(row) != size for row in mask):
            raise ValueError("valid_mask shape must match reference_positions")
        object.__setattr__(self, "valid_mask", mask)

    def to_json(self) -> dict[str, JSONValue]:
        return {
            "reference_positions": list(self.reference_positions),
            "reference_distances_angstrom": [list(row) for row in self.reference_distances_angstrom],
            "target_distances_angstrom": [list(row) for row in self.target_distances_angstrom],
            "delta_angstrom": [list(row) for row in self.delta_angstrom],
            "valid_mask": [list(row) for row in self.valid_mask],
        }


@dataclass(frozen=True, slots=True)
class InputQualityBundle:
    """Quality reports for the reference and target input structures."""

    reference: StructureQualityReport
    target: StructureQualityReport

    def __post_init__(self) -> None:
        if not isinstance(self.reference, StructureQualityReport) or not isinstance(self.target, StructureQualityReport):
            raise TypeError("reference and target must be StructureQualityReport values")

    def to_json(self) -> dict[str, JSONValue]:
        return {"reference": cast(JSONValue, self.reference.to_json()), "target": cast(JSONValue, self.target.to_json())}


def _availability(value: Availability | str, name: str) -> Availability:
    try:
        return value if isinstance(value, Availability) else Availability(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"unknown {name}: {value!r}") from exc


@dataclass(frozen=True, slots=True)
class SectionAvailability:
    """Typed availability state for every optional report section."""

    input_quality: Availability = Availability.NOT_APPLICABLE
    analysis: Availability = Availability.NOT_APPLICABLE
    msa: Availability = Availability.NOT_APPLICABLE
    interactions: Availability = Availability.NOT_APPLICABLE
    sites: Availability = Availability.NOT_APPLICABLE
    distance_map: Availability = Availability.NOT_APPLICABLE
    displacement_vectors: Availability = Availability.NOT_APPLICABLE
    evidence_cards: Availability = Availability.NOT_APPLICABLE
    pockets: Availability = Availability.NOT_APPLICABLE

    def __post_init__(self) -> None:
        for name in (
            "input_quality",
            "analysis",
            "msa",
            "interactions",
            "sites",
            "distance_map",
            "displacement_vectors",
            "evidence_cards",
            "pockets",
        ):
            object.__setattr__(self, name, _availability(getattr(self, name), name))

    def to_json(self) -> dict[str, JSONValue]:
        return {name: getattr(self, name).value for name in self.__dataclass_fields__}

    @property
    def vectors(self) -> Availability:
        """Compatibility spelling used by the v0.3 presentation layer."""

        return self.displacement_vectors


@dataclass(frozen=True, slots=True)
class AnalysisReport:
    """The single immutable scientific artifact consumed by all presenters."""

    reference_selection: InputSelection
    target_selection: InputSelection
    input_quality: InputQualityBundle
    analysis: AnalysisSnapshot | None = None
    msa: MultipleSequenceAlignment | None = None
    interactions: tuple[InteractionDifference | InteractionRecord, ...] | InteractionEvidence | None = None
    sites: tuple[SiteMetrics, ...] | SiteEvidence | None = None
    distance_map: DistanceMapSnapshot | None = None
    displacement_vectors: tuple[ResidueDisplacementVector, ...] = ()
    evidence_cards: tuple[EvidenceCard, ...] = ()
    site_definitions: tuple[SiteDefinition, ...] = ()
    schema_version: str = "4.0"
    comparison_mode: ComparisonMode = ComparisonMode.PAIRWISE
    availability: SectionAvailability = field(default_factory=SectionAvailability)
    diagnostics: tuple[Diagnostic, ...] = ()
    provenance: MethodProvenance | None = None
    report_id: str = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.reference_selection, InputSelection) or not isinstance(self.target_selection, InputSelection):
            raise TypeError("reference_selection and target_selection must be InputSelection values")
        if not isinstance(self.input_quality, InputQualityBundle):
            raise TypeError("input_quality must be an InputQualityBundle")
        if self.analysis is not None and not isinstance(self.analysis, AnalysisSnapshot):
            raise TypeError("analysis must be an AnalysisSnapshot or None")
        if self.msa is not None:
            object.__setattr__(self, "msa", copy_msa(self.msa))
        if self.interactions is not None:
            if isinstance(self.interactions, InteractionEvidence):
                pass
            else:
                interactions = tuple(self.interactions)
                if any(not isinstance(item, (InteractionDifference, InteractionRecord)) for item in interactions):
                    raise TypeError("interactions must contain typed interaction values")
                object.__setattr__(self, "interactions", interactions)
        if self.sites is not None:
            if isinstance(self.sites, SiteEvidence):
                pass
            else:
                sites = tuple(self.sites)
                if any(not isinstance(item, SiteMetrics) for item in sites):
                    raise TypeError("sites must contain SiteMetrics values")
                object.__setattr__(self, "sites", sites)
        if self.distance_map is not None and not isinstance(self.distance_map, DistanceMapSnapshot):
            raise TypeError("distance_map must be a DistanceMapSnapshot or None")
        vectors = tuple(self.displacement_vectors)
        if any(not isinstance(item, ResidueDisplacementVector) for item in vectors):
            raise TypeError("vectors must contain ResidueDisplacementVector values")
        object.__setattr__(self, "displacement_vectors", vectors)
        cards = tuple(self.evidence_cards)
        if any(not isinstance(item, EvidenceCard) for item in cards):
            raise TypeError("evidence_cards must contain EvidenceCard values")
        object.__setattr__(self, "evidence_cards", cards)
        definitions = tuple(self.site_definitions)
        if any(not isinstance(item, SiteDefinition) for item in definitions):
            raise TypeError("site_definitions must contain SiteDefinition values")
        if len({item.site_id for item in definitions}) != len(definitions):
            raise ValueError("site_definitions must have unique site_id values")
        object.__setattr__(self, "site_definitions", definitions)
        object.__setattr__(self, "schema_version", _text(self.schema_version, "schema_version"))
        if not isinstance(self.comparison_mode, ComparisonMode):
            object.__setattr__(self, "comparison_mode", ComparisonMode(self.comparison_mode))
        if not isinstance(self.availability, SectionAvailability):
            raise TypeError("availability must be SectionAvailability")
        diagnostics = tuple(self.diagnostics)
        if any(not isinstance(item, Diagnostic) for item in diagnostics):
            raise TypeError("diagnostics must contain Diagnostic values")
        object.__setattr__(self, "diagnostics", diagnostics)
        if self.provenance is not None and not isinstance(self.provenance, MethodProvenance):
            raise TypeError("provenance must be MethodProvenance or None")
        self._validate_availability_coherence()
        self._validate_content_hashes()
        scientific = self._payload(include_report_id=False)
        object.__setattr__(self, "report_id", hashlib.sha256(canonical_bytes(scientific)).hexdigest())

    def _validate_content_hashes(self) -> None:
        if self.provenance is None:
            return
        hashes = {key.casefold().replace("-", "_"): value for key, value in self.provenance.input_hashes.items()}
        for label, expected in (
            ("reference", self.reference_selection.content_id),
            ("target", self.target_selection.content_id),
        ):
            for key, value in hashes.items():
                if key in {f"{label}_content", f"{label}_content_id", f"{label}_selection_content"}:
                    if value.lower() != expected.lower():
                        raise ValueError(f"{label} selection content hash is incoherent with the selected content")

    def _validate_availability_coherence(self) -> None:
        checks = (
            ("analysis", self.analysis is not None),
            ("msa", self.msa is not None),
            ("interactions", self.interactions is not None),
            ("sites", sites_present(self.sites)),
            ("distance_map", self.distance_map is not None),
            ("evidence_cards", bool(self.evidence_cards)),
        )
        for name, present in checks:
            state = getattr(self.availability, name)
            if state is Availability.AVAILABLE and not present:
                raise ValueError(f"{name} availability is available but its typed payload is absent")
            if state is not Availability.AVAILABLE and present:
                raise ValueError(f"{name} availability is {state.value} but its typed payload is present")
        if self.availability.input_quality is Availability.AVAILABLE and (
            self.input_quality.reference.availability is not Availability.AVAILABLE
            or self.input_quality.target.availability is not Availability.AVAILABLE
        ):
            raise ValueError("input_quality availability is available but an input quality report is unavailable")

    def _payload(self, *, include_report_id: bool) -> dict[str, JSONValue]:
        interactions: JSONValue
        if isinstance(self.interactions, InteractionEvidence):
            interactions = {
                "differences": [interaction_json(item) for item in self.interactions.differences],
                "reference_interactions": [interaction_record_json(item) for item in self.interactions.reference_interactions],
                "target_interactions": [interaction_record_json(item) for item in self.interactions.target_interactions],
            }
        elif self.interactions is None:
            interactions = None
        else:
            interactions = [interaction_json(item) for item in self.interactions]
        if isinstance(self.sites, SiteEvidence):
            site_payload: JSONValue = {"metrics": [site_json(item) for item in self.sites.metrics]}
        elif self.sites is None:
            site_payload = None
        else:
            site_payload = [site_json(item) for item in self.sites]
        payload: dict[str, JSONValue] = {
            "schema_version": self.schema_version,
            "comparison_mode": self.comparison_mode.value,
            "reference_selection": selection_json(self.reference_selection),
            "target_selection": selection_json(self.target_selection),
            "input_quality": cast(JSONValue, self.input_quality.to_json()),
            "analysis": cast(JSONValue, self.analysis.to_json()) if self.analysis is not None else None,
            "msa": msa_json(self.msa),
            "interactions": interactions,
            "sites": site_payload,
            "distance_map": cast(JSONValue, self.distance_map.to_json()) if self.distance_map is not None else None,
            "displacement_vectors": [vector_json(item) for item in self.displacement_vectors],
            "evidence_cards": [evidence_card_json(item) for item in self.evidence_cards],
            "site_definitions": [site_definition_json(item) for item in self.site_definitions],
            "availability": self.availability.to_json(),
            "diagnostics": [cast(JSONValue, item.to_json()) for item in self.diagnostics],
            "provenance": cast(JSONValue, self.provenance.to_json()) if self.provenance is not None else None,
        }
        if include_report_id:
            payload["report_id"] = self.report_id
        return payload

    def to_json(self) -> dict[str, JSONValue]:
        return self._payload(include_report_id=True)

    def canonical_json_bytes(self) -> bytes:
        return canonical_bytes(self.to_json())

    @property
    def vectors(self) -> tuple[ResidueDisplacementVector, ...]:
        """Compatibility spelling used by the v0.3 presentation layer."""

        return self.displacement_vectors


__all__ = [
    "AnalysisReport",
    "AnalysisSnapshot",
    "CorrespondenceSnapshot",
    "DistanceMapSnapshot",
    "InputQualityBundle",
    "SectionAvailability",
]
