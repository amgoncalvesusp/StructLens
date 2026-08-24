"""Residue Evidence Card contracts.

Evidence cards are a presentation-neutral aggregation of stored scientific
evidence.  They intentionally contain no impact, damage, functional,
pathogenic, or binding-affinity score.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum

from structlens.core.interactions import InteractionDifference, InteractionRecord
from structlens.core.models import ResidueId
from structlens.core.msa import SequenceResidueRef
from structlens.core.sites import SiteMetrics


class EvidenceAvailability(str, Enum):
    AVAILABLE = "available"
    PARTIAL = "partial"
    UNAVAILABLE = "unavailable"


def _finite(value: float | None, name: str, *, minimum: float | None = None) -> None:
    if value is None:
        return
    if not math.isfinite(value) or (minimum is not None and value < minimum):
        suffix = " and non-negative" if minimum is not None else ""
        raise ValueError(f"{name} must be finite{suffix}")


@dataclass(frozen=True, slots=True)
class SequenceEvidence:
    """Sequence and MSA facts for one residue position."""

    reference_one_letter: str | None = None
    target_one_letter: str | None = None
    alignment_index: int | None = None
    sequence_identity: float | None = None
    conservation_fraction: float | None = None
    entropy_bits: float | None = None
    gap_fraction: float | None = None
    ambiguous_fraction: float | None = None
    source_refs: tuple[SequenceResidueRef, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        if self.alignment_index is not None and self.alignment_index < 0:
            raise ValueError("alignment_index must be non-negative or None")
        for name in ("sequence_identity", "conservation_fraction", "gap_fraction", "ambiguous_fraction"):
            value = getattr(self, name)
            _finite(value, name)
            if value is not None and not 0.0 <= value <= 1.0:
                raise ValueError(f"{name} must be between 0 and 1")
        _finite(self.entropy_bits, "entropy_bits", minimum=0.0)
        refs = tuple(self.source_refs)
        if any(not isinstance(item, SequenceResidueRef) for item in refs):
            raise TypeError("source_refs must contain SequenceResidueRef values")
        object.__setattr__(self, "source_refs", refs)


@dataclass(frozen=True, slots=True)
class StructureEvidence:
    """Structure-alignment and local geometry facts."""

    ca_displacement_angstrom: float | None = None
    backbone_rmsd_angstrom: float | None = None
    sidechain_rmsd_angstrom: float | None = None
    all_heavy_atom_rmsd_angstrom: float | None = None
    sasa_reference_angstrom2: float | None = None
    sasa_target_angstrom2: float | None = None
    available: bool = False

    def __post_init__(self) -> None:
        for name in (
            "ca_displacement_angstrom",
            "backbone_rmsd_angstrom",
            "sidechain_rmsd_angstrom",
            "all_heavy_atom_rmsd_angstrom",
            "sasa_reference_angstrom2",
            "sasa_target_angstrom2",
        ):
            _finite(getattr(self, name), name, minimum=0.0)


@dataclass(frozen=True, slots=True)
class InteractionEvidence:
    """Interaction observations and reference-normalized differences."""

    differences: tuple[InteractionDifference, ...] = field(default_factory=tuple)
    reference_interactions: tuple[InteractionRecord, ...] = field(default_factory=tuple)
    target_interactions: tuple[InteractionRecord, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        differences = tuple(self.differences)
        reference = tuple(self.reference_interactions)
        target = tuple(self.target_interactions)
        if any(not isinstance(item, InteractionDifference) for item in differences):
            raise TypeError("differences must contain InteractionDifference values")
        if any(not isinstance(item, InteractionRecord) for item in reference + target):
            raise TypeError("interaction collections must contain InteractionRecord values")
        object.__setattr__(self, "differences", differences)
        object.__setattr__(self, "reference_interactions", reference)
        object.__setattr__(self, "target_interactions", target)


@dataclass(frozen=True, slots=True)
class SiteEvidence:
    """Site definitions' computed metrics for a residue card."""

    metrics: tuple[SiteMetrics, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        values = tuple(self.metrics)
        if any(not isinstance(item, SiteMetrics) for item in values):
            raise TypeError("metrics must contain SiteMetrics values")
        object.__setattr__(self, "metrics", values)

    @property
    def site_metrics(self) -> tuple[SiteMetrics, ...]:
        return self.metrics


@dataclass(frozen=True, slots=True)
class EvidenceQuality:
    """Categorical quality/availability metadata for an Evidence Card."""

    overall_status: str = EvidenceAvailability.UNAVAILABLE.value
    available_sections: tuple[str, ...] = field(default_factory=tuple)
    unavailable_sections: tuple[str, ...] = field(default_factory=tuple)
    warnings: tuple[str, ...] = field(default_factory=tuple)
    coverage_fraction: float | None = None
    source_count: int | None = None

    def __post_init__(self) -> None:
        if self.overall_status not in {item.value for item in EvidenceAvailability}:
            raise ValueError("overall_status must be available, partial, or unavailable")
        _finite(self.coverage_fraction, "coverage_fraction")
        if self.coverage_fraction is not None and not 0.0 <= self.coverage_fraction <= 1.0:
            raise ValueError("coverage_fraction must be between 0 and 1")
        if self.source_count is not None and self.source_count < 0:
            raise ValueError("source_count must be non-negative or None")
        object.__setattr__(self, "available_sections", tuple(self.available_sections))
        object.__setattr__(self, "unavailable_sections", tuple(self.unavailable_sections))
        object.__setattr__(self, "warnings", tuple(self.warnings))


@dataclass(frozen=True, slots=True, init=False)
class EvidenceCard:
    """All evidence sections for one authoritative sequence residue."""

    reference_residue: ResidueId | SequenceResidueRef
    target_id: str | None
    residue_ref: SequenceResidueRef
    sequence: SequenceEvidence = field(default_factory=SequenceEvidence)
    structure: StructureEvidence = field(default_factory=StructureEvidence)
    interactions: InteractionEvidence = field(default_factory=InteractionEvidence)
    site: SiteEvidence = field(default_factory=SiteEvidence)
    quality: EvidenceQuality = field(default_factory=EvidenceQuality)
    schema_version: str = "3.0"
    provenance: tuple[str, ...] = field(default_factory=tuple)

    def __init__(
        self,
        reference_residue: ResidueId | SequenceResidueRef,
        target_id: str | None = None,
        sequence: SequenceEvidence | None = None,
        structure: StructureEvidence | None = None,
        interactions: InteractionEvidence | None = None,
        site: SiteEvidence | None = None,
        quality: EvidenceQuality | None = None,
        schema_version: str = "3.0",
        provenance: tuple[str, ...] = (),
        **legacy: object,
    ) -> None:
        if "residue_ref" in legacy:
            raw_reference = legacy.pop("residue_ref")
            if not isinstance(raw_reference, (ResidueId, SequenceResidueRef)):
                raise TypeError("residue_ref must be ResidueId or SequenceResidueRef")
            reference_residue = raw_reference
        if legacy:
            raise TypeError(f"unexpected EvidenceCard fields: {', '.join(sorted(legacy))}")
        if isinstance(reference_residue, SequenceResidueRef):
            ref = reference_residue
        else:
            if not isinstance(reference_residue, ResidueId):
                raise TypeError("reference_residue must be ResidueId or SequenceResidueRef")
            ref = SequenceResidueRef(0, reference_residue.residue_name[:1] or "X", reference_residue)
        object.__setattr__(self, "reference_residue", ref.residue_id or ref)
        object.__setattr__(self, "target_id", target_id)
        object.__setattr__(self, "residue_ref", ref)
        object.__setattr__(self, "sequence", sequence or SequenceEvidence())
        object.__setattr__(self, "structure", structure or StructureEvidence())
        object.__setattr__(self, "interactions", interactions or InteractionEvidence())
        object.__setattr__(self, "site", site or SiteEvidence())
        object.__setattr__(self, "quality", quality or EvidenceQuality())
        object.__setattr__(self, "schema_version", schema_version)
        object.__setattr__(self, "provenance", tuple(provenance))
        self.__post_init__()

    def __post_init__(self) -> None:
        if not isinstance(self.residue_ref, SequenceResidueRef):
            raise TypeError("residue_ref must be a SequenceResidueRef")
        for name, expected in (
            ("sequence", SequenceEvidence),
            ("structure", StructureEvidence),
            ("interactions", InteractionEvidence),
            ("site", SiteEvidence),
            ("quality", EvidenceQuality),
        ):
            if not isinstance(getattr(self, name), expected):
                raise TypeError(f"{name} must be a {expected.__name__}")
        if not self.schema_version:
            raise ValueError("schema_version must not be empty")
        object.__setattr__(self, "provenance", tuple(self.provenance))

    @property
    def residue_id(self) -> ResidueId | None:
        return self.residue_ref.residue_id

    @property
    def sequence_evidence(self) -> SequenceEvidence:
        return self.sequence

    @property
    def structure_evidence(self) -> StructureEvidence:
        return self.structure

    @property
    def interaction_evidence(self) -> InteractionEvidence:
        return self.interactions

    @property
    def site_evidence(self) -> SiteEvidence:
        return self.site


ResidueEvidenceCard = EvidenceCard

from .builder import EvidenceCardBuilder, build_evidence_card  # noqa: E402
from .completeness import quality_for_sections  # noqa: E402
from .formatting import format_evidence_card  # noqa: E402

__all__ = [
    "EvidenceAvailability",
    "EvidenceCard",
    "EvidenceQuality",
    "InteractionEvidence",
    "ResidueEvidenceCard",
    "SequenceEvidence",
    "SiteEvidence",
    "StructureEvidence",
    "EvidenceCardBuilder",
    "build_evidence_card",
    "format_evidence_card",
    "quality_for_sections",
]
