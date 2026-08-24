"""Evidence Card builder: aggregates authoritative result objects only."""

from __future__ import annotations

from collections.abc import Sequence

from structlens.core.models import ResidueId
from structlens.core.msa import SequenceResidueRef

from . import EvidenceCard, EvidenceQuality, InteractionEvidence, SequenceEvidence, SiteEvidence, StructureEvidence


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
    ) -> EvidenceCard:
        return EvidenceCard(
            self.reference_residue,
            target_id,
            sequence,
            structure,
            interactions,
            site,
            quality,
            provenance=self.provenance,
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
    provenance: Sequence[str] = (),
) -> EvidenceCard:
    return EvidenceCardBuilder(reference_residue, provenance=provenance).build(
        target_id=target_id,
        sequence=sequence,
        structure=structure,
        interactions=interactions,
        site=site,
        quality=quality,
    )


__all__ = ["EvidenceCardBuilder", "build_evidence_card"]
