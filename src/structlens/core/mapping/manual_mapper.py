"""Explicit locked correspondence builder."""

from __future__ import annotations

from collections.abc import Iterable

from structlens.core.errors import MappingError
from structlens.core.models import (
    CorrespondenceStatus,
    ProteinChain,
    ResidueCorrespondence,
    ResidueId,
)


class ManualResidueMapper:
    def build_correspondence(
        self,
        reference: ProteinChain,
        target: ProteinChain,
        pairs: Iterable[tuple[ResidueId, ResidueId]],
    ) -> list[ResidueCorrespondence]:
        reference_records = {
            record.residue_id: record for record in reference.residue_records
        }
        target_records = {
            record.residue_id: record for record in target.residue_records
        }
        result: list[ResidueCorrespondence] = []
        for index, (reference_id, target_id) in enumerate(pairs):
            if reference_id not in reference_records or target_id not in target_records:
                raise MappingError(
                    "Manual mapping references a residue not present in the selected chains"
                )
            reference_record = reference_records[reference_id]
            target_record = target_records[target_id]
            status = (
                CorrespondenceStatus.CONSERVED
                if reference_record.one_letter == target_record.one_letter
                else CorrespondenceStatus.SUBSTITUTION
            )
            if not reference_record.is_standard or not target_record.is_standard:
                status = CorrespondenceStatus.NONSTANDARD
            result.append(
                ResidueCorrespondence(
                    alignment_index=index,
                    reference=reference_id,
                    target=target_id,
                    reference_one_letter=reference_record.one_letter,
                    target_one_letter=target_record.one_letter,
                    status=status,
                    mapping_source="manual",
                    mapping_locked=True,
                )
            )
        return result


__all__ = ["ManualResidueMapper"]
