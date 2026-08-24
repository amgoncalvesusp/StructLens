"""Reference-normalized interaction change classification."""

from __future__ import annotations

from collections.abc import Callable, Iterable

from structlens.core.models import ResidueId

from . import InteractionChange, InteractionDifference, InteractionRecord, ReferenceInteractionKey


def compare_interactions(
    reference: Iterable[InteractionRecord],
    target: Iterable[InteractionRecord],
    position_for: Callable[[ResidueId], str | None],
) -> tuple[InteractionDifference, ...]:
    """Compare interactions using mapped reference positions, never raw resi."""

    def key(record: InteractionRecord) -> ReferenceInteractionKey | None:
        first = position_for(record.residue_a)
        second = position_for(record.residue_b) if record.residue_b is not None else None
        if first is None or (record.residue_b is not None and second is None):
            return None
        ordered = tuple(sorted((first, second))) if second is not None else (first, None)
        return ReferenceInteractionKey(record.interaction_type, ordered[0], ordered[1], record.ligand_or_metal_id)

    ref_by_key = {item: record for record in reference if (item := key(record)) is not None}
    target_by_key: dict[ReferenceInteractionKey, InteractionRecord] = {}
    output: list[InteractionDifference] = []
    for record in target:
        item = key(record)
        if item is None:
            output.append(InteractionDifference(ReferenceInteractionKey(record.interaction_type, "unmapped", None, record.ligand_or_metal_id), InteractionChange.TARGET_ONLY_UNMAPPED, None, record))
        else:
            target_by_key[item] = record
    for item in sorted(set(ref_by_key) | set(target_by_key), key=repr):
        ref = ref_by_key.get(item)
        tar = target_by_key.get(item)
        change = InteractionChange.CONSERVED if ref is not None and tar is not None else InteractionChange.LOST if ref is not None else InteractionChange.GAINED
        output.append(InteractionDifference(item, change, ref, tar))
    return tuple(output)


__all__ = ["compare_interactions"]
