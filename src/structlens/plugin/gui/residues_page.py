from collections.abc import Sequence

from structlens.core.evidence import EvidenceCard
from structlens.core.models import MutationEvent, ResidueCorrespondence
from structlens.core.reports import CorrespondenceSnapshot

from .qt_pages import number, residue_label
from .widgets import PageDescriptor

PAGE = PageDescriptor("Residues", "Inspect residue-by-residue correspondence and structural differences.")

def residue_details(
    item: ResidueCorrespondence | CorrespondenceSnapshot,
    mutations: Sequence[MutationEvent],
    card: EvidenceCard | None = None,
) -> str:
    """Describe one correspondence using measurements already in its report."""
    mutation = next((event for event in mutations if event.alignment_index == item.alignment_index), None)
    reference = residue_label(item.reference) if item.reference else "gap"
    target = residue_label(item.target) if item.target else "gap"
    lines = [f"Reference {reference} → Target {target}", f"Alignment status: {item.status.value}"]
    if mutation is not None:
        lines.append(f"{mutation.kind.value.capitalize()}: {mutation.canonical_notation}")
    for label, value in (("Cα displacement", item.ca_displacement_angstrom), ("Backbone RMSD", item.backbone_rmsd_angstrom)):
        lines.append(f"{label}: {number(value)} Å" if value is not None else f"{label}: unavailable")
    if card is None:
        lines.append("Evidence card unavailable for this position.")
    else:
        contacts = card.interactions
        if "interactions" in card.quality.available_sections:
            target_contacts = f"{len(contacts.target_interactions)} target" if item.target else "target gap"
            lines.append(
                f"Contacts recorded: {len(contacts.reference_interactions)} reference, "
                f"{target_contacts}"
            )
        else:
            lines.append("Contacts: unavailable")
        if "sites" in card.quality.available_sections:
            sites = ", ".join(metric.site_id for metric in card.site.metrics)
            lines.append(f"Sites recorded: {sites or 'no association reported in this card'}")
        else:
            lines.append("Sites: unavailable")
        lines.append(f"Evidence quality: {card.quality.overall_status}")
        if card.quality.warnings:
            lines.append(f"Report warnings: {', '.join(dict.fromkeys(card.quality.warnings))}")
        if card.provenance:
            lines.append(f"Source: {' · '.join(card.provenance)}")
    return "\n".join(lines)


__all__ = ["PAGE", "residue_details"]
