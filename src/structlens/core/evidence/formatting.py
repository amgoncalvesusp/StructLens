"""Safe human-readable evidence card formatting."""

from __future__ import annotations

from . import EvidenceCard


def format_evidence_card(card: EvidenceCard) -> str:
    """Format descriptive evidence without introducing an impact score."""
    reference = card.reference_residue
    label = getattr(reference, "residue_name", "unknown")
    target = card.target_id or "unavailable"
    status = card.quality.overall_status
    return f"Reference {label} · target {target} · evidence {status}"


__all__ = ["format_evidence_card"]
