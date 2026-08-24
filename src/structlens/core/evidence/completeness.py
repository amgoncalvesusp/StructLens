"""Categorical completeness helpers for evidence cards."""

from __future__ import annotations

from collections.abc import Iterable

from . import EvidenceAvailability, EvidenceQuality


def quality_for_sections(available: Iterable[str], unavailable: Iterable[str], *, warnings: Iterable[str] = ()) -> EvidenceQuality:
    present = tuple(available)
    missing = tuple(unavailable)
    if present and not missing:
        status = EvidenceAvailability.AVAILABLE.value
    elif present:
        status = EvidenceAvailability.PARTIAL.value
    else:
        status = EvidenceAvailability.UNAVAILABLE.value
    return EvidenceQuality(status, present, missing, tuple(warnings), None, len(present))


__all__ = ["quality_for_sections"]
