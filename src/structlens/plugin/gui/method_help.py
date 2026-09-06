"""Short, scientifically bounded method explanations for GUI controls."""

from __future__ import annotations

METHOD_HELP = {
    "quality": "Coordinate QC reports parser-preserved diagnostics, counts, and geometry checks. It is descriptive evidence, not a biological effect score.",
    "pockets": "Pocket detection uses validated alpha-sphere geometry on the selected polymer scope. Candidates, thresholds, and diagnostics remain traceable.",
    "pocket_volume": "Volume is a dual-resolution atomic-envelope estimate in Å³. Coarse/fine disagreement is reported as sensitivity; it is not a prediction of biological effect.",
    "pocket_matching": "Candidate matching combines declared lining overlap and transformed centroid distance. Ambiguous and unmatched assignments remain explicit.",
}


def method_explanation(topic: str) -> str:
    """Return a bounded explanation, with a safe fallback for new controls."""

    return METHOD_HELP.get(topic, "This control displays evidence produced by the canonical report and does not recalculate it.")


def pocket_volume_explanation() -> str:
    return METHOD_HELP["pocket_volume"]


__all__ = ["METHOD_HELP", "method_explanation", "pocket_volume_explanation"]
