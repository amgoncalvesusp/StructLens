"""Safe, namespaced PyMOL selection-string helpers."""

from __future__ import annotations

import re

from structlens.core.models import ResidueId


def selection_name(project_id: str, target_id: str, purpose: str) -> str:
    parts = [project_id, target_id, purpose]
    safe = "_".join(re.sub(r"[^A-Za-z0-9_]+", "_", part).strip("_") for part in parts)
    return f"structlens_{safe.lower()}"


def residue_selection(residue: ResidueId) -> str:
    insertion = residue.insertion_code or ""
    if not re.fullmatch(r"-?\d+", str(residue.auth_seq_id)) or not re.fullmatch(r"[A-Za-z0-9]?", insertion):
        raise ValueError("Unsupported residue identifier for a PyMOL selection")
    return f"chain {quote_identifier(residue.chain_id)} and resi {residue.auth_seq_id}{insertion}"


def quote_identifier(value: str) -> str:
    if not value.strip():
        return '""'
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", value):
        raise ValueError("Unsupported identifier for a PyMOL selection")
    return value


__all__ = ["quote_identifier", "residue_selection", "selection_name"]
