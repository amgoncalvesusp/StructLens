"""Reference-relative labels for alignment insertion columns."""

from __future__ import annotations

from structlens.core.models import ResidueId


def format_reference_label(residue: ResidueId | None, previous: str | None = None, insertion_index: int = 0) -> str:
    if residue is not None:
        return f"{residue.chain_id}:{residue.auth_seq_id}{residue.insertion_code or ''}"
    if previous is not None:
        return f"{previous}+{insertion_index}"
    return f"N+{insertion_index}"


__all__ = ["format_reference_label"]
