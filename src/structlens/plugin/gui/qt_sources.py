"""Pure source-selection helpers for the Qt workflow."""

from __future__ import annotations

import re
from typing import Any

from structlens.core.models import ProteinChain, ProteinStructure, ResidueId

_RESIDUE_TOKEN = re.compile(r"^(?P<chain>[^:]+):(?P<number>-?\d+)(?P<insertion>[A-Za-z]?)$")


def combo_data(combo: Any) -> str | None:
    """Read one user selection as normalized text without scientific logic."""

    value = combo.currentData()
    return str(value) if value is not None else None


def selected_chain(
    structure: ProteinStructure | None,
    combo: Any,
    model_combo: Any | None = None,
) -> ProteinChain | None:
    """Resolve the exact chain/model selected by the user."""

    if structure is None:
        return None
    chain_id = combo_data(combo)
    model_id = combo_data(model_combo) if model_combo is not None else None
    for chain in structure.chains:
        if (chain_id is None or chain.chain_id == chain_id) and (model_id is None or str(chain.model_id) == model_id):
            return chain
    return None


def structure_meta(structure: ProteinStructure) -> str:
    """Format source metadata for the project page."""

    count = sum(len(chain.residue_records or chain.residues) for chain in structure.chains)
    return f"{len(structure.chains)} chain(s) · {count} residues · {structure.structure_id}"


def find_residue(chain: ProteinChain, token: str) -> ResidueId | None:
    """Resolve a displayed chain/number/insertion token to its source identity."""

    match = _RESIDUE_TOKEN.match(token)
    if match is None:
        return None
    chain_id = match.group("chain")
    auth_seq_id = match.group("number")
    insertion_code = match.group("insertion") or None
    for residue in chain.residue_records:
        residue_id = residue.residue_id
        if (
            residue_id.chain_id == chain_id
            and residue_id.auth_seq_id == auth_seq_id
            and residue_id.insertion_code == insertion_code
        ):
            return residue_id
    for residue_id in chain.residues:
        if (
            residue_id.chain_id == chain_id
            and residue_id.auth_seq_id == auth_seq_id
            and residue_id.insertion_code == insertion_code
        ):
            return residue_id
    return None


__all__ = ["combo_data", "find_residue", "selected_chain", "structure_meta"]
