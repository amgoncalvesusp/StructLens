"""Stable residue identity and source-numbering value objects."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class ResidueId:
    """A globally meaningful residue locator.

    A residue number alone is not an identity: model, chain, and insertion
    code are all retained so that PDB and mmCIF numbering remains lossless.
    ``residue_name`` is kept in the locator because non-standard residues must
    remain distinguishable from canonical amino acids.
    """

    structure_id: str
    model_id: str
    chain_id: str
    auth_seq_id: str
    insertion_code: str | None
    residue_name: str


@dataclass(frozen=True, slots=True)
class ResidueNumbering:
    """Author and label numbering captured from an mmCIF source."""

    auth_seq_id: str
    label_seq_id: str | None
    insertion_code: str | None


__all__ = ["ResidueId", "ResidueNumbering"]
