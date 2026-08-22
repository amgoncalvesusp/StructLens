"""Explicit, non-functional physicochemical residue descriptors."""

from __future__ import annotations

from .blosum import CANONICAL_AMINO_ACIDS

_CATEGORIES = {
    "A": "nonpolar",
    "V": "nonpolar",
    "L": "nonpolar",
    "I": "nonpolar",
    "M": "nonpolar",
    "P": "nonpolar",
    "G": "nonpolar",
    "F": "aromatic",
    "W": "aromatic",
    "Y": "aromatic",
    "S": "polar",
    "T": "polar",
    "N": "polar",
    "Q": "polar",
    "C": "polar",
    "K": "positively_charged",
    "R": "positively_charged",
    "H": "positively_charged",
    "D": "negatively_charged",
    "E": "negatively_charged",
}


def amino_acid_category(amino_acid: str | None) -> str | None:
    """Return an explicit residue chemistry category for canonical symbols."""

    if amino_acid is None:
        return None
    symbol = amino_acid.upper()
    if len(symbol) != 1 or symbol not in CANONICAL_AMINO_ACIDS:
        return None
    return _CATEGORIES[symbol]


def classify_substitution(
    reference_aa: str | None, target_aa: str | None
) -> str | None:
    """Classify a substitution as a chemistry descriptor, never an effect claim."""

    reference_category = amino_acid_category(reference_aa)
    target_category = amino_acid_category(target_aa)
    if reference_category is None or target_category is None:
        return None
    if reference_category == target_category:
        return "conservative"
    return "non_conservative"


__all__ = ["amino_acid_category", "classify_substitution"]
