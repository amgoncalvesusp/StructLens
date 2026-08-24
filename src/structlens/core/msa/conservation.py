"""Alignment conservation and frequency calculations."""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Iterable, Sequence

from . import MSAColumn

CANONICAL_AMINO_ACIDS = tuple("ACDEFGHIKLMNPQRSTVWY")
AMBIGUOUS_AMINO_ACIDS = frozenset("XBZJ")


def column_statistics(characters: Iterable[str]) -> tuple[float | None, float | None, float, float, dict[str, float]]:
    values = tuple(character.upper() for character in characters)
    if not values:
        return None, None, 1.0, 0.0, {}
    canonical = [value for value in values if value in CANONICAL_AMINO_ACIDS]
    frequencies = {amino: count / len(canonical) for amino, count in Counter(canonical).items()} if canonical else {}
    entropy = -sum(value * math.log2(value) for value in frequencies.values()) if len(canonical) >= 2 else None
    conservation = 1.0 - entropy / math.log2(20) if entropy is not None else None
    gap_fraction = sum(value in {"-", "."} for value in values) / len(values)
    ambiguous_fraction = sum(value in AMBIGUOUS_AMINO_ACIDS for value in values) / len(values)
    return conservation, entropy, gap_fraction, ambiguous_fraction, frequencies


def conservation_profile(columns: Sequence[MSAColumn]) -> tuple[float | None, ...]:
    return tuple(column.conservation_score for column in columns)


def amino_acid_frequency_matrix(columns: Sequence[MSAColumn]) -> tuple[dict[str, float], ...]:
    return tuple(column_statistics(cell.character for cell in column.cells)[4] for column in columns)


def sequence_logo_heights(columns: Sequence[MSAColumn]) -> tuple[dict[str, float], ...]:
    """Return p_i * (log2(20) - H) with gaps/ambiguity excluded."""
    output: list[dict[str, float]] = []
    for column in columns:
        conservation, entropy, _, _, frequencies = column_statistics(cell.character for cell in column.cells)
        information = math.log2(20) - entropy if entropy is not None else None
        output.append({key: value * information for key, value in frequencies.items()} if information is not None else {})
    return tuple(output)


__all__ = ["AMBIGUOUS_AMINO_ACIDS", "CANONICAL_AMINO_ACIDS", "amino_acid_frequency_matrix", "column_statistics", "conservation_profile", "sequence_logo_heights"]
