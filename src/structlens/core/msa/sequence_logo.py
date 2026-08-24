"""Presentation-neutral sequence-logo data."""

from __future__ import annotations

from collections.abc import Sequence

from . import MSAColumn
from .conservation import sequence_logo_heights


def build_sequence_logo(columns: Sequence[MSAColumn]) -> tuple[dict[str, float], ...]:
    return sequence_logo_heights(columns)


__all__ = ["build_sequence_logo"]
