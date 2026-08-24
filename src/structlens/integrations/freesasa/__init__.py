"""Optional FreeSASA boundary; scientific callers receive explicit None when unavailable."""

from .adapter import FreeSASAAdapter, calculate_sasa

__all__ = ["FreeSASAAdapter", "calculate_sasa"]
