"""Legend metadata for continuous Å scales."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class Legend:
    title: str
    unit: str
    minimum: float
    maximum: float


def displacement_legend(minimum: float, maximum: float) -> Legend:
    return Legend("Cα displacement", "Å", minimum, maximum)


def backbone_rmsd_legend(minimum: float, maximum: float) -> Legend:
    return Legend("Backbone RMSD", "Å", minimum, maximum)


__all__ = ["Legend", "backbone_rmsd_legend", "displacement_legend"]
