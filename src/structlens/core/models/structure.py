"""Minimal normalized structure and chain containers.

The containers deliberately carry source identity and residue order without
coupling the domain layer to a parser or to PyMOL. Parsers can add richer atom
records in later layers while keeping these stable identifiers intact.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from types import MappingProxyType

from .residue import ResidueId, ResidueNumbering


@dataclass(frozen=True, slots=True)
class AtomRecord:
    """An immutable atom record normalized from a structure source.

    Coordinates are accepted as a NumPy-like sequence (including an
    ``ndarray``) at the boundary and stored as a tuple, preventing accidental
    in-place edits from leaking into scientific calculations.
    """

    name: str
    element: str
    coordinate: tuple[float, float, float] | Sequence[float]
    altloc: str | None = None
    occupancy: float | None = None

    def __post_init__(self) -> None:
        values = tuple(float(value) for value in self.coordinate)
        if len(values) != 3:
            raise ValueError("atom coordinate must contain exactly three values")
        object.__setattr__(self, "coordinate", values)


@dataclass(frozen=True, slots=True)
class ResidueRecord:
    """A residue with source numbering and its normalized atom records."""

    residue_id: ResidueId
    numbering: ResidueNumbering
    residue_name: str
    one_letter: str | None
    atoms: tuple[AtomRecord, ...] = field(default_factory=tuple)
    is_standard: bool = True


@dataclass(frozen=True, slots=True)
class ProteinChain:
    """An ordered protein chain in one model of a structure."""

    structure_id: str
    model_id: str
    chain_id: str
    residues: tuple[ResidueId, ...] = field(default_factory=tuple)
    sequence: str = ""
    residue_records: tuple[ResidueRecord, ...] = field(default_factory=tuple)
    source_path: str | None = None
    metadata: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))


@dataclass(frozen=True, slots=True)
class ProteinStructure:
    """A normalized structure containing one or more protein chains."""

    structure_id: str
    chains: tuple[ProteinChain, ...] = field(default_factory=tuple)
    source_path: str | None = None
    metadata: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))


__all__ = [
    "AtomRecord",
    "ProteinChain",
    "ProteinStructure",
    "ResidueRecord",
]
