"""Minimal normalized structure and chain containers.

The containers deliberately carry source identity and residue order without
coupling the domain layer to a parser or to PyMOL. Parsers can add richer atom
records in later layers while keeping these stable identifiers intact.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any

from .residue import ResidueId, ResidueNumbering


def _freeze_value(value: Any) -> Any:
    """Recursively copy container metadata into immutable values."""

    if isinstance(value, Mapping):
        frozen: dict[str, Any] = {}
        for key, item in value.items():
            if not isinstance(key, str) or not key.strip():
                raise ValueError("metadata keys must be non-empty strings")
            frozen[key] = _freeze_value(item)
        return MappingProxyType(frozen)
    if isinstance(value, (list, tuple)):
        return tuple(_freeze_value(item) for item in value)
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("metadata float values must be finite")
        return value
    raise TypeError(f"metadata value of type {type(value).__name__} is not JSON-like")


def _freeze_metadata(metadata: Mapping[str, object]) -> Mapping[str, object]:
    if not isinstance(metadata, Mapping):
        raise TypeError("metadata must be a mapping")
    frozen: dict[str, object] = {}
    for key, value in metadata.items():
        if not isinstance(key, str) or not key.strip():
            raise ValueError("metadata keys must be non-empty strings")
        frozen[key] = _freeze_value(value)
    return MappingProxyType(frozen)


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
    b_factor: float | None = None
    source_atom_id: str | None = None
    source_serial: int | str | None = None
    formal_charge: int | None = None

    def __post_init__(self) -> None:
        try:
            values = tuple(float(value) for value in self.coordinate)
        except (TypeError, ValueError) as exc:
            raise ValueError("atom coordinate must contain numeric values") from exc
        if len(values) != 3:
            raise ValueError("atom coordinate must contain exactly three values")
        if any(not math.isfinite(value) for value in values):
            raise ValueError("atom coordinate values must be finite")
        name = str(self.name).strip()
        element = str(self.element).strip().upper()
        if not name:
            raise ValueError("atom name must not be empty")
        if not element:
            raise ValueError("atom element must not be empty")
        if self.altloc is not None:
            altloc = str(self.altloc).strip()
            object.__setattr__(self, "altloc", altloc or None)
        if self.occupancy is not None:
            occupancy = float(self.occupancy)
            if not math.isfinite(occupancy) or not 0.0 <= occupancy <= 1.0:
                raise ValueError("atom occupancy must be finite and between 0 and 1")
            object.__setattr__(self, "occupancy", occupancy)
        if self.b_factor is not None:
            b_factor = float(self.b_factor)
            if not math.isfinite(b_factor) or b_factor < 0.0:
                raise ValueError("atom B-factor must be finite and non-negative")
            object.__setattr__(self, "b_factor", b_factor)
        if self.source_atom_id is not None:
            source_atom_id = str(self.source_atom_id).strip()
            if not source_atom_id:
                raise ValueError("source_atom_id must not be empty")
            object.__setattr__(self, "source_atom_id", source_atom_id)
        if self.source_serial is not None:
            if isinstance(self.source_serial, bool) or not isinstance(self.source_serial, (int, str)):
                raise ValueError("source_serial must be an integer or string")
            source_serial = self.source_serial
            if isinstance(source_serial, str):
                source_serial = source_serial.strip()
                if not source_serial:
                    raise ValueError("source_serial must not be empty")
            object.__setattr__(self, "source_serial", source_serial)
        if self.formal_charge is not None:
            if isinstance(self.formal_charge, bool):
                raise ValueError("formal_charge must be an integer")
            try:
                formal_charge = int(self.formal_charge)
            except (TypeError, ValueError) as exc:
                raise ValueError("formal_charge must be an integer") from exc
            if isinstance(self.formal_charge, float) and formal_charge != self.formal_charge:
                raise ValueError("formal_charge must be an integer")
            object.__setattr__(self, "formal_charge", formal_charge)
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "element", element)
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

    def __post_init__(self) -> None:
        if not isinstance(self.residue_id, ResidueId):
            raise TypeError("residue_id must be a ResidueId")
        if not isinstance(self.numbering, ResidueNumbering):
            raise TypeError("numbering must be a ResidueNumbering")
        if not isinstance(self.residue_name, str):
            raise TypeError("residue_name must be a string")
        residue_name = self.residue_name.strip().upper()
        if not residue_name or residue_name != self.residue_id.residue_name.strip().upper():
            raise ValueError("residue_name must match residue_id.residue_name")
        if str(self.numbering.auth_seq_id) != self.residue_id.auth_seq_id:
            raise ValueError("numbering auth_seq_id must match residue_id")
        if self.numbering.insertion_code != self.residue_id.insertion_code:
            raise ValueError("numbering insertion_code must match residue_id")
        atoms = tuple(self.atoms)
        if any(not isinstance(atom, AtomRecord) for atom in atoms):
            raise TypeError("atoms must contain AtomRecord values")
        object.__setattr__(self, "residue_name", residue_name)
        object.__setattr__(self, "atoms", atoms)


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
    author_chain_id: str | None = None
    label_chain_id: str | None = None
    entity_id: str | None = None

    def __post_init__(self) -> None:
        if self.chain_id is not None:
            object.__setattr__(self, "chain_id", str(self.chain_id).strip())
        residues = tuple(self.residues)
        if any(not isinstance(residue_id, ResidueId) for residue_id in residues):
            raise TypeError("residues must contain ResidueId values")
        residue_records = tuple(self.residue_records)
        if any(not isinstance(record, ResidueRecord) for record in residue_records):
            raise TypeError("residue_records must contain ResidueRecord values")
        if residues and residue_records and tuple(record.residue_id for record in residue_records) != residues:
            raise ValueError("residue IDs and residue_records IDs must correspond")
        for residue_id in residues + tuple(record.residue_id for record in residue_records):
            if (
                residue_id.structure_id != self.structure_id
                or residue_id.model_id != self.model_id
                or residue_id.chain_id != self.chain_id
            ):
                raise ValueError("residue ID is not coherent with structure/model/chain")
        object.__setattr__(self, "residues", residues)
        object.__setattr__(self, "residue_records", residue_records)
        object.__setattr__(self, "metadata", _freeze_metadata(self.metadata))
        for field_name in ("author_chain_id", "label_chain_id"):
            value = getattr(self, field_name)
            if value is not None:
                value = str(value).strip()
                object.__setattr__(self, field_name, value)
        if self.entity_id is not None:
            value = str(self.entity_id).strip()
            if not value:
                raise ValueError("entity_id must not be empty")
            object.__setattr__(self, "entity_id", value)


@dataclass(frozen=True, slots=True)
class ProteinStructure:
    """A normalized structure containing one or more protein chains."""

    structure_id: str
    chains: tuple[ProteinChain, ...] = field(default_factory=tuple)
    source_path: str | None = None
    metadata: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        chains = tuple(self.chains)
        if any(not isinstance(chain, ProteinChain) for chain in chains):
            raise TypeError("chains must contain ProteinChain values")
        object.__setattr__(self, "chains", chains)
        object.__setattr__(self, "metadata", _freeze_metadata(self.metadata))


__all__ = [
    "AtomRecord",
    "ProteinChain",
    "ProteinStructure",
    "ResidueRecord",
]
