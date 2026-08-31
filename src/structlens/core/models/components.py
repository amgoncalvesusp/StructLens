"""Immutable retained non-polymer components from a structure source."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import Enum

from .residue import ResidueId
from .structure import AtomRecord, _freeze_metadata


class ComponentKind(str, Enum):
    """Closed vocabulary for source components retained beside the polymer."""

    POLYMER_RESIDUE = "polymer_residue"
    LIGAND = "ligand"
    WATER = "water"
    ION = "ion"
    OTHER = "other"


@dataclass(frozen=True, slots=True)
class StructureComponent:
    """A source component with explicit identity and retained atom evidence.

    Classification is supplied by a boundary parser.  This model deliberately
    does not infer that an unknown HETATM is a ligand; callers must use
    :attr:`ComponentKind.OTHER` when source rules cannot establish a type.
    """

    component_id: str
    kind: ComponentKind
    atoms: tuple[AtomRecord, ...] = field(default_factory=tuple)
    metadata: Mapping[str, object] = field(default_factory=dict)
    residue_id: ResidueId | None = None
    model_id: str | None = None
    author_chain_id: str | None = None
    label_chain_id: str | None = None
    entity_id: str | None = None
    residue_name: str | None = None
    auth_seq_id: str | None = None
    insertion_code: str | None = None

    def __post_init__(self) -> None:
        component_id = str(self.component_id).strip()
        if not component_id:
            raise ValueError("component_id must not be empty")
        kind = self.kind
        if not isinstance(kind, ComponentKind):
            try:
                kind = ComponentKind(kind)
            except (TypeError, ValueError) as exc:
                raise ValueError(f"unknown component kind: {self.kind!r}") from exc
            object.__setattr__(self, "kind", kind)
        atoms = tuple(self.atoms)
        if any(not isinstance(atom, AtomRecord) for atom in atoms):
            raise TypeError("atoms must contain AtomRecord values")
        object.__setattr__(self, "component_id", component_id)
        object.__setattr__(self, "atoms", atoms)
        object.__setattr__(self, "metadata", _freeze_metadata(self.metadata))
        if self.residue_id is not None and not isinstance(self.residue_id, ResidueId):
            raise TypeError("residue_id must be a ResidueId or None")
        for field_name in (
            "model_id",
            "author_chain_id",
            "label_chain_id",
            "entity_id",
            "residue_name",
            "auth_seq_id",
            "insertion_code",
        ):
            value = getattr(self, field_name)
            if value is not None:
                value = str(value).strip()
                if not value and field_name not in {"author_chain_id", "label_chain_id"}:
                    raise ValueError(f"{field_name} must not be empty")
                object.__setattr__(self, field_name, value)

    @property
    def identity(self) -> str:
        """Stable source identity retained by this component."""

        return self.component_id

    @property
    def is_ligand(self) -> bool:
        """Whether the parser explicitly classified this component as ligand."""

        return self.kind is ComponentKind.LIGAND


__all__ = ["ComponentKind", "StructureComponent"]
