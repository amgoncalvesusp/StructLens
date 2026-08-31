"""Known-ligand support for detected geometric pocket candidates."""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, replace

from structlens.core.models import AtomRecord, ComponentKind, ResidueId, StructureComponent

from .models import PocketCandidate

POCKET_LIGAND_RULES_VERSION = "structlens-pocket-known-ligands-1"
POCKET_LIGAND_ATOM_SCOPE = "primary_deposited_heavy_atoms"
_ALTLOC_POLICIES = frozenset({"highest_occupancy", "first"})
_HYDROGEN_ELEMENTS = frozenset({"H", "D", "T"})
_EXCLUDED_LIGAND_RESIDUE_NAMES = frozenset(
    {
        "HOH",
        "WAT",
        "DOD",
        "SO4",
        "PO4",
        "PEG",
        "EDO",
        "GOL",
        "MPD",
        "MES",
        "TRS",
        "BME",
        "ACT",
        "ACE",
        "FMT",
        "CIT",
    }
)


@dataclass(frozen=True, slots=True)
class PocketLigandSupport:
    """Deterministic support metrics for one eligible retained ligand."""

    component_id: str
    residue_name: str
    ligand_center_distance_angstrom: float
    atom_coverage_fraction: float
    lining_residue_overlap_fraction: float
    covered_atom_count: int
    atom_count: int
    atom_scope: str = POCKET_LIGAND_ATOM_SCOPE

    def __post_init__(self) -> None:
        component_id = str(self.component_id).strip()
        residue_name = str(self.residue_name).strip()
        if not component_id:
            raise ValueError("component_id must not be empty")
        if not residue_name:
            raise ValueError("residue_name must not be empty")
        object.__setattr__(self, "component_id", component_id)
        object.__setattr__(self, "residue_name", residue_name)
        for name in ("ligand_center_distance_angstrom", "atom_coverage_fraction", "lining_residue_overlap_fraction"):
            value = float(getattr(self, name))
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
            if value < 0.0:
                raise ValueError(f"{name} must be non-negative")
            if name.endswith("fraction") and value > 1.0:
                raise ValueError(f"{name} must not exceed 1")
            object.__setattr__(self, name, value)
        for name in ("covered_atom_count", "atom_count"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ValueError(f"{name} must be a non-negative integer")
        if self.atom_count == 0:
            raise ValueError("atom_count must be positive")
        if self.covered_atom_count > self.atom_count:
            raise ValueError("covered_atom_count must not exceed atom_count")
        expected_fraction = self.covered_atom_count / self.atom_count
        if not math.isclose(self.atom_coverage_fraction, expected_fraction, rel_tol=1e-9, abs_tol=1e-9):
            raise ValueError("atom_coverage_fraction must match covered_atom_count / atom_count")
        atom_scope = str(self.atom_scope).strip()
        if atom_scope != POCKET_LIGAND_ATOM_SCOPE:
            raise ValueError("atom_scope must identify the primary deposited heavy-atom representation")
        object.__setattr__(self, "atom_scope", atom_scope)

    def to_json(self) -> dict[str, object]:
        return {
            "component_id": self.component_id,
            "residue_name": self.residue_name,
            "ligand_center_distance_angstrom": self.ligand_center_distance_angstrom,
            "atom_coverage_fraction": self.atom_coverage_fraction,
            "lining_residue_overlap_fraction": self.lining_residue_overlap_fraction,
            "covered_atom_count": self.covered_atom_count,
            "atom_count": self.atom_count,
            "atom_scope": self.atom_scope,
        }


def eligible_pocket_ligands(
    components: Sequence[StructureComponent],
    *,
    altloc_policy: str = "highest_occupancy",
) -> tuple[StructureComponent, ...]:
    """Return selected retained ligands eligible for ligand-support evaluation."""

    policy = _altloc_policy(altloc_policy)
    values = tuple(components)
    if any(not isinstance(component, StructureComponent) for component in values):
        raise TypeError("components must contain StructureComponent values")
    eligible: list[StructureComponent] = []
    for component in values:
        if component.kind is not ComponentKind.LIGAND:
            continue
        if component.metadata.get("selected_for_analysis") is not True:
            continue
        residue_name = (component.residue_name or "").strip().upper()
        if not residue_name or residue_name in _EXCLUDED_LIGAND_RESIDUE_NAMES:
            continue
        atoms = _primary_heavy_atoms(component.atoms, policy)
        if not atoms:
            continue
        eligible.append(replace(component, atoms=atoms))
    return tuple(sorted(eligible, key=lambda item: item.component_id))


def measure_ligand_support(
    candidate: PocketCandidate,
    ligand: StructureComponent,
    *,
    ligand_contact_residues: Sequence[ResidueId] = (),
    altloc_policy: str = "highest_occupancy",
) -> PocketLigandSupport:
    """Measure support between one pocket candidate and one known ligand."""

    if not isinstance(candidate, PocketCandidate):
        raise TypeError("candidate must be a PocketCandidate")
    if not isinstance(ligand, StructureComponent):
        raise TypeError("ligand must be a StructureComponent")
    if ligand.kind is not ComponentKind.LIGAND:
        raise ValueError("ligand must have component kind 'ligand'")
    contacts = tuple(ligand_contact_residues)
    if any(not isinstance(residue, ResidueId) for residue in contacts):
        raise TypeError("ligand_contact_residues must contain ResidueId values")
    atoms = _primary_heavy_atoms(ligand.atoms, _altloc_policy(altloc_policy))
    if not atoms:
        raise ValueError("ligand must contain at least one atom")
    atom_count = len(atoms)
    covered_atom_count = 0
    for atom in atoms:
        if any(_is_inside_sphere(atom.coordinate, sphere.center_xyz, sphere.radius_angstrom) for sphere in candidate.alpha_spheres):
            covered_atom_count += 1
    centroid = tuple(sum(atom.coordinate[index] for atom in atoms) / atom_count for index in range(3))
    distance = math.sqrt(sum((candidate.centroid_xyz[index] - centroid[index]) ** 2 for index in range(3)))
    contact_set = set(contacts)
    lining_set = set(candidate.lining_residues)
    union = contact_set | lining_set
    overlap_fraction = len(contact_set & lining_set) / len(union) if union else 0.0
    return PocketLigandSupport(
        component_id=ligand.component_id,
        residue_name=ligand.residue_name or ligand.component_id,
        ligand_center_distance_angstrom=distance,
        atom_coverage_fraction=covered_atom_count / atom_count,
        lining_residue_overlap_fraction=overlap_fraction,
        covered_atom_count=covered_atom_count,
        atom_count=atom_count,
    )


def _altloc_policy(value: str) -> str:
    policy = str(value).strip().lower()
    if policy not in _ALTLOC_POLICIES:
        raise ValueError("altloc_policy must be one of: highest_occupancy, first")
    return policy


def _primary_heavy_atoms(atoms: Sequence[AtomRecord], policy: str) -> tuple[AtomRecord, ...]:
    grouped: dict[str, list[AtomRecord]] = {}
    for atom in atoms:
        grouped.setdefault(atom.name, []).append(atom)
    selected: list[AtomRecord] = []
    for choices in grouped.values():
        if len(choices) == 1 or not any(atom.altloc for atom in choices):
            selected.extend(choices)
        elif policy == "first":
            selected.append(choices[0])
        else:
            selected.append(min(choices, key=_alternate_location_sort_key))
    return tuple(
        atom for atom in selected if atom.element.strip().upper() not in _HYDROGEN_ELEMENTS
    )


def _alternate_location_sort_key(atom: AtomRecord) -> tuple[float, int, str]:
    occupancy = atom.occupancy
    altloc = atom.altloc or ""
    return (-(occupancy if occupancy is not None else -1.0), 0 if not altloc else 1, altloc)


def _is_inside_sphere(
    point: Sequence[float],
    center: Sequence[float],
    radius_angstrom: float,
) -> bool:
    return sum((point[index] - center[index]) ** 2 for index in range(3)) <= radius_angstrom**2


__all__ = [
    "POCKET_LIGAND_RULES_VERSION",
    "POCKET_LIGAND_ATOM_SCOPE",
    "PocketLigandSupport",
    "eligible_pocket_ligands",
    "measure_ligand_support",
]
