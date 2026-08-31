"""Deterministic residue classification and alternate-location policy."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from structlens.core.models import ComponentKind

from .models import AltlocPolicy, InputSelection

WATER_NAMES = frozenset({"HOH", "WAT", "DOD", "H2O", "SOL"})
# An unlisted HETATM is ``OTHER`` rather than a guessed drug-like ligand.
# ``KNOWN_LIGANDS`` is intentionally limited to cofactors/nucleotides. Common
# carbohydrate components (NAG/BMA/MAN/GLC/FUC) require CCD/connectivity and
# remain OTHER until that evidence is available.
LIGAND_TABLE_VERSION = "2026-08-30"
KNOWN_LIGANDS = frozenset(
    {
        "ATP",
        "ADP",
        "AMP",
        "GTP",
        "GDP",
        "GMP",
        "HEM",
        "HEC",
        "NAD",
        "NAI",
        "FAD",
        "FMN",
        "SAM",
        "SAH",
        "PLP",
        "COA",
    }
)
ION_NAMES = frozenset(
    {
        "LI",
        "NA",
        "K",
        "RB",
        "CS",
        "MG",
        "CA",
        "SR",
        "BA",
        "ZN",
        "MN",
        "FE",
        "CO",
        "NI",
        "CU",
        "CD",
        "HG",
        "CL",
        "BR",
        "IOD",
        "I",
        "F",
    }
)
MODIFIED_POLYMER_NAMES = frozenset({"MSE"})


def classify_residue(residue_name: str, group: str, atoms: Sequence[Any]) -> ComponentKind:
    """Classify a source residue using only explicit, versioned rules."""

    if group.upper() == "ATOM" or residue_name in MODIFIED_POLYMER_NAMES:
        return ComponentKind.POLYMER_RESIDUE
    if residue_name in WATER_NAMES:
        return ComponentKind.WATER
    if residue_name in ION_NAMES:
        return ComponentKind.ION
    elements = {str(getattr(atom, "element", "")).strip().upper() for atom in atoms}
    if len(atoms) == 1 and elements & ION_NAMES:
        return ComponentKind.ION
    if residue_name in KNOWN_LIGANDS:
        return ComponentKind.LIGAND
    return ComponentKind.OTHER


def selected_atoms(residue: Any, policy: AltlocPolicy) -> list[Any]:
    """Select primary atoms while retaining all choices in component metadata."""

    selected: list[Any] = []
    for atom in residue.child_dict.values():
        if atom.is_disordered() != 2:
            selected.append(atom)
            continue
        choices = list(atom.disordered_get_list())
        if policy is AltlocPolicy.ALL:
            selected.extend(choices)
        elif policy is AltlocPolicy.FIRST:
            selected.append(choices[0])
        else:
            selected.append(min(choices, key=alternate_location_sort_key))
    return selected


def alternate_location_sort_key(atom: Any) -> tuple[float, int, str]:
    occupancy = atom.get_occupancy()
    altloc = str(atom.get_altloc()).strip()
    return (-(float(occupancy) if occupancy is not None else -1.0), 0 if not altloc else 1, altloc)


def altloc_inventory(residue: Any) -> dict[str, tuple[str, ...]]:
    inventory: dict[str, tuple[str, ...]] = {}
    for atom in residue.child_dict.values():
        if atom.is_disordered() == 2:
            locations = tuple(
                sorted(
                    {
                        str(choice.get_altloc()).strip()
                        for choice in atom.disordered_get_list()
                        if clean_optional(str(choice.get_altloc())) is not None
                    }
                )
            )
            if locations:
                inventory[str(atom.get_name()).strip()] = locations
    return inventory


def component_id(
    model_id: str,
    author_chain_id: str,
    label_chain_id: str | None,
    entity_id: str | None,
    auth_seq_id: str,
    insertion_code: str | None,
    residue_name: str,
) -> str:
    return ":".join(
        (
            model_id,
            author_chain_id,
            label_chain_id or "-",
            entity_id or "-",
            auth_seq_id,
            insertion_code or "-",
            residue_name,
        )
    )


def chain_selected(selection: InputSelection, author: str, label: str | None, entity: str | None) -> bool:
    if selection.author_chain_ids and author not in selection.author_chain_ids:
        return False
    if selection.label_chain_ids and label not in selection.label_chain_ids:
        return False
    if selection.chain_locators and not any(
        (locator.author_chain_id is None or locator.author_chain_id == author)
        and (locator.label_chain_id is None or locator.label_chain_id == label)
        and (locator.entity_id is None or locator.entity_id == entity)
        for locator in selection.chain_locators
    ):
        return False
    return True


def clean_optional(value: str) -> str | None:
    value = value.strip().strip("'\"")
    return None if value in {"", ".", "?"} else value


__all__ = [
    "ION_NAMES",
    "KNOWN_LIGANDS",
    "LIGAND_TABLE_VERSION",
    "MODIFIED_POLYMER_NAMES",
    "WATER_NAMES",
    "altloc_inventory",
    "alternate_location_sort_key",
    "chain_selected",
    "classify_residue",
    "component_id",
    "selected_atoms",
]
