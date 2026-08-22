"""Chemically valid atom-name permutations for symmetric side chains."""

from __future__ import annotations

from collections.abc import Iterable

_SWAP_GROUPS: dict[str, tuple[tuple[str, ...], ...]] = {
    "ASP": (("OD1", "OD2"),),
    "GLU": (("OE1", "OE2"),),
    "ARG": (("NH1", "NH2"),),
    # Aromatic rings can rotate by 180 degrees: both pairs must swap together.
    "PHE": (("CD1", "CD2", "CE1", "CE2"),),
    "TYR": (("CD1", "CD2", "CE1", "CE2"),),
}


def symmetry_permutations(
    residue_name: str, atom_names: Iterable[str]
) -> tuple[dict[str, str], ...]:
    """Return target-name mappings that preserve the named residue chemistry.

    A mapping maps each reference atom name to the target atom name to compare
    with it.  The identity mapping is always included.
    """

    names = frozenset(atom_names)
    identity = {name: name for name in names}
    permutations: list[dict[str, str]] = [identity]
    for group in _SWAP_GROUPS.get(residue_name.upper(), ()):
        if not set(group).issubset(names):
            continue
        swapped = identity.copy()
        if len(group) == 2:
            first, second = group
            swapped[first], swapped[second] = second, first
        else:
            swapped.update({"CD1": "CD2", "CD2": "CD1", "CE1": "CE2", "CE2": "CE1"})
        permutations.append(swapped)
    return tuple(permutations)


__all__ = ["symmetry_permutations"]
