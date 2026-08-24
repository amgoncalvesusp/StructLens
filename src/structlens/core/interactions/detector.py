"""Descriptive geometry-based interaction detection."""

from __future__ import annotations

import math
from collections.abc import Iterable

from structlens.core.models import AtomRecord, ResidueRecord

from . import InteractionRecord, InteractionType
from .chemistry import (
    METAL_ELEMENTS,
    atom_is_acceptor,
    atom_is_donor,
    residue_is_cationic,
    residue_is_hydrophobic,
)
from .thresholds import InteractionThresholds


def _distance(a: AtomRecord, b: AtomRecord) -> float:
    return math.sqrt(sum((left - right) ** 2 for left, right in zip(a.coordinate, b.coordinate, strict=True)))


def detect_interactions(
    residues: Iterable[ResidueRecord],
    thresholds: InteractionThresholds | None = None,
    *,
    structure_id: str | None = None,
) -> tuple[InteractionRecord, ...]:
    """Detect pairwise contacts with explicit, descriptive evidence modes."""

    values = tuple(residues)
    limits = thresholds or InteractionThresholds()
    records: list[InteractionRecord] = []
    for index, first in enumerate(values):
        for second in values[index + 1 :]:
            best: tuple[InteractionType, AtomRecord, AtomRecord, float, str] | None = None
            for atom_a in first.atoms:
                for atom_b in second.atoms:
                    distance = _distance(atom_a, atom_b)
                    if distance <= limits.hbond_distance_angstrom and (
                        (atom_is_donor(atom_a) and atom_is_acceptor(atom_b))
                        or (atom_is_donor(atom_b) and atom_is_acceptor(atom_a))
                    ):
                        best = (InteractionType.HBOND_GEOMETRIC, atom_a, atom_b, distance, "heavy_atom_geometry")
                        break
                    if distance <= limits.hydrophobic_distance_angstrom and residue_is_hydrophobic(first) and residue_is_hydrophobic(second):
                        best = (InteractionType.HYDROPHOBIC, atom_a, atom_b, distance, "heavy_atom_geometry")
                        break
                if best:
                    break
            if best is None and residue_is_cationic(first) and residue_is_cationic(second):
                # Cation-cation is not a supported salt bridge; leave it out.
                pass
            elif best is None:
                # Residue-level salt-bridge classification uses the closest
                # heavy atoms and explicit acidic/basic residue names.
                charged_pair = {first.residue_name.upper(), second.residue_name.upper()}
                if charged_pair & {"ARG", "LYS", "HIS"} and charged_pair & {"ASP", "GLU"}:
                    closest = min((_distance(a, b), a, b) for a in first.atoms for b in second.atoms)
                    if closest[0] <= limits.salt_bridge_distance_angstrom:
                        best = (InteractionType.SALT_BRIDGE, closest[1], closest[2], closest[0], "heavy_atom_geometry")
            if best is not None:
                kind, atom_a, atom_b, distance, evidence = best
                records.append(
                    InteractionRecord(
                        structure_id or first.residue_id.structure_id,
                        kind,
                        first.residue_id,
                        second.residue_id,
                        atom_a.name,
                        atom_b.name,
                        distance,
                        None,
                        None,
                        evidence,
                    )
                )
    # Explicit metal contacts are independent of residue-residue chemistry.
    for residue in values:
        for atom in residue.atoms:
            if atom.element.upper() not in METAL_ELEMENTS:
                continue
            for partner in values:
                if partner is residue:
                    continue
                for other in partner.atoms:
                    distance = _distance(atom, other)
                    if distance <= limits.metal_distance_angstrom:
                        records.append(
                            InteractionRecord(
                                structure_id or residue.residue_id.structure_id,
                                InteractionType.METAL_CONTACT,
                                residue.residue_id,
                                partner.residue_id,
                                atom.name,
                                other.name,
                                distance,
                                None,
                                atom.name,
                                "heavy_atom_geometry",
                            )
                        )
                        break
    return tuple(records)


__all__ = ["detect_interactions"]
