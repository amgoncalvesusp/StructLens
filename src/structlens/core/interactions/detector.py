"""Descriptive geometry-based interaction detection."""

from __future__ import annotations

import math
from collections.abc import Iterable

from structlens.core.models import AtomRecord, ResidueRecord

from . import InteractionRecord, InteractionType
from .chemistry import (
    METAL_ELEMENTS,
    aromatic_ring_geometries,
    atom_is_acceptor,
    atom_is_donor,
    cationic_group_geometries,
    residue_is_cationic,
    residue_is_hydrophobic,
)
from .thresholds import InteractionThresholds


def _distance(a: AtomRecord, b: AtomRecord) -> float:
    return math.sqrt(sum((left - right) ** 2 for left, right in zip(a.coordinate, b.coordinate, strict=True)))


def _point_distance(first: tuple[float, float, float], second: tuple[float, float, float]) -> float:
    return math.sqrt(sum((left - right) ** 2 for left, right in zip(first, second, strict=True)))


def _angle_between(first: tuple[float, float, float], second: tuple[float, float, float]) -> float | None:
    first_norm = math.sqrt(sum(value * value for value in first))
    second_norm = math.sqrt(sum(value * value for value in second))
    if first_norm <= 1e-12 or second_norm <= 1e-12:
        return None
    cosine = sum(left * right for left, right in zip(first, second, strict=True)) / (first_norm * second_norm)
    cosine = max(-1.0, min(1.0, cosine))
    return math.degrees(math.acos(cosine))


def _normal_angle(first: tuple[float, float, float], second: tuple[float, float, float]) -> float | None:
    """Return the unoriented angle between two planes' normals."""

    angle = _angle_between(first, second)
    return None if angle is None else min(angle, 180.0 - angle)


def _pi_stacking_geometry(
    first: ResidueRecord,
    second: ResidueRecord,
    limits: InteractionThresholds,
) -> tuple[float, float, str, str] | None:
    """Return best ring distance/angle and representative atom names."""

    first_rings = aromatic_ring_geometries(first)
    second_rings = aromatic_ring_geometries(second)
    candidates: list[tuple[float, float, str, str]] = []
    for first_ring in first_rings:
        for second_ring in second_rings:
            distance = _point_distance(first_ring.centroid, second_ring.centroid)
            if distance > limits.pi_centroid_distance_angstrom:
                continue
            full_angle = _angle_between(first_ring.normal, second_ring.normal)
            if full_angle is None:
                continue
            angle = min(full_angle, 180.0 - full_angle)
            parallel = angle <= limits.pi_parallel_angle_tolerance_degrees
            t_shaped = limits.pi_t_shape_min_angle_degrees <= full_angle <= limits.pi_t_shape_max_angle_degrees
            if parallel or t_shaped:
                candidates.append((distance, angle, first_ring.atom_names[0], second_ring.atom_names[0]))
    return min(candidates, key=lambda candidate: (candidate[0], candidate[1])) if candidates else None


def _cation_pi_geometry(
    ring_residue: ResidueRecord,
    cation_residue: ResidueRecord,
    limits: InteractionThresholds,
) -> tuple[float, float, str, str] | None:
    """Return best cation-to-ring geometry and representative atom names."""

    candidates: list[tuple[float, float, str, str]] = []
    for ring in aromatic_ring_geometries(ring_residue):
        for cation in cationic_group_geometries(cation_residue):
            vector_values = tuple(cation.centroid[index] - ring.centroid[index] for index in range(3))
            vector = (vector_values[0], vector_values[1], vector_values[2])
            distance = math.sqrt(sum(value * value for value in vector))
            if distance > limits.cation_pi_distance_angstrom:
                continue
            full_angle = _angle_between(ring.normal, vector)
            if full_angle is None:
                continue
            angle = min(full_angle, 180.0 - full_angle)
            if angle > limits.cation_pi_normal_tolerance_degrees:
                continue
            # Plane normals are unoriented: parallel approaches at 0° and
            # 180° are scientifically equivalent for a cation–π contact.
            candidates.append((distance, min(angle, 180.0 - angle), ring.atom_names[0], cation.atom_names[0]))
    return min(candidates, key=lambda candidate: (candidate[0], candidate[1])) if candidates else None


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
            pi_geometry = _pi_stacking_geometry(first, second, limits)
            if pi_geometry is not None:
                distance, angle, ring_atom_a, ring_atom_b = pi_geometry
                records.append(
                    InteractionRecord(
                        structure_id or first.residue_id.structure_id,
                        InteractionType.PI_STACKING,
                        first.residue_id,
                        second.residue_id,
                        ring_atom_a,
                        ring_atom_b,
                        distance,
                        angle,
                        None,
                        "aromatic_ring_geometry",
                    )
                )
            if residue_is_cationic(first):
                cation_pi = _cation_pi_geometry(second, first, limits)
                if cation_pi is not None:
                    distance, angle, ring_atom, cation_atom = cation_pi
                    records.append(
                        InteractionRecord(
                            structure_id or first.residue_id.structure_id,
                            InteractionType.CATION_PI,
                            first.residue_id,
                            second.residue_id,
                            cation_atom,
                            ring_atom,
                            distance,
                            angle,
                            None,
                            "cation_ring_geometry",
                        )
                    )
            elif residue_is_cationic(second):
                cation_pi = _cation_pi_geometry(first, second, limits)
                if cation_pi is not None:
                    distance, angle, ring_atom, cation_atom = cation_pi
                    records.append(
                        InteractionRecord(
                            structure_id or first.residue_id.structure_id,
                            InteractionType.CATION_PI,
                            first.residue_id,
                            second.residue_id,
                            ring_atom,
                            cation_atom,
                            distance,
                            angle,
                            None,
                            "cation_ring_geometry",
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
