"""Coordinate and correspondence helpers used by report sections."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np

from structlens.application.difference_map_service import (
    build_displacement_vectors,
    calculate_distance_difference,
)
from structlens.application.dto import AnalysisReportRequest
from structlens.core.difference_maps import ResidueDisplacementVector
from structlens.core.geometry.kabsch import apply_transform
from structlens.core.models import (
    AnalysisResult,
    ProteinChain,
    ResidueCorrespondence,
    ResidueId,
    ResidueRecord,
    StructuralTransform,
)
from structlens.core.msa import AnalysisSequence, SequenceResidueRef
from structlens.core.reports import DistanceMapSnapshot


def analysis_sequence(chain: ProteinChain) -> AnalysisSequence:
    """Build the sequence object consumed by the MSA engine."""

    residues = tuple(
        SequenceResidueRef(index, record.one_letter or "X", record.residue_id)
        for index, record in enumerate(chain.residue_records)
    )
    return AnalysisSequence(chain.structure_id, chain.chain_id, chain.sequence, residues, "structure")


def reference_position(residue: ResidueId) -> str:
    insertion = residue.insertion_code or ""
    return f"{residue.chain_id}:{residue.auth_seq_id}{insertion}"


def position_lookup(correspondences: Sequence[ResidueCorrespondence]) -> dict[ResidueId, str]:
    """Map both sides of a correspondence to the stable reference position."""

    output: dict[ResidueId, str] = {}
    for item in correspondences:
        if item.reference is None:
            continue
        position = reference_position(item.reference)
        output[item.reference] = position
        if item.target is not None:
            output[item.target] = position
    return output


def ca(record: ResidueRecord | None) -> tuple[float, float, float] | None:
    if record is None:
        return None
    coordinate = next((atom.coordinate for atom in record.atoms if atom.name.upper() == "CA"), None)
    if coordinate is None:
        return None
    return (float(coordinate[0]), float(coordinate[1]), float(coordinate[2]))


def coordinate_maps(
    reference: ProteinChain,
    target: ProteinChain,
    correspondences: Sequence[ResidueCorrespondence],
) -> tuple[
    tuple[str, ...],
    dict[str, ResidueId],
    dict[str, ResidueId],
    dict[str, tuple[float, float, float]],
    dict[str, tuple[float, float, float]],
]:
    """Collect matched residue IDs and observed C-alpha coordinates."""

    reference_by_id = {record.residue_id: record for record in reference.residue_records}
    target_by_id = {record.residue_id: record for record in target.residue_records}
    labels: list[str] = []
    reference_ids: dict[str, ResidueId] = {}
    target_ids: dict[str, ResidueId] = {}
    reference_ca: dict[str, tuple[float, float, float]] = {}
    target_ca: dict[str, tuple[float, float, float]] = {}
    for item in correspondences:
        if item.reference is None:
            continue
        label = reference_position(item.reference)
        labels.append(label)
        reference_ids[label] = item.reference
        coordinate = ca(reference_by_id.get(item.reference))
        if coordinate is not None:
            reference_ca[label] = coordinate
        if item.target is None:
            continue
        target_ids[label] = item.target
        coordinate = ca(target_by_id.get(item.target))
        if coordinate is not None:
            target_ca[label] = coordinate
    return tuple(labels), reference_ids, target_ids, reference_ca, target_ca


def distance_map(
    reference: ProteinChain,
    target: ProteinChain,
    correspondences: Sequence[ResidueCorrespondence],
) -> DistanceMapSnapshot:
    labels, _, _, reference_ca, target_ca = coordinate_maps(reference, target, correspondences)
    matrix = calculate_distance_difference(labels, reference_ca, target_ca)
    return DistanceMapSnapshot.from_matrix(matrix)


def displacement_vectors(
    reference: ProteinChain,
    target: ProteinChain,
    result: AnalysisResult,
    request: AnalysisReportRequest,
) -> tuple[ResidueDisplacementVector, ...]:
    if result.transform is None:
        raise ValueError("displacement vectors require a structural transform")
    labels, reference_ids, target_ids, reference_ca, target_ca = coordinate_maps(
        reference,
        target,
        result.correspondences,
    )
    transformed = transform_coordinates(target_ca, result.transform)
    return build_displacement_vectors(
        labels,
        reference_ids,
        target_ids,
        reference_ca,
        transformed,
        minimum_magnitude_angstrom=request.minimum_vector_magnitude_angstrom,
        maximum_vectors=request.maximum_vectors,
    )


def transform_coordinates(
    coordinates: Mapping[str, Sequence[float]],
    transform: StructuralTransform,
) -> dict[str, tuple[float, float, float]]:
    if not coordinates:
        return {}
    labels = tuple(coordinates)
    values = np.asarray([coordinates[label] for label in labels], dtype=np.float64)
    fitted = apply_transform(values, transform.rotation, transform.translation)
    return {
        label: (float(row[0]), float(row[1]), float(row[2]))
        for label, row in zip(labels, fitted, strict=True)
    }


def homogeneous_transform(transform: StructuralTransform | None) -> np.ndarray | None:
    """Convert row-vector transform semantics to the site-service matrix form."""

    if transform is None:
        return None
    matrix = np.eye(4, dtype=np.float64)
    matrix[:3, :3] = np.asarray(transform.rotation, dtype=np.float64).T
    matrix[:3, 3] = np.asarray(transform.translation, dtype=np.float64)
    return matrix


__all__ = [
    "analysis_sequence",
    "ca",
    "coordinate_maps",
    "displacement_vectors",
    "distance_map",
    "homogeneous_transform",
    "position_lookup",
    "reference_position",
    "transform_coordinates",
]
