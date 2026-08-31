"""Serialization helpers for the immutable report value objects.

Keeping the JSON projection code separate from the report models makes the
contracts easier to review.  These helpers only consume typed core objects;
they do not accept arbitrary mappings or filesystem paths as scientific data.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from typing import TypeAlias, cast

import numpy as np

from structlens.core.difference_maps import ResidueDisplacementVector
from structlens.core.evidence import EvidenceCard, SiteEvidence
from structlens.core.interactions import InteractionDifference, InteractionRecord
from structlens.core.models import MutationEvent, ResidueId, StructuralTransform
from structlens.core.msa import (
    AnalysisSequence,
    MSAColumn,
    MSAResidueCell,
    MultipleSequenceAlignment,
    SequenceResidueRef,
)
from structlens.core.parsing import ChainLocator, InputSelection
from structlens.core.sites import SiteDefinition, SiteMetrics

JSONScalar: TypeAlias = str | int | float | bool | None
JSONValue: TypeAlias = JSONScalar | list["JSONValue"] | dict[str, "JSONValue"]


def canonical_bytes(payload: Mapping[str, JSONValue]) -> bytes:
    """Encode a JSON-ready payload deterministically."""

    return json.dumps(
        payload,
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")


def residue_json(residue: ResidueId | None) -> JSONValue:
    if residue is None:
        return None
    if not isinstance(residue, ResidueId):
        raise TypeError("residue values must be ResidueId or None")
    return {
        "structure_id": residue.structure_id,
        "model_id": residue.model_id,
        "chain_id": residue.chain_id,
        "auth_seq_id": residue.auth_seq_id,
        "insertion_code": residue.insertion_code,
        "residue_name": residue.residue_name,
    }


def _locator_json(locator: ChainLocator) -> JSONValue:
    return {
        "author_chain_id": locator.author_chain_id,
        "label_chain_id": locator.label_chain_id,
        "entity_id": locator.entity_id,
    }


def selection_json(selection: InputSelection) -> dict[str, JSONValue]:
    """Project a selection without leaking its local path."""

    return {
        "content_id": selection.content_id,
        "format": selection.format.value,
        "model_id": selection.model_id,
        "author_chain_ids": list(selection.author_chain_ids),
        "label_chain_ids": list(selection.label_chain_ids),
        "altloc_policy": selection.altloc_policy.value,
        "assembly_scope": selection.assembly_scope.value,
        "chain_locators": [_locator_json(item) for item in selection.chain_locators],
        "selection_id": selection.selection_id,
    }


def mutation_json(event: MutationEvent) -> dict[str, JSONValue]:
    return {
        "alignment_index": event.alignment_index,
        "kind": event.kind.value,
        "reference": residue_json(event.reference),
        "target": residue_json(event.target),
        "reference_aa": event.reference_aa,
        "target_aa": event.target_aa,
        "reference_label": event.reference_label,
        "target_label": event.target_label,
        "canonical_notation": event.canonical_notation,
        "blosum62_score": event.blosum62_score,
        "grantham_distance": event.grantham_distance,
        "physicochemical_class": event.physicochemical_class,
    }


def transform_json(transform: StructuralTransform | None) -> JSONValue:
    if transform is None:
        return None
    return {"rotation": [list(row) for row in transform.rotation], "translation": list(transform.translation)}


def matrix_rows(
    value: object,
    name: str,
    *,
    boolean: bool = False,
) -> tuple[tuple[float, ...] | tuple[bool, ...], ...]:
    array = np.asarray(value, dtype=bool if boolean else np.float64)
    if array.ndim != 2 or array.shape[0] != array.shape[1]:
        raise ValueError(f"{name} must be a square matrix")
    if not boolean and not np.isfinite(array).all():
        raise ValueError(f"{name} must contain finite values")
    if boolean:
        return tuple(tuple(bool(item) for item in row) for row in array.tolist())
    return tuple(tuple(float(item) for item in row) for row in array.tolist())


def copy_msa(msa: MultipleSequenceAlignment) -> MultipleSequenceAlignment:
    """Copy all nested MSA values so a report does not alias caller data."""

    if not isinstance(msa, MultipleSequenceAlignment):
        raise TypeError("msa must be a MultipleSequenceAlignment or None")
    sequences = tuple(
        AnalysisSequence(
            structure_id=row.structure_id,
            chain_id=row.chain_id,
            sequence=row.sequence,
            residues=tuple(
                SequenceResidueRef(item.sequence_index, item.one_letter, item.residue_id) for item in row.residues
            ),
            source=row.source,
        )
        for row in tuple(msa.sequences)
    )
    columns = tuple(
        MSAColumn(
            index=column.index,
            reference_label=column.reference_label,
            reference_residue=column.reference_residue,
            cells=tuple(
                MSAResidueCell(cell.structure_id, cell.alignment_column, cell.residue, cell.character)
                for cell in tuple(column.cells)
            ),
            non_gap_count=column.non_gap_count,
            gap_fraction=column.gap_fraction,
            ambiguous_fraction=column.ambiguous_fraction,
            conservation_score=column.conservation_score,
            entropy_bits=column.entropy_bits,
        )
        for column in tuple(msa.columns)
    )
    return MultipleSequenceAlignment(
        sequences=sequences,
        aligned_rows=tuple((str(item[0]), str(item[1])) for item in tuple(msa.aligned_rows)),
        columns=columns,
        reference_structure_id=msa.reference_structure_id,
        algorithm=msa.algorithm,
        provenance=tuple(msa.provenance),
    )


def interaction_record_json(value: InteractionRecord | None) -> JSONValue:
    if value is None:
        return None
    return {
        "structure_id": value.structure_id,
        "interaction_type": value.interaction_type.value,
        "residue_a": residue_json(value.residue_a),
        "residue_b": residue_json(value.residue_b),
        "atom_a": value.atom_a,
        "atom_b": value.atom_b,
        "distance_angstrom": value.distance_angstrom,
        "angle_degrees": value.angle_degrees,
        "ligand_or_metal_id": value.ligand_or_metal_id,
        "evidence_mode": value.evidence_mode,
    }


def interaction_json(value: InteractionDifference | InteractionRecord) -> JSONValue:
    if isinstance(value, InteractionDifference):
        return {
            "key": {
                "interaction_type": value.key.interaction_type.value,
                "reference_position_a": value.key.reference_position_a,
                "reference_position_b": value.key.reference_position_b,
                "external_partner_id": value.key.external_partner_id,
            },
            "change": value.change.value,
            "reference_record": interaction_record_json(value.reference_record),
            "target_record": interaction_record_json(value.target_record),
        }
    return interaction_record_json(value)


def site_json(value: SiteMetrics) -> dict[str, JSONValue]:
    names = (
        "site_id",
        "structure_id",
        "mapped_residue_count",
        "coverage_fraction",
        "global_frame_backbone_rmsd_angstrom",
        "site_fitted_backbone_rmsd_angstrom",
        "centroid_displacement_angstrom",
        "radius_of_gyration_angstrom",
        "atomic_envelope_volume_angstrom3",
        "sasa_angstrom2",
        "polar_residue_fraction",
        "charged_residue_fraction",
    )
    return {name: cast(JSONValue, getattr(value, name)) for name in names}


def site_definition_json(value: SiteDefinition) -> dict[str, JSONValue]:
    return {
        "site_id": value.site_id,
        "name": value.name,
        "mode": value.mode.value,
        "reference_residues": [residue_json(item) for item in value.reference_residues],
        "center_residue": residue_json(value.center_residue),
        "ligand_id": value.ligand_id,
        "radius_angstrom": value.radius_angstrom,
    }


def sites_present(value: tuple[SiteMetrics, ...] | SiteEvidence | None) -> bool:
    if isinstance(value, SiteEvidence):
        return bool(value.metrics)
    return bool(value)


def evidence_card_json(value: EvidenceCard) -> JSONValue:
    """Project an evidence card through its typed nested contracts."""

    return {
        "reference_residue": residue_json(value.residue_id),
        "target_id": value.target_id,
        "residue_ref": {
            "sequence_index": value.residue_ref.sequence_index,
            "one_letter": value.residue_ref.one_letter,
            "residue_id": residue_json(value.residue_ref.residue_id),
        },
        "sequence": {
            "reference_one_letter": value.sequence.reference_one_letter,
            "target_one_letter": value.sequence.target_one_letter,
            "alignment_index": value.sequence.alignment_index,
            "sequence_identity": value.sequence.sequence_identity,
            "conservation_fraction": value.sequence.conservation_fraction,
            "entropy_bits": value.sequence.entropy_bits,
            "gap_fraction": value.sequence.gap_fraction,
            "ambiguous_fraction": value.sequence.ambiguous_fraction,
            "source_refs": [
                {
                    "sequence_index": ref.sequence_index,
                    "one_letter": ref.one_letter,
                    "residue_id": residue_json(ref.residue_id),
                }
                for ref in value.sequence.source_refs
            ],
        },
        "structure": {
            name: cast(JSONValue, getattr(value.structure, name))
            for name in (
                "ca_displacement_angstrom",
                "backbone_rmsd_angstrom",
                "sidechain_rmsd_angstrom",
                "all_heavy_atom_rmsd_angstrom",
                "sasa_reference_angstrom2",
                "sasa_target_angstrom2",
                "available",
            )
        },
        "interactions": {
            "differences": [interaction_json(item) for item in value.interactions.differences],
            "reference_interactions": [interaction_record_json(item) for item in value.interactions.reference_interactions],
            "target_interactions": [interaction_record_json(item) for item in value.interactions.target_interactions],
        },
        "site": {"metrics": [site_json(item) for item in value.site.metrics]},
        "quality": {
            "overall_status": value.quality.overall_status,
            "available_sections": list(value.quality.available_sections),
            "unavailable_sections": list(value.quality.unavailable_sections),
            "warnings": list(value.quality.warnings),
            "coverage_fraction": value.quality.coverage_fraction,
            "source_count": value.quality.source_count,
        },
        "schema_version": value.schema_version,
        "provenance": list(value.provenance),
    }


def vector_json(value: ResidueDisplacementVector) -> dict[str, JSONValue]:
    return {
        "reference_position": value.reference_position,
        "reference_residue": residue_json(value.reference_residue),
        "target_residue": residue_json(value.target_residue),
        "start_xyz": list(value.start_xyz),
        "end_xyz": list(value.end_xyz),
        "vector_xyz": list(value.vector_xyz),
        "magnitude_angstrom": value.magnitude_angstrom,
    }


def msa_json(msa: MultipleSequenceAlignment | None) -> JSONValue:
    if msa is None:
        return None
    return {
        "sequences": [
            {
                "structure_id": row.structure_id,
                "chain_id": row.chain_id,
                "sequence": row.sequence,
                "source": row.source,
                "residues": [
                    {
                        "sequence_index": item.sequence_index,
                        "one_letter": item.one_letter,
                        "residue_id": residue_json(item.residue_id),
                    }
                    for item in row.residues
                ],
            }
            for row in msa.sequences
        ],
        "aligned_rows": [list(row) for row in msa.aligned_rows],
        "columns": [
            {
                "index": column.index,
                "reference_label": column.reference_label,
                "reference_residue": residue_json(column.reference_residue),
                "cells": [
                    {
                        "structure_id": cell.structure_id,
                        "alignment_column": cell.alignment_column,
                        "residue": (
                            {
                                "sequence_index": cell.residue.sequence_index,
                                "one_letter": cell.residue.one_letter,
                                "residue_id": residue_json(cell.residue.residue_id),
                            }
                            if cell.residue is not None
                            else None
                        ),
                        "character": cell.character,
                    }
                    for cell in column.cells
                ],
                "non_gap_count": column.non_gap_count,
                "gap_fraction": column.gap_fraction,
                "ambiguous_fraction": column.ambiguous_fraction,
                "conservation_score": column.conservation_score,
                "entropy_bits": column.entropy_bits,
            }
            for column in msa.columns
        ],
        "reference_structure_id": msa.reference_structure_id,
        "algorithm": msa.algorithm,
        "provenance": list(msa.provenance),
    }


__all__ = [
    "JSONScalar",
    "JSONValue",
    "canonical_bytes",
    "copy_msa",
    "evidence_card_json",
    "interaction_json",
    "interaction_record_json",
    "matrix_rows",
    "msa_json",
    "mutation_json",
    "residue_json",
    "selection_json",
    "site_definition_json",
    "site_json",
    "sites_present",
    "transform_json",
    "vector_json",
]
