"""Stable provenance serialization for canonical pairwise reports."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import fields
from importlib.metadata import PackageNotFoundError, version

from structlens.application.dto import AnalysisReportRequest
from structlens.core.models import ResidueId
from structlens.core.provenance import FrozenJSON, MethodProvenance
from structlens.core.sites import SiteDefinition

_METHOD_ID = "structlens.analysis_report"
_METHOD_VERSION = "1"


def report_provenance(request: AnalysisReportRequest) -> MethodProvenance:
    """Build deterministic provenance without serializing local file paths."""

    analysis = request.analysis_settings
    thresholds = request.interaction_thresholds
    parameters: dict[str, FrozenJSON] = {
        "transform_direction": "target_to_reference",
        "coordinate_convention": "row_vector_target_at_rotation_plus_translation",
        "reference_selection_id": request.reference_selection.selection_id,
        "target_selection_id": request.target_selection.selection_id,
        "alignment_mode": analysis.alignment_mode.value,
        "minimum_sequence_identity": analysis.minimum_sequence_identity,
        "minimum_sequence_coverage": analysis.minimum_sequence_coverage,
        "substitution_matrix": analysis.substitution_matrix,
        "gap_open": analysis.gap_open,
        "gap_extend": analysis.gap_extend,
        "refined_rmsd": analysis.refined_rmsd,
        "refinement_cutoff_angstrom": analysis.refinement_cutoff_angstrom,
        "refinement_max_iterations": analysis.refinement_max_iterations,
        "msa_algorithm": request.msa_settings.algorithm,
        "msa_mode": request.msa_settings.mode,
        "site_count": len(request.site_definitions),
        "site_definitions": tuple(site_definition_parameter(item) for item in request.site_definitions),
        "manual_pairs": tuple(
            {
                "reference": residue_parameter(canonical_residue_id(reference, "reference")),
                "target": residue_parameter(canonical_residue_id(target, "target")),
            }
            for reference, target in request.manual_pairs
        ),
        "interaction_scope": "selected_polymer_residue_heavy_atoms",
        "minimum_vector_magnitude_angstrom": request.minimum_vector_magnitude_angstrom,
        "maximum_vectors": request.maximum_vectors,
    }
    parameters.update(
        {f"interaction_{field.name}": float(getattr(thresholds, field.name)) for field in fields(thresholds)}
    )
    units = {
        "refinement_cutoff_angstrom": "angstrom",
        "minimum_vector_magnitude_angstrom": "angstrom",
        "site_definitions.radius_angstrom": "angstrom",
        **{
            f"interaction_{field.name}": "degree" if field.name.endswith("_degrees") else "angstrom"
            for field in fields(thresholds)
        },
    }
    return MethodProvenance(
        _METHOD_ID,
        _METHOD_VERSION,
        parameters=parameters,
        units=units,
        backend_versions={
            "biopython": package_version("biopython"),
            "numpy": package_version("numpy"),
            "scipy": package_version("scipy"),
        },
        input_hashes={
            "reference_raw": request.reference_snapshot.raw_sha256,
            "reference_content": request.reference_snapshot.content_id,
            "target_raw": request.target_snapshot.raw_sha256,
            "target_content": request.target_snapshot.content_id,
        },
        analyzed_representation="asymmetric_unit_pair",
    )


def canonical_residue_id(residue: ResidueId, role: str = "reference") -> ResidueId:
    return ResidueId(
        role,
        residue.model_id,
        residue.chain_id,
        residue.auth_seq_id,
        residue.insertion_code,
        residue.residue_name,
    )


def canonical_optional_residue_id(
    residue: ResidueId | None,
    role: str = "reference",
) -> ResidueId | None:
    return None if residue is None else canonical_residue_id(residue, role)


def site_definition_parameter(definition: SiteDefinition) -> Mapping[str, FrozenJSON]:
    return {
        "site_id": definition.site_id,
        "name": definition.name,
        "mode": definition.mode.value,
        "reference_residues": tuple(
            residue_parameter(canonical_residue_id(item)) for item in definition.reference_residues
        ),
        "center_residue": residue_parameter(canonical_optional_residue_id(definition.center_residue)),
        "ligand_id": definition.ligand_id,
        "radius_angstrom": definition.radius_angstrom,
    }


def residue_parameter(residue: ResidueId | None) -> Mapping[str, FrozenJSON] | None:
    if residue is None:
        return None
    return {
        "structure_id": residue.structure_id,
        "model_id": residue.model_id,
        "chain_id": residue.chain_id,
        "auth_seq_id": residue.auth_seq_id,
        "insertion_code": residue.insertion_code,
        "residue_name": residue.residue_name,
    }


def package_version(package: str) -> str:
    try:
        return version(package)
    except PackageNotFoundError:
        return "unavailable"


__all__ = [
    "canonical_optional_residue_id",
    "canonical_residue_id",
    "package_version",
    "report_provenance",
    "residue_parameter",
    "site_definition_parameter",
]
