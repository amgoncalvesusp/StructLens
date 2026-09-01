"""Versioned validation for the canonical analysis-report JSON contract.

The report model is intentionally richer than a JSON schema can express.  The
schema therefore checks the wire shape and this module performs the remaining
finite-value and compatibility checks before a payload is persisted/exported.
"""

from __future__ import annotations

import importlib
import json
import math
import re
from collections.abc import Mapping
from importlib import resources
from typing import Any, cast

JSONValue = str | int | float | bool | None | list["JSONValue"] | dict[str, "JSONValue"]
SUPPORTED_REPORT_SCHEMA_VERSION = "4.0"
_RESOURCE_NAME = "analysis-report-v1.json"
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_REPORT_FIELDS = frozenset(
    {
        "schema_version",
        "comparison_mode",
        "reference_selection",
        "target_selection",
        "input_quality",
        "analysis",
        "msa",
        "interactions",
        "sites",
        "distance_map",
        "displacement_vectors",
        "evidence_cards",
        "site_definitions",
        "availability",
        "diagnostics",
        "provenance",
        "pockets",
        "report_id",
    }
)
_SELECTION_FIELDS = frozenset(
    {
        "content_id",
        "format",
        "model_id",
        "author_chain_ids",
        "label_chain_ids",
        "altloc_policy",
        "assembly_scope",
        "chain_locators",
        "selection_id",
    }
)
_CHAIN_LOCATOR_FIELDS = frozenset({"author_chain_id", "label_chain_id", "entity_id"})
_AVAILABILITY_FIELDS = frozenset(
    {
        "input_quality",
        "analysis",
        "msa",
        "interactions",
        "sites",
        "distance_map",
        "displacement_vectors",
        "evidence_cards",
        "pockets",
    }
)
_QUALITY_FIELDS = frozenset({"availability", "diagnostics", "counts", "settings", "provenance"})
_QC_SETTINGS_FIELDS = frozenset(
    {
        "altloc_occupancy_tolerance",
        "b_factor_minimum",
        "occupancy_minimum",
        "occupancy_maximum",
        "min_cn_distance_angstrom",
        "max_cn_distance_angstrom",
        "min_ca_distance_angstrom",
        "max_ca_distance_angstrom",
        "heavy_atom_overlap_tolerance_angstrom",
    }
)
_DIAGNOSTIC_FIELDS = frozenset({"code", "severity", "message", "source_id", "atom_id", "residue_id", "remediation"})
_PROVENANCE_FIELDS = frozenset(
    {
        "method_id",
        "method_version",
        "parameters",
        "units",
        "backend_versions",
        "input_hashes",
        "analyzed_representation",
        "artifact_id",
    }
)
_POCKET_FIELDS = frozenset(
    {
        "detections",
        "volumes",
        "matching",
        "comparisons",
        "lining_residues",
        "concordance",
        "diagnostics",
        "provenance",
        "units",
        "availability",
    }
)
# fmt: off
_ANALYSIS_FIELDS = frozenset({"reference_id", "target_id", "correspondences", "mutations", "sequence_identity", "sequence_coverage", "alignment_decision", "sequence_similarity", "strict_rmsd_angstrom", "refined_rmsd_angstrom", "mapped_residue_count", "refined_residue_count", "excluded_alignment_indices", "tm_score", "legacy_provenance", "transform", "method_provenance"})
_CORRESPONDENCE_FIELDS = frozenset({"alignment_index", "reference", "target", "reference_one_letter", "target_one_letter", "status", "sequence_score", "ca_displacement_angstrom", "backbone_rmsd_angstrom", "sidechain_rmsd_angstrom", "all_heavy_atom_rmsd_angstrom", "is_outlier", "is_key_residue", "mapping_source", "mapping_locked"})
_MUTATION_FIELDS = frozenset({"alignment_index", "kind", "reference", "target", "reference_aa", "target_aa", "reference_label", "target_label", "canonical_notation", "blosum62_score", "grantham_distance", "physicochemical_class"})
_TRANSFORM_FIELDS = frozenset({"rotation", "translation"})
_VECTOR_FIELDS = frozenset({"reference_position", "reference_residue", "target_residue", "start_xyz", "end_xyz", "vector_xyz", "magnitude_angstrom"})
_SITE_DEFINITION_FIELDS = frozenset({"site_id", "name", "mode", "reference_residues", "center_residue", "ligand_id", "radius_angstrom"})
_DISTANCE_MAP_FIELDS = frozenset({"reference_positions", "reference_distances_angstrom", "target_distances_angstrom", "delta_angstrom", "valid_mask"})
_EVIDENCE_CARD_FIELDS = frozenset({"reference_residue", "target_id", "residue_ref", "sequence", "structure", "interactions", "site", "quality", "schema_version", "provenance", "pocket_concordance"})
_POCKET_DETECTION_FIELDS = frozenset({"role", "availability", "detection", "candidates", "settings", "diagnostics", "counts", "provenance"})
_POCKET_RESULT_FIELDS = frozenset({"availability", "coarse_voxel_count", "fine_voxel_count", "coarse_volume_angstrom3", "fine_volume_angstrom3", "coarse", "fine", "sensitivity", "absolute_sensitivity_angstrom3", "relative_sensitivity", "rotation_error_bound_angstrom3", "grid", "units", "provenance_parameters", "provenance", "settings", "candidate_id", "compatibility_signature", "diagnostics"})
_POCKET_VOLUME_FIELDS = frozenset({"role", "availability", "result"})
_POCKET_MATCHING_FIELDS = frozenset({"availability", "matches", "diagnostics", "settings", "transform", "optimal_score", "second_best_score"})
_POCKET_MATCH_FIELDS = frozenset({"state", "status", "reference_candidate", "target_candidate", "lining_jaccard", "lining_overlap_fraction", "centroid_distance_angstrom", "score", "alternative_score", "ambiguity_margin", "assignment_score", "alternative_assignment_score", "units", "diagnostics"})
_POCKET_CANDIDATE_FIELDS = frozenset({"alpha_spheres", "centroid_xyz", "lining_residues", "touching_atom_ids", "candidate_id", "lineage"})
_POCKET_SPHERE_FIELDS = frozenset({"center_xyz", "radius_angstrom", "touching_atom_ids", "lining_residues", "source_simplex_atom_ids", "sphere_id"})
_POCKET_LINEAGE_FIELDS = frozenset({"source_content_id", "selection_id"})
_POCKET_COMPARISON_FIELDS = frozenset({"availability", "match", "volume", "surface_delta_angstrom2", "relative_surface_delta_fraction", "reference_surface_area_angstrom2", "target_surface_area_angstrom2", "reference_surface_method", "target_surface_method", "reference_surface_provenance", "target_surface_provenance", "reference_surface_units", "target_surface_units", "lining_residue_conserved", "lining_residue_gains", "lining_residue_losses", "associated_mutations", "interaction_changes", "local_displacement_angstrom", "ca_displacement_angstrom", "sidechain_displacement_angstrom", "local_displacements", "qc_availability", "qc_diagnostics", "diagnostics", "units"})
_POCKET_VOLUME_COMPARISON_FIELDS = frozenset({"availability", "delta_angstrom3", "relative_delta_fraction", "relative_delta", "reference_volume_angstrom3", "target_volume_angstrom3", "units", "compatibility_signature", "reference_sensitivity", "target_sensitivity", "provenance", "diagnostics"})
_POCKET_LINING_FIELDS = frozenset({"role", "candidate_id", "residues", "availability", "diagnostics", "units"})
_POCKET_CHANNEL_FIELDS = frozenset({"availability", "measure", "units", "diagnostics", "provenance"})
_RESIDUE_FIELDS = frozenset({"structure_id", "model_id", "chain_id", "auth_seq_id", "insertion_code", "residue_name"})
_LOCATOR_FIELDS = frozenset({"author_chain_id", "label_chain_id", "entity_id"})
_MSA_FIELDS = frozenset({"sequences", "aligned_rows", "columns", "reference_structure_id", "algorithm", "provenance"})
_MSA_SEQUENCE_FIELDS = frozenset({"structure_id", "chain_id", "sequence", "source", "residues"})
_SEQUENCE_REF_FIELDS = frozenset({"sequence_index", "one_letter", "residue_id"})
_MSA_COLUMN_FIELDS = frozenset({"index", "reference_label", "reference_residue", "cells", "non_gap_count", "gap_fraction", "ambiguous_fraction", "conservation_score", "entropy_bits"})
_MSA_CELL_FIELDS = frozenset({"structure_id", "alignment_column", "residue", "character"})
_INTERACTION_RECORD_FIELDS = frozenset({"structure_id", "interaction_type", "residue_a", "residue_b", "atom_a", "atom_b", "distance_angstrom", "angle_degrees", "ligand_or_metal_id", "evidence_mode"})
_INTERACTION_KEY_FIELDS = frozenset(
    {"interaction_type", "reference_position_a", "reference_position_b", "external_partner_id"}
)
_INTERACTION_DIFFERENCE_FIELDS = frozenset({"key", "change", "reference_record", "target_record"})
_INTERACTION_EVIDENCE_FIELDS = frozenset({"differences", "reference_interactions", "target_interactions"})
_SITE_METRIC_FIELDS = frozenset({"site_id", "structure_id", "mapped_residue_count", "coverage_fraction", "global_frame_backbone_rmsd_angstrom", "site_fitted_backbone_rmsd_angstrom", "centroid_displacement_angstrom", "radius_of_gyration_angstrom", "atomic_envelope_volume_angstrom3", "sasa_angstrom2", "polar_residue_fraction", "charged_residue_fraction"})
_SITE_EVIDENCE_FIELDS = frozenset({"metrics"})
_SEQUENCE_EVIDENCE_FIELDS = frozenset({"reference_one_letter", "target_one_letter", "alignment_index", "sequence_identity", "conservation_fraction", "entropy_bits", "gap_fraction", "ambiguous_fraction", "source_refs"})
_STRUCTURE_EVIDENCE_FIELDS = frozenset({"ca_displacement_angstrom", "backbone_rmsd_angstrom", "sidechain_rmsd_angstrom", "all_heavy_atom_rmsd_angstrom", "sasa_reference_angstrom2", "sasa_target_angstrom2", "available"})
_EVIDENCE_INTERACTIONS_FIELDS = _INTERACTION_EVIDENCE_FIELDS
_EVIDENCE_QUALITY_FIELDS = frozenset(
    {"overall_status", "available_sections", "unavailable_sections", "warnings", "coverage_fraction", "source_count"}
)
_POCKET_SETTINGS_FIELDS = frozenset({"geometry", "cluster_distance_padding_angstrom", "solvent_grid_spacing_angstrom", "solvent_boundary_margin_angstrom", "minimum_cluster_size", "maximum_candidates", "max_atom_count", "max_estimated_simplices", "max_solvent_grid_cells", "max_clearance_atom_checks", "max_solvent_raster_cells"})
_POCKET_GEOMETRY_FIELDS = frozenset({"minimum_alpha_sphere_radius_angstrom", "maximum_alpha_sphere_radius_angstrom", "probe_radius_angstrom", "lining_contact_slack_angstrom"})
_POCKET_DETECTION_RESULT_FIELDS = frozenset({"availability", "spheres", "diagnostics"})
_POCKET_VOLUME_SETTINGS_FIELDS = frozenset({"coarse_grid_spacing_angstrom", "fine_grid_spacing_angstrom", "boundary_margin_angstrom", "component_exclusion_policy", "max_voxel_count", "voxel_chunk_size", "max_sphere_count", "max_exclusion_atom_count", "max_sphere_voxel_checks", "max_exclusion_neighbor_checks", "radii_version", "grid_phase", "grid_sampling"})
_POCKET_VOLUME_GRID_FIELDS = frozenset({"shape", "coarse_shape", "fine_shape", "origin_xyz", "phase", "sampling"})
_POCKET_VOLUME_BLOCK_FIELDS = frozenset({"voxel_count", "volume_angstrom3", "grid_shape"})
_POCKET_SENSITIVITY_FIELDS = frozenset({"absolute_angstrom3", "relative_fraction"})
_POCKET_MATCHING_SETTINGS_FIELDS = frozenset({"minimum_lining_jaccard", "maximum_centroid_distance_angstrom", "ambiguity_margin", "lining_weight", "centroid_weight", "maximum_candidate_count"})
_POCKET_COMPARISON_LOCAL_DISPLACEMENT_FIELDS = frozenset({"residue", "magnitude_angstrom"})
_AVAILABILITY_VALUES = frozenset(
    {"available", "not_applicable", "not_detected", "invalid_input", "dependency_unavailable", "numerical_failure"}
)
_MAX_STRING_LENGTH = 1_000_000


class AnalysisReportSchemaError(ValueError):
    """Raised when a report is malformed or uses an unsupported schema."""


def load_analysis_report_schema() -> dict[str, Any]:
    """Load the bundled report schema without depending on the working path."""

    try:
        text = resources.files("structlens.resources.schema").joinpath(_RESOURCE_NAME).read_text(encoding="utf-8")
        decoded = json.loads(text)
    except (ImportError, OSError, json.JSONDecodeError) as exc:
        raise AnalysisReportSchemaError("Unable to load the analysis-report schema") from exc
    if not isinstance(decoded, dict):
        raise AnalysisReportSchemaError("Analysis-report schema must be a JSON object")
    return cast(dict[str, Any], decoded)


def validate_analysis_report_payload(payload: Mapping[str, Any] | str | bytes) -> dict[str, Any]:
    """Validate and return a fresh report payload.

    ``jsonschema`` is deliberately optional at runtime.  The bounded checks
    below remain fail-closed when the optional validator is absent.
    """

    if isinstance(payload, bytes) and len(payload) > 100 * 1024 * 1024:
        raise AnalysisReportSchemaError("analysis-report JSON exceeds the maximum supported size")
    if isinstance(payload, str) and len(payload.encode("utf-8")) > 100 * 1024 * 1024:
        raise AnalysisReportSchemaError("analysis-report JSON exceeds the maximum supported size")
    decoded: object = payload
    if isinstance(payload, (str, bytes)):
        try:
            decoded = json.loads(payload)
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            message = getattr(exc, "msg", str(exc))
            raise AnalysisReportSchemaError(f"Invalid analysis-report JSON: {message}") from exc
    if not isinstance(decoded, Mapping):
        raise AnalysisReportSchemaError("Analysis-report JSON must contain an object")
    candidate = _copy_json(decoded, "report")
    version = candidate.get("schema_version")
    if version != SUPPORTED_REPORT_SCHEMA_VERSION:
        raise AnalysisReportSchemaError(
            f"Unsupported analysis-report schema version {version!r}; supported version is {SUPPORTED_REPORT_SCHEMA_VERSION}"
        )
    _assert_finite(candidate, "report")
    unknown = sorted(set(candidate).difference(_REPORT_FIELDS))
    if unknown:
        raise AnalysisReportSchemaError(f"Analysis-report payload contains unknown fields: {', '.join(unknown)}")
    missing = sorted(_REPORT_FIELDS.difference(candidate))
    if missing:
        raise AnalysisReportSchemaError(f"Analysis-report payload is missing required fields: {', '.join(missing)}")
    if not isinstance(candidate["report_id"], str) or not _SHA256.fullmatch(candidate["report_id"]):
        raise AnalysisReportSchemaError("report_id must be a SHA-256 hexadecimal digest")
    try:
        jsonschema = importlib.import_module("jsonschema")
        validator = jsonschema.Draft202012Validator(load_analysis_report_schema())
        errors = sorted(validator.iter_errors(candidate), key=str)
        if errors:
            error = errors[0]
            location = ".".join(str(item) for item in error.path)
            suffix = f" at {location}" if location else ""
            message = error.message
            nested_errors = list(getattr(error, "context", ()))
            while nested_errors:
                nested = nested_errors.pop(0)
                if getattr(nested, "validator", None) == "additionalProperties":
                    message = f"unknown nested field: {nested.message}"
                    break
                if getattr(nested, "validator", None) == "required" and "pocket_concordance" in location:
                    message = f"must contain all seven channels: {nested.message}"
                    break
                nested_errors.extend(getattr(nested, "context", ()))
            raise AnalysisReportSchemaError(f"Analysis-report schema validation failed{suffix}: {message}")
    except ImportError:
        # The structural checks above are the required fallback for minimal
        # installations; the package itself does not require jsonschema.
        pass
    # Validate the same wire contract without relying on jsonschema.  This is
    # intentionally run even when the optional validator is installed so the
    # fallback path remains continuously exercised against the bundled schema.
    from .schema_interpreter import SchemaInterpreterError, validate_instance

    try:
        validate_instance(candidate, load_analysis_report_schema())
    except SchemaInterpreterError as exc:
        raise AnalysisReportSchemaError(f"Analysis-report schema validation failed: {exc}") from exc
    # Keep the public contract identical with and without the optional
    # jsonschema dependency.  In particular, its deliberately generic object
    # branches (MSA, interactions, and evidence measures) are closed here.
    _validate_fallback_shape(candidate)
    return cast(dict[str, Any], candidate)


def validate_report_schema(payload: Mapping[str, Any] | str | bytes) -> dict[str, Any]:
    """Compatibility alias for callers that use the shorter name."""

    return validate_analysis_report_payload(payload)


def _copy_json(value: object, path: str, *, depth: int = 0, nodes: list[int] | None = None) -> Any:
    counter = nodes if nodes is not None else [0]
    counter[0] += 1
    if depth > 64 or counter[0] > 1_000_000:
        raise AnalysisReportSchemaError(f"{path} exceeds report JSON bounds")
    if isinstance(value, str):
        if len(value) > _MAX_STRING_LENGTH:
            raise AnalysisReportSchemaError(f"{path} contains an oversized string")
        return value
    if value is None or isinstance(value, (bool, int, float)):
        return value
    if isinstance(value, Mapping):
        if len(value) > 100_000:
            raise AnalysisReportSchemaError(f"{path} contains too many fields")
        result: dict[str, Any] = {}
        for key, child in value.items():
            if not isinstance(key, str) or not key or len(key) > 1_024:
                raise AnalysisReportSchemaError(f"{path} has an invalid field name")
            result[key] = _copy_json(child, f"{path}.{key}", depth=depth + 1, nodes=counter)
        return result
    if isinstance(value, (list, tuple)):
        if len(value) > 100_000:
            raise AnalysisReportSchemaError(f"{path} contains too many items")
        return [
            _copy_json(child, f"{path}[{index}]", depth=depth + 1, nodes=counter) for index, child in enumerate(value)
        ]
    raise AnalysisReportSchemaError(f"{path} contains unsupported value {type(value).__name__}")


def _assert_finite(value: object, path: str) -> None:
    if isinstance(value, float) and not math.isfinite(value):
        raise AnalysisReportSchemaError(f"{path} must contain only finite JSON numbers")
    if isinstance(value, Mapping):
        for key, child in value.items():
            _assert_finite(child, f"{path}.{key}")
    elif isinstance(value, (list, tuple)):
        for index, child in enumerate(value):
            _assert_finite(child, f"{path}[{index}]")


def _validate_fallback_shape(payload: Mapping[str, Any]) -> None:
    """Apply the schema's important closed-world checks without jsonschema."""

    if payload["comparison_mode"] not in {"pairwise", "reference_vs_many", "all_vs_all"}:
        raise AnalysisReportSchemaError("comparison_mode is not a supported value")
    for name in ("reference_selection", "target_selection"):
        _validate_fallback_selection(payload[name], name)
    quality_payload = _validate_fallback_object(
        payload["input_quality"], "input_quality", {"reference", "target"}, required={"reference", "target"}
    )
    for role in ("reference", "target"):
        quality = _validate_fallback_object(
            quality_payload[role], f"input_quality.{role}", _QUALITY_FIELDS, required=_QUALITY_FIELDS
        )
        _validate_fallback_availability(quality["availability"], f"input_quality.{role}.availability")
        _validate_fallback_object(quality["settings"], f"input_quality.{role}.settings", _QC_SETTINGS_FIELDS)
        _validate_fallback_diagnostics(quality["diagnostics"], f"input_quality.{role}.diagnostics")
        if quality["provenance"] is not None:
            _validate_fallback_provenance(quality["provenance"], f"input_quality.{role}.provenance")
    availability = _validate_fallback_object(
        payload["availability"], "availability", _AVAILABILITY_FIELDS, required=_AVAILABILITY_FIELDS
    )
    for name, value in availability.items():
        _validate_fallback_availability(value, f"availability.{name}")
    for name in ("displacement_vectors", "evidence_cards", "site_definitions", "diagnostics"):
        if not isinstance(payload[name], list):
            raise AnalysisReportSchemaError(f"{name} must be an array")
    for name in ("analysis", "msa", "distance_map", "provenance", "pockets"):
        if payload[name] is not None and not isinstance(payload[name], Mapping):
            raise AnalysisReportSchemaError(f"{name} must be an object or null")
    if payload["analysis"] is not None:
        _validate_fallback_analysis(payload["analysis"])
    if payload["msa"] is not None:
        _validate_fallback_msa(payload["msa"])
    if payload["distance_map"] is not None:
        _validate_fallback_distance_map(payload["distance_map"])
    _validate_fallback_interactions(payload["interactions"])
    _validate_fallback_sites(payload["sites"])
    for index, item in enumerate(payload["displacement_vectors"]):
        vector = _validate_fallback_object(
            item, f"displacement_vectors[{index}]", _VECTOR_FIELDS, required=_VECTOR_FIELDS
        )
        _validate_fallback_residue(vector["reference_residue"], f"displacement_vectors[{index}].reference_residue")
        _validate_fallback_residue(vector["target_residue"], f"displacement_vectors[{index}].target_residue")
    if payload["evidence_cards"]:
        for index, item in enumerate(payload["evidence_cards"]):
            _validate_fallback_evidence_card(item, index)
    for index, item in enumerate(payload["site_definitions"]):
        definition = _validate_fallback_object(
            item, f"site_definitions[{index}]", _SITE_DEFINITION_FIELDS, required=_SITE_DEFINITION_FIELDS
        )
        for residue_index, residue in enumerate(definition["reference_residues"]):
            _validate_fallback_residue(residue, f"site_definitions[{index}].reference_residues[{residue_index}]")
        if definition["center_residue"] is not None:
            _validate_fallback_residue(definition["center_residue"], f"site_definitions[{index}].center_residue")
    if payload["provenance"] is not None:
        _validate_fallback_provenance(payload["provenance"], "provenance")
    if payload["pockets"] is not None:
        _validate_fallback_pockets(payload["pockets"])
    _validate_fallback_diagnostics(payload["diagnostics"], "diagnostics")
    for name in ("interactions", "sites"):
        if payload[name] is not None and not isinstance(payload[name], (list, Mapping)):
            raise AnalysisReportSchemaError(f"{name} must be an array, object, or null")


def _validate_fallback_selection(value: object, name: str) -> None:
    if not isinstance(value, Mapping):
        raise AnalysisReportSchemaError(f"{name} must be an object")
    unknown = sorted(set(value).difference(_SELECTION_FIELDS))
    if unknown:
        raise AnalysisReportSchemaError(f"{name} contains unknown fields: {', '.join(unknown)}")
    missing = sorted(_SELECTION_FIELDS.difference(value))
    if missing:
        raise AnalysisReportSchemaError(f"{name} is missing required fields: {', '.join(missing)}")
    for field in ("content_id", "selection_id"):
        if not isinstance(value[field], str) or not _SHA256.fullmatch(value[field]):
            raise AnalysisReportSchemaError(f"{name}.{field} must be a SHA-256 hexadecimal digest")
    if value["format"] not in {"pdb", "mmcif"}:
        raise AnalysisReportSchemaError(f"{name}.format is not supported")
    if not isinstance(value["model_id"], str) or not value["model_id"]:
        raise AnalysisReportSchemaError(f"{name}.model_id must be a non-empty string")
    for field in ("author_chain_ids", "label_chain_ids", "chain_locators"):
        if not isinstance(value[field], list):
            raise AnalysisReportSchemaError(f"{name}.{field} must be an array")
    for index, item in enumerate(value["chain_locators"]):
        locator_name = f"{name}.chain_locators[{index}]"
        locator = _validate_fallback_object(
            item, locator_name, _CHAIN_LOCATOR_FIELDS, required=_CHAIN_LOCATOR_FIELDS
        )
        for field, locator_value in locator.items():
            if locator_value is not None and not isinstance(locator_value, str):
                raise AnalysisReportSchemaError(f"{locator_name}.{field} must be a string or null")
        if locator["entity_id"] == "":
            raise AnalysisReportSchemaError(f"{locator_name}.entity_id must not be empty")
        if all(locator_value is None for locator_value in locator.values()):
            raise AnalysisReportSchemaError(f"{locator_name} must identify at least one chain")


def _validate_fallback_object(
    value: object,
    name: str,
    fields: set[str] | frozenset[str],
    *,
    required: set[str] | frozenset[str] = frozenset(),
) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise AnalysisReportSchemaError(f"{name} must be an object")
    unknown = sorted(set(value).difference(fields))
    if unknown:
        raise AnalysisReportSchemaError(f"{name} contains unknown fields: {', '.join(unknown)}")
    missing = sorted(set(required).difference(value))
    if missing:
        raise AnalysisReportSchemaError(f"{name} is missing required fields: {', '.join(missing)}")
    return value


def _validate_fallback_availability(value: object, name: str) -> None:
    if value not in _AVAILABILITY_VALUES:
        raise AnalysisReportSchemaError(f"{name} is not a supported availability value")


def _validate_fallback_diagnostics(value: object, name: str) -> None:
    if not isinstance(value, list):
        raise AnalysisReportSchemaError(f"{name} must be an array")
    for index, item in enumerate(value):
        _validate_fallback_object(item, f"{name}[{index}]", _DIAGNOSTIC_FIELDS)


def _validate_fallback_provenance(value: object, name: str) -> None:
    _validate_fallback_object(value, name, _PROVENANCE_FIELDS)


def _validate_fallback_residue(value: object, name: str) -> None:
    _validate_fallback_object(value, name, _RESIDUE_FIELDS, required=_RESIDUE_FIELDS)


def _validate_fallback_residue_or_none(value: object, name: str) -> None:
    if value is not None:
        _validate_fallback_residue(value, name)


def _validate_fallback_number_or_none(value: object, name: str) -> None:
    if value is not None and (not isinstance(value, (int, float)) or isinstance(value, bool)):
        raise AnalysisReportSchemaError(f"{name} must be a number or null")


def _validate_fallback_msa(value: object) -> None:
    msa = _validate_fallback_object(value, "msa", _MSA_FIELDS, required=_MSA_FIELDS)
    for index, raw_sequence in enumerate(msa["sequences"]):
        sequence = _validate_fallback_object(
            raw_sequence, f"msa.sequences[{index}]", _MSA_SEQUENCE_FIELDS, required=_MSA_SEQUENCE_FIELDS
        )
        for residue_index, raw_ref in enumerate(sequence["residues"]):
            ref = _validate_fallback_object(
                raw_ref,
                f"msa.sequences[{index}].residues[{residue_index}]",
                _SEQUENCE_REF_FIELDS,
                required=_SEQUENCE_REF_FIELDS,
            )
            _validate_fallback_residue_or_none(
                ref["residue_id"], f"msa.sequences[{index}].residues[{residue_index}].residue_id"
            )
    for index, raw_row in enumerate(msa["aligned_rows"]):
        row = _validate_fallback_array(raw_row, f"msa.aligned_rows[{index}]")
        if len(row) != 2:
            raise AnalysisReportSchemaError(f"msa.aligned_rows[{index}] must contain two values")
    for index, raw_column in enumerate(msa["columns"]):
        column = _validate_fallback_object(
            raw_column, f"msa.columns[{index}]", _MSA_COLUMN_FIELDS, required=_MSA_COLUMN_FIELDS
        )
        _validate_fallback_number_or_none(column["conservation_score"], f"msa.columns[{index}].conservation_score")
        _validate_fallback_number_or_none(column["entropy_bits"], f"msa.columns[{index}].entropy_bits")
        _validate_fallback_residue_or_none(column["reference_residue"], f"msa.columns[{index}].reference_residue")
        for cell_index, raw_cell in enumerate(column["cells"]):
            cell = _validate_fallback_object(
                raw_cell,
                f"msa.columns[{index}].cells[{cell_index}]",
                _MSA_CELL_FIELDS,
                required=_MSA_CELL_FIELDS,
            )
            if cell["residue"] is not None:
                residue = _validate_fallback_object(
                    cell["residue"],
                    f"msa.columns[{index}].cells[{cell_index}].residue",
                    _SEQUENCE_REF_FIELDS,
                    required=_SEQUENCE_REF_FIELDS,
                )
                _validate_fallback_residue_or_none(
                    residue["residue_id"], f"msa.columns[{index}].cells[{cell_index}].residue.residue_id"
                )


def _validate_fallback_array(value: object, name: str) -> list[Any]:
    if not isinstance(value, list):
        raise AnalysisReportSchemaError(f"{name} must be an array")
    return value


def _validate_fallback_interaction_record(value: object, name: str) -> None:
    record = _validate_fallback_object(value, name, _INTERACTION_RECORD_FIELDS, required=_INTERACTION_RECORD_FIELDS)
    _validate_fallback_residue(record["residue_a"], f"{name}.residue_a")
    _validate_fallback_residue_or_none(record["residue_b"], f"{name}.residue_b")


def _validate_fallback_interaction_difference(value: object, name: str) -> None:
    difference = _validate_fallback_object(
        value, name, _INTERACTION_DIFFERENCE_FIELDS, required=_INTERACTION_DIFFERENCE_FIELDS
    )
    key = _validate_fallback_object(
        difference["key"], f"{name}.key", _INTERACTION_KEY_FIELDS, required=_INTERACTION_KEY_FIELDS
    )
    if not isinstance(key["interaction_type"], str) or not isinstance(key["reference_position_a"], str):
        raise AnalysisReportSchemaError(f"{name}.key contains invalid values")
    if difference["reference_record"] is not None:
        _validate_fallback_interaction_record(difference["reference_record"], f"{name}.reference_record")
    if difference["target_record"] is not None:
        _validate_fallback_interaction_record(difference["target_record"], f"{name}.target_record")


def _validate_fallback_interaction_item(value: object, name: str) -> None:
    if not isinstance(value, Mapping):
        raise AnalysisReportSchemaError(f"{name} must be an object")
    if "key" in value:
        _validate_fallback_interaction_difference(value, name)
    else:
        _validate_fallback_interaction_record(value, name)


def _validate_fallback_interactions(value: object) -> None:
    if value is None:
        return
    if isinstance(value, list):
        for index, item in enumerate(value):
            _validate_fallback_interaction_item(item, f"interactions[{index}]")
        return
    payload = _validate_fallback_object(
        value, "interactions", _INTERACTION_EVIDENCE_FIELDS, required=_INTERACTION_EVIDENCE_FIELDS
    )
    for name in ("differences", "reference_interactions", "target_interactions"):
        values = _validate_fallback_array(payload[name], f"interactions.{name}")
        for index, item in enumerate(values):
            if name == "differences":
                _validate_fallback_interaction_difference(item, f"interactions.{name}[{index}]")
            else:
                _validate_fallback_interaction_record(item, f"interactions.{name}[{index}]")


def _validate_fallback_site_metric(value: object, name: str) -> None:
    _validate_fallback_object(value, name, _SITE_METRIC_FIELDS, required=_SITE_METRIC_FIELDS)


def _validate_fallback_sites(value: object) -> None:
    if value is None:
        return
    if isinstance(value, list):
        for index, item in enumerate(value):
            _validate_fallback_site_metric(item, f"sites[{index}]")
        return
    payload = _validate_fallback_object(value, "sites", _SITE_EVIDENCE_FIELDS, required=_SITE_EVIDENCE_FIELDS)
    for index, item in enumerate(_validate_fallback_array(payload["metrics"], "sites.metrics")):
        _validate_fallback_site_metric(item, f"sites.metrics[{index}]")


def _validate_fallback_distance_map(value: object) -> None:
    distance_map = _validate_fallback_object(value, "distance_map", _DISTANCE_MAP_FIELDS, required=_DISTANCE_MAP_FIELDS)
    for name in (
        "reference_distances_angstrom",
        "target_distances_angstrom",
        "delta_angstrom",
        "valid_mask",
    ):
        for index, row in enumerate(_validate_fallback_array(distance_map[name], f"distance_map.{name}")):
            _validate_fallback_array(row, f"distance_map.{name}[{index}]")


def _validate_fallback_evidence_card(value: object, index: int) -> None:
    name = f"evidence_cards[{index}]"
    card = _validate_fallback_object(value, name, _EVIDENCE_CARD_FIELDS, required=_EVIDENCE_CARD_FIELDS)
    _validate_fallback_residue_or_none(card["reference_residue"], f"{name}.reference_residue")
    residue_ref = _validate_fallback_object(
        card["residue_ref"], f"{name}.residue_ref", _SEQUENCE_REF_FIELDS, required=_SEQUENCE_REF_FIELDS
    )
    _validate_fallback_residue_or_none(residue_ref["residue_id"], f"{name}.residue_ref.residue_id")
    sequence = _validate_fallback_object(
        card["sequence"], f"{name}.sequence", _SEQUENCE_EVIDENCE_FIELDS, required=_SEQUENCE_EVIDENCE_FIELDS
    )
    for ref_index, ref in enumerate(_validate_fallback_array(sequence["source_refs"], f"{name}.sequence.source_refs")):
        source_ref = _validate_fallback_object(
            ref, f"{name}.sequence.source_refs[{ref_index}]", _SEQUENCE_REF_FIELDS, required=_SEQUENCE_REF_FIELDS
        )
        _validate_fallback_residue_or_none(
            source_ref["residue_id"], f"{name}.sequence.source_refs[{ref_index}].residue_id"
        )
    _validate_fallback_object(
        card["structure"], f"{name}.structure", _STRUCTURE_EVIDENCE_FIELDS, required=_STRUCTURE_EVIDENCE_FIELDS
    )
    _validate_fallback_interactions(card["interactions"])
    site = _validate_fallback_object(
        card["site"], f"{name}.site", _SITE_EVIDENCE_FIELDS, required=_SITE_EVIDENCE_FIELDS
    )
    for metric_index, metric in enumerate(_validate_fallback_array(site["metrics"], f"{name}.site.metrics")):
        _validate_fallback_site_metric(metric, f"{name}.site.metrics[{metric_index}]")
    _validate_fallback_object(
        card["quality"], f"{name}.quality", _EVIDENCE_QUALITY_FIELDS, required=_EVIDENCE_QUALITY_FIELDS
    )
    if card["pocket_concordance"] is not None:
        _validate_fallback_concordance(card["pocket_concordance"], f"{name}.pocket_concordance")


def _validate_fallback_concordance(value: object, name: str) -> None:
    from .schema_nested import validate_concordance

    validate_concordance(value, name)


def _validate_fallback_analysis(value: object) -> None:
    from .schema_nested import validate_analysis

    validate_analysis(value)


def _validate_fallback_pockets(value: object) -> None:
    from .schema_nested import validate_pockets

    validate_pockets(value)


def _validate_fallback_candidate(value: object, name: str) -> None:
    candidate = _validate_fallback_object(value, name, _POCKET_CANDIDATE_FIELDS, required=_POCKET_CANDIDATE_FIELDS)
    _validate_fallback_object(
        candidate["lineage"], f"{name}.lineage", _POCKET_LINEAGE_FIELDS, required=_POCKET_LINEAGE_FIELDS
    )
    for index, sphere in enumerate(candidate["alpha_spheres"]):
        _validate_fallback_object(
            sphere, f"{name}.alpha_spheres[{index}]", _POCKET_SPHERE_FIELDS, required=_POCKET_SPHERE_FIELDS
        )


__all__ = [
    "AnalysisReportSchemaError",
    "SUPPORTED_REPORT_SCHEMA_VERSION",
    "load_analysis_report_schema",
    "validate_analysis_report_payload",
    "validate_report_schema",
]
