"""Closed-world validation for nested analysis-report objects.

Kept separate from the public schema loader so the loader remains small and
the optional jsonschema dependency cannot make nested objects permissive.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from . import schema as _schema


def validate_concordance(value: object, name: str) -> None:
    expected = {
        "geometry",
        "ligand_support",
        "parameter_persistence",
        "volume_sensitivity",
        "match_ambiguity",
        "qc",
        "interactions",
    }
    try:
        concordance = _schema._validate_fallback_object(value, name, expected, required=expected)
    except _schema.AnalysisReportSchemaError as exc:
        if "missing required fields" in str(exc):
            raise _schema.AnalysisReportSchemaError(f"{name} must contain all seven channels") from exc
        raise
    for channel_name, channel_value in concordance.items():
        channel = _schema._validate_fallback_object(
            channel_value,
            f"{name}.{channel_name}",
            _schema._POCKET_CHANNEL_FIELDS,
            required=_schema._POCKET_CHANNEL_FIELDS,
        )
        _schema._validate_fallback_availability(channel["availability"], f"{name}.{channel_name}.availability")
        _schema._validate_fallback_diagnostics(channel["diagnostics"], f"{name}.{channel_name}.diagnostics")


def validate_analysis(value: object) -> None:
    analysis = _schema._validate_fallback_object(
        value, "analysis", _schema._ANALYSIS_FIELDS, required=_schema._ANALYSIS_FIELDS
    )
    for index, item in enumerate(analysis["correspondences"]):
        correspondence = _schema._validate_fallback_object(
            item,
            f"analysis.correspondences[{index}]",
            _schema._CORRESPONDENCE_FIELDS,
            required=_schema._CORRESPONDENCE_FIELDS,
        )
        _schema._validate_fallback_residue_or_none(
            correspondence["reference"], f"analysis.correspondences[{index}].reference"
        )
        _schema._validate_fallback_residue_or_none(
            correspondence["target"], f"analysis.correspondences[{index}].target"
        )
    for index, item in enumerate(analysis["mutations"]):
        mutation = _schema._validate_fallback_object(
            item, f"analysis.mutations[{index}]", _schema._MUTATION_FIELDS, required=_schema._MUTATION_FIELDS
        )
        _schema._validate_fallback_residue_or_none(mutation["reference"], f"analysis.mutations[{index}].reference")
        _schema._validate_fallback_residue_or_none(mutation["target"], f"analysis.mutations[{index}].target")
    if analysis["transform"] is not None:
        transform = _schema._validate_fallback_object(
            analysis["transform"], "analysis.transform", _schema._TRANSFORM_FIELDS, required=_schema._TRANSFORM_FIELDS
        )
        _schema._validate_fallback_array(transform["rotation"], "analysis.transform.rotation")
        for index, row in enumerate(transform["rotation"]):
            if len(_schema._validate_fallback_array(row, f"analysis.transform.rotation[{index}]")) != 3:
                raise _schema.AnalysisReportSchemaError("analysis.transform.rotation rows must contain three values")
        if len(_schema._validate_fallback_array(transform["translation"], "analysis.transform.translation")) != 3:
            raise _schema.AnalysisReportSchemaError("analysis.transform.translation must contain three values")
    if analysis["method_provenance"] is not None:
        _schema._validate_fallback_provenance(analysis["method_provenance"], "analysis.method_provenance")


def _array(value: object, name: str) -> list[Any]:
    return _schema._validate_fallback_array(value, name)


def _obj(
    value: object,
    name: str,
    fields: set[str] | frozenset[str],
    required: set[str] | frozenset[str] = frozenset(),
) -> Mapping[str, Any]:
    return _schema._validate_fallback_object(value, name, fields, required=required)


def _diagnostics(value: object, name: str) -> None:
    _schema._validate_fallback_diagnostics(value, name)


def _provenance(value: object, name: str) -> None:
    _schema._validate_fallback_provenance(value, name)


def _candidate(value: object, name: str) -> None:
    candidate = _obj(value, name, _schema._POCKET_CANDIDATE_FIELDS, _schema._POCKET_CANDIDATE_FIELDS)
    _obj(candidate["lineage"], f"{name}.lineage", _schema._POCKET_LINEAGE_FIELDS, _schema._POCKET_LINEAGE_FIELDS)
    for index, sphere_value in enumerate(_array(candidate["alpha_spheres"], f"{name}.alpha_spheres")):
        sphere = _obj(
            sphere_value, f"{name}.alpha_spheres[{index}]", _schema._POCKET_SPHERE_FIELDS, _schema._POCKET_SPHERE_FIELDS
        )
        for residue_index, residue in enumerate(
            _array(sphere["lining_residues"], f"{name}.alpha_spheres[{index}].lining_residues")
        ):
            _schema._validate_fallback_residue(
                residue, f"{name}.alpha_spheres[{index}].lining_residues[{residue_index}]"
            )


def _candidate_array(value: object, name: str) -> None:
    for index, candidate in enumerate(_array(value, name)):
        _candidate(candidate, f"{name}[{index}]")


def _validate_detection(value: object, name: str) -> None:
    detection = _obj(value, name, _schema._POCKET_DETECTION_RESULT_FIELDS, _schema._POCKET_DETECTION_RESULT_FIELDS)
    _schema._validate_fallback_availability(detection["availability"], f"{name}.availability")
    for index, sphere in enumerate(_array(detection["spheres"], f"{name}.spheres")):
        sphere_value = _obj(
            sphere, f"{name}.spheres[{index}]", _schema._POCKET_SPHERE_FIELDS, _schema._POCKET_SPHERE_FIELDS
        )
        for residue_index, residue in enumerate(
            _array(sphere_value["lining_residues"], f"{name}.spheres[{index}].lining_residues")
        ):
            _schema._validate_fallback_residue(residue, f"{name}.spheres[{index}].lining_residues[{residue_index}]")
    _diagnostics(detection["diagnostics"], f"{name}.diagnostics")


def _validate_volume_result(value: object, name: str) -> None:
    result = _obj(value, name, _schema._POCKET_RESULT_FIELDS, _schema._POCKET_RESULT_FIELDS)
    _schema._validate_fallback_availability(result["availability"], f"{name}.availability")
    for block_name in ("coarse", "fine"):
        block = _obj(
            result[block_name],
            f"{name}.{block_name}",
            _schema._POCKET_VOLUME_BLOCK_FIELDS,
            _schema._POCKET_VOLUME_BLOCK_FIELDS,
        )
        _array(block["grid_shape"], f"{name}.{block_name}.grid_shape")
    grid = _obj(result["grid"], f"{name}.grid", _schema._POCKET_VOLUME_GRID_FIELDS, _schema._POCKET_VOLUME_GRID_FIELDS)
    for shape_name in ("shape", "coarse_shape", "fine_shape"):
        if grid[shape_name] is not None:
            _array(grid[shape_name], f"{name}.grid.{shape_name}")
    if grid["origin_xyz"] is not None:
        _array(grid["origin_xyz"], f"{name}.grid.origin_xyz")
    if result["sensitivity"] is not None:
        _obj(
            result["sensitivity"],
            f"{name}.sensitivity",
            _schema._POCKET_SENSITIVITY_FIELDS,
            _schema._POCKET_SENSITIVITY_FIELDS,
        )
    if result["settings"] is not None:
        _obj(
            result["settings"],
            f"{name}.settings",
            _schema._POCKET_VOLUME_SETTINGS_FIELDS,
            _schema._POCKET_VOLUME_SETTINGS_FIELDS,
        )
    _diagnostics(result["diagnostics"], f"{name}.diagnostics")
    if result["provenance"] is not None:
        _provenance(result["provenance"], f"{name}.provenance")


def _validate_match(value: object, name: str) -> None:
    match = _obj(value, name, _schema._POCKET_MATCH_FIELDS, _schema._POCKET_MATCH_FIELDS)
    for role in ("reference_candidate", "target_candidate"):
        if match[role] is not None:
            _candidate(match[role], f"{name}.{role}")
    _diagnostics(match["diagnostics"], f"{name}.diagnostics")


def _validate_comparison(value: object, name: str) -> None:
    comparison = _obj(value, name, _schema._POCKET_COMPARISON_FIELDS, _schema._POCKET_COMPARISON_FIELDS)
    _validate_match(comparison["match"], f"{name}.match")
    volume_value = comparison["volume"]
    if volume_value is not None:
        volume = _obj(
            volume_value,
            f"{name}.volume",
            _schema._POCKET_VOLUME_COMPARISON_FIELDS,
            _schema._POCKET_VOLUME_COMPARISON_FIELDS,
        )
        for sensitivity_name in ("reference_sensitivity", "target_sensitivity"):
            if volume[sensitivity_name] is not None:
                _obj(
                    volume[sensitivity_name],
                    f"{name}.volume.{sensitivity_name}",
                    _schema._POCKET_SENSITIVITY_FIELDS,
                    _schema._POCKET_SENSITIVITY_FIELDS,
                )
        _diagnostics(volume["diagnostics"], f"{name}.volume.diagnostics")
        for index, provenance in enumerate(_array(volume["provenance"], f"{name}.volume.provenance")):
            _provenance(provenance, f"{name}.volume.provenance[{index}]")
    for field in ("lining_residue_conserved", "lining_residue_gains", "lining_residue_losses"):
        for index, residue in enumerate(_array(comparison[field], f"{name}.{field}")):
            _schema._validate_fallback_residue(residue, f"{name}.{field}[{index}]")
    for index, mutation in enumerate(_array(comparison["associated_mutations"], f"{name}.associated_mutations")):
        mutation_name = f"{name}.associated_mutations[{index}]"
        item = _obj(mutation, mutation_name, _schema._MUTATION_FIELDS, _schema._MUTATION_FIELDS)
        _schema._validate_fallback_residue_or_none(item["reference"], f"{mutation_name}.reference")
        _schema._validate_fallback_residue_or_none(item["target"], f"{mutation_name}.target")
    for index, interaction in enumerate(_array(comparison["interaction_changes"], f"{name}.interaction_changes")):
        _schema._validate_fallback_interaction_difference(interaction, f"{name}.interaction_changes[{index}]")
    for index, displacement in enumerate(_array(comparison["local_displacements"], f"{name}.local_displacements")):
        item = _obj(
            displacement,
            f"{name}.local_displacements[{index}]",
            _schema._POCKET_COMPARISON_LOCAL_DISPLACEMENT_FIELDS,
            _schema._POCKET_COMPARISON_LOCAL_DISPLACEMENT_FIELDS,
        )
        _schema._validate_fallback_residue(item["residue"], f"{name}.local_displacements[{index}].residue")
    _diagnostics(comparison["qc_diagnostics"], f"{name}.qc_diagnostics")
    _diagnostics(comparison["diagnostics"], f"{name}.diagnostics")


def validate_pockets(value: object) -> None:
    pockets = _obj(value, "pockets", _schema._POCKET_FIELDS, _schema._POCKET_FIELDS)
    _schema._validate_fallback_availability(pockets["availability"], "pockets.availability")
    for index, raw in enumerate(_array(pockets["detections"], "pockets.detections")):
        detection = _obj(
            raw, f"pockets.detections[{index}]", _schema._POCKET_DETECTION_FIELDS, _schema._POCKET_DETECTION_FIELDS
        )
        _schema._validate_fallback_availability(detection["availability"], f"pockets.detections[{index}].availability")
        _candidate_array(detection["candidates"], f"pockets.detections[{index}].candidates")
        if detection["detection"] is not None:
            _validate_detection(detection["detection"], f"pockets.detections[{index}].detection")
        if detection["settings"] is not None:
            settings = _obj(
                detection["settings"],
                f"pockets.detections[{index}].settings",
                _schema._POCKET_SETTINGS_FIELDS,
                _schema._POCKET_SETTINGS_FIELDS,
            )
            _obj(
                settings["geometry"],
                f"pockets.detections[{index}].settings.geometry",
                _schema._POCKET_GEOMETRY_FIELDS,
                _schema._POCKET_GEOMETRY_FIELDS,
            )
        _diagnostics(detection["diagnostics"], f"pockets.detections[{index}].diagnostics")
        if detection["provenance"] is not None:
            _provenance(detection["provenance"], f"pockets.detections[{index}].provenance")
    for index, raw in enumerate(_array(pockets["volumes"], "pockets.volumes")):
        volume = _obj(raw, f"pockets.volumes[{index}]", _schema._POCKET_VOLUME_FIELDS, _schema._POCKET_VOLUME_FIELDS)
        if volume["result"] is not None:
            _validate_volume_result(volume["result"], f"pockets.volumes[{index}].result")
    matching = pockets["matching"]
    if matching is not None:
        result = _obj(matching, "pockets.matching", _schema._POCKET_MATCHING_FIELDS, _schema._POCKET_MATCHING_FIELDS)
        for index, match in enumerate(_array(result["matches"], "pockets.matching.matches")):
            _validate_match(match, f"pockets.matching.matches[{index}]")
        _diagnostics(result["diagnostics"], "pockets.matching.diagnostics")
        _obj(
            result["settings"],
            "pockets.matching.settings",
            _schema._POCKET_MATCHING_SETTINGS_FIELDS,
            _schema._POCKET_MATCHING_SETTINGS_FIELDS,
        )
        transform = _obj(
            result["transform"], "pockets.matching.transform", _schema._TRANSFORM_FIELDS, _schema._TRANSFORM_FIELDS
        )
        _array(transform["rotation"], "pockets.matching.transform.rotation")
        _array(transform["translation"], "pockets.matching.transform.translation")
    for index, raw in enumerate(_array(pockets["comparisons"], "pockets.comparisons")):
        _validate_comparison(raw, f"pockets.comparisons[{index}]")
    for index, raw in enumerate(_array(pockets["lining_residues"], "pockets.lining_residues")):
        lining = _obj(
            raw, f"pockets.lining_residues[{index}]", _schema._POCKET_LINING_FIELDS, _schema._POCKET_LINING_FIELDS
        )
        for residue_index, residue in enumerate(
            _array(lining["residues"], f"pockets.lining_residues[{index}].residues")
        ):
            _schema._validate_fallback_residue(residue, f"pockets.lining_residues[{index}].residues[{residue_index}]")
        _diagnostics(lining["diagnostics"], f"pockets.lining_residues[{index}].diagnostics")
    if pockets["concordance"] is not None:
        validate_concordance(pockets["concordance"], "pockets.concordance")
    _diagnostics(pockets["diagnostics"], "pockets.diagnostics")
    if pockets["provenance"] is not None:
        _provenance(pockets["provenance"], "pockets.provenance")


__all__ = ["validate_analysis", "validate_concordance", "validate_pockets"]
