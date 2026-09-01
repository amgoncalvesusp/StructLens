from __future__ import annotations

import importlib
import json
from copy import deepcopy
from typing import Any

import pytest

import structlens.core.reports.schema as schema_module
from structlens.application.report_serialization import report_to_dict
from structlens.core.reports.schema import (
    AnalysisReportSchemaError,
    load_analysis_report_schema,
    validate_analysis_report_payload,
    validate_report_schema,
)

fixtures = importlib.import_module("tests.unit.application.test_task11_round1")


def _valid_payload() -> dict[str, Any]:
    reference = fixtures._snapshot("reference.pdb", b"reference")
    target = fixtures._snapshot("target.pdb", b"target")
    return report_to_dict(fixtures._report(reference, target))


def test_analysis_report_schema_is_packaged_and_versioned() -> None:
    schema = load_analysis_report_schema()
    assert schema["$id"].endswith("analysis-report-v1.json")
    assert schema["properties"]["schema_version"]["const"] == "4.0"


def test_schema_rejects_future_schema_and_non_finite_numbers() -> None:
    with pytest.raises(AnalysisReportSchemaError, match="schema"):
        validate_analysis_report_payload({"schema_version": "5.0"})

    with pytest.raises(AnalysisReportSchemaError, match="finite|JSON"):
        validate_analysis_report_payload({"schema_version": "4.0", "value": float("nan")})


def test_schema_validation_does_not_accept_non_object_payload() -> None:
    with pytest.raises(AnalysisReportSchemaError, match="object"):
        validate_analysis_report_payload(json.dumps([]))

    with pytest.raises(AnalysisReportSchemaError, match="Invalid"):
        validate_analysis_report_payload(b"\xff")


def test_schema_fallback_is_fail_closed(monkeypatch: pytest.MonkeyPatch) -> None:
    """The optional validator cannot make malformed wire data acceptable."""

    def missing_validator(_name: str) -> object:
        raise ImportError("jsonschema is not installed")

    monkeypatch.setattr("structlens.core.reports.schema.importlib.import_module", missing_validator)
    with pytest.raises(AnalysisReportSchemaError, match="unknown fields"):
        validate_analysis_report_payload({"schema_version": "4.0", "unexpected": True})


def test_schema_rejects_missing_and_invalid_digest_and_alias_accepts_valid_payload() -> None:
    with pytest.raises(AnalysisReportSchemaError, match="missing required fields"):
        validate_analysis_report_payload({"schema_version": "4.0"})

    invalid_digest = {field: None for field in schema_module._REPORT_FIELDS}
    invalid_digest["schema_version"] = "4.0"
    invalid_digest["report_id"] = "invalid"
    with pytest.raises(AnalysisReportSchemaError, match="SHA-256"):
        validate_analysis_report_payload(invalid_digest)

    payload = _valid_payload()
    assert validate_report_schema(payload) == payload


def test_copy_json_enforces_mapping_size_and_field_name_bounds() -> None:
    with pytest.raises(AnalysisReportSchemaError, match="too many fields"):
        schema_module._copy_json({str(index): index for index in range(100_001)}, "payload")
    with pytest.raises(AnalysisReportSchemaError, match="invalid field name"):
        schema_module._copy_json({1: "value"}, "payload")


def test_fallback_shape_covers_optional_values_and_type_boundaries() -> None:
    payload = _valid_payload()
    payload["input_quality"]["reference"]["provenance"] = {}
    residue = {field: None for field in schema_module._RESIDUE_FIELDS}
    site = {field: None for field in schema_module._SITE_DEFINITION_FIELDS}
    site["reference_residues"] = []
    site["center_residue"] = residue
    payload["site_definitions"] = [site]
    schema_module._validate_fallback_shape(payload)

    invalid_sites = _valid_payload()
    invalid_sites["sites"] = 1
    with pytest.raises(AnalysisReportSchemaError, match="object|array, object, or null"):
        schema_module._validate_fallback_shape(invalid_sites)

    with pytest.raises(AnalysisReportSchemaError, match="missing required fields"):
        schema_module._validate_fallback_object({}, "object", {"required"}, required={"required"})
    with pytest.raises(AnalysisReportSchemaError, match="availability"):
        schema_module._validate_fallback_availability("invalid", "availability")
    with pytest.raises(AnalysisReportSchemaError, match="array"):
        schema_module._validate_fallback_diagnostics({}, "diagnostics")


def test_fallback_nested_collections_cover_msa_interactions_sites_and_candidates() -> None:
    msa = {field: [] for field in schema_module._MSA_FIELDS}
    msa["aligned_rows"] = [["one"]]
    with pytest.raises(AnalysisReportSchemaError, match="two values"):
        schema_module._validate_fallback_msa(msa)

    residue = {field: None for field in schema_module._RESIDUE_FIELDS}
    record = {field: None for field in schema_module._INTERACTION_RECORD_FIELDS}
    record["residue_a"] = residue
    key = {field: None for field in schema_module._INTERACTION_KEY_FIELDS}
    key["interaction_type"] = "contact"
    key["reference_position_a"] = "10"
    difference = {
        "key": key,
        "change": None,
        "reference_record": record,
        "target_record": record,
    }
    schema_module._validate_fallback_interaction_difference(difference, "difference")
    schema_module._validate_fallback_interaction_item(difference, "difference")
    schema_module._validate_fallback_interaction_item(record, "record")
    schema_module._validate_fallback_interactions([record])
    schema_module._validate_fallback_interactions(
        {"differences": [difference], "reference_interactions": [record], "target_interactions": [record]}
    )

    metric = {field: None for field in schema_module._SITE_METRIC_FIELDS}
    schema_module._validate_fallback_sites([metric])
    schema_module._validate_fallback_sites({"metrics": [metric]})

    sphere = {field: None for field in schema_module._POCKET_SPHERE_FIELDS}
    candidate = {field: None for field in schema_module._POCKET_CANDIDATE_FIELDS}
    candidate["alpha_spheres"] = [sphere]
    candidate["lineage"] = {field: None for field in schema_module._POCKET_LINEAGE_FIELDS}
    schema_module._validate_fallback_candidate(candidate, "candidate")


def test_nested_schema_validates_optional_analysis_concordance_and_pocket_paths() -> None:
    nested = importlib.import_module("structlens.core.reports.schema_nested")
    residue = {field: None for field in schema_module._RESIDUE_FIELDS}
    correspondence = {field: None for field in schema_module._CORRESPONDENCE_FIELDS}
    correspondence["reference"] = residue
    correspondence["target"] = residue
    mutation = {field: None for field in schema_module._MUTATION_FIELDS}
    mutation["reference"] = residue
    mutation["target"] = residue
    analysis = {field: None for field in schema_module._ANALYSIS_FIELDS}
    analysis["correspondences"] = [correspondence]
    analysis["mutations"] = [mutation]
    analysis["transform"] = {"rotation": [[1, 2, 3]], "translation": [1, 2, 3]}
    analysis["method_provenance"] = {}
    nested.validate_analysis(analysis)

    channel = {field: None for field in schema_module._POCKET_CHANNEL_FIELDS}
    channel["availability"] = "available"
    channel["diagnostics"] = []
    concordance = {
        name: deepcopy(channel)
        for name in (
            "geometry",
            "ligand_support",
            "parameter_persistence",
            "volume_sensitivity",
            "match_ambiguity",
            "qc",
            "interactions",
        )
    }
    nested.validate_concordance(concordance, "concordance")
    with pytest.raises(AnalysisReportSchemaError, match="seven channels"):
        nested.validate_concordance({}, "concordance")
    malformed_concordance = deepcopy(concordance)
    malformed_concordance["unexpected"] = channel
    with pytest.raises(AnalysisReportSchemaError, match="unknown"):
        nested.validate_concordance(malformed_concordance, "concordance")

    reference = fixtures._snapshot("reference.pdb", b"reference")
    target = fixtures._snapshot("target.pdb", b"target")
    report_payload = report_to_dict(
        fixtures._report(reference, target, pockets=fixtures._pocket_report(reference, target))
    )
    pockets = report_payload["pockets"]
    pockets["volumes"][0]["result"]["grid"]["origin_xyz"] = [0, 0, 0]
    pockets["comparisons"][0]["volume"]["provenance"] = [{}]
    difference = {
        "key": {
            "interaction_type": "contact",
            "reference_position_a": "10",
            "reference_position_b": None,
            "external_partner_id": None,
        },
        "change": None,
        "reference_record": None,
        "target_record": None,
    }
    pockets["comparisons"][0]["interaction_changes"] = [difference]
    nested.validate_pockets(pockets)
