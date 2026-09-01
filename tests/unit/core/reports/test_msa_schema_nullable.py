from __future__ import annotations

import importlib
from copy import deepcopy
from typing import Any

import pytest

import structlens.core.reports.schema as schema_module
from structlens.application.report_serialization import report_to_dict
from structlens.core.reports.schema import AnalysisReportSchemaError, validate_analysis_report_payload

fixtures = importlib.import_module("tests.unit.application.test_task11_round1")


def _ambiguous_msa_payload() -> dict[str, Any]:
    reference = fixtures._snapshot("reference.pdb", b"reference")
    target = fixtures._snapshot("target.pdb", b"target")
    payload = deepcopy(report_to_dict(fixtures._report(reference, target)))
    payload["msa"] = {
        "sequences": [],
        "aligned_rows": [["reference", "X"], ["target", "X"]],
        "columns": [
            {
                "index": 0,
                "reference_label": "A:1",
                "reference_residue": None,
                "cells": [],
                "non_gap_count": 2,
                "gap_fraction": 0.0,
                "ambiguous_fraction": 1.0,
                "conservation_score": None,
                "entropy_bits": None,
            }
        ],
        "reference_structure_id": "reference",
        "algorithm": "fallback",
        "provenance": ["test"],
    }
    return payload


@pytest.mark.parametrize("without_jsonschema", [False, True])
def test_msa_ambiguous_metrics_are_nullable_in_bundled_and_fallback_schema(
    without_jsonschema: bool, monkeypatch: pytest.MonkeyPatch
) -> None:
    if without_jsonschema:
        real_import = importlib.import_module

        def missing_jsonschema(name: str) -> object:
            if name == "jsonschema":
                raise ImportError("optional dependency absent")
            return real_import(name)

        monkeypatch.setattr(schema_module.importlib, "import_module", missing_jsonschema)

    payload = _ambiguous_msa_payload()
    assert validate_analysis_report_payload(payload) == payload

    for field in ("conservation_score", "entropy_bits"):
        missing = deepcopy(payload)
        del missing["msa"]["columns"][0][field]
        with pytest.raises(AnalysisReportSchemaError, match="schema|missing|required"):
            validate_analysis_report_payload(missing)

        for invalid_value in ("unknown", True):
            wrong_type = deepcopy(payload)
            wrong_type["msa"]["columns"][0][field] = invalid_value
            with pytest.raises(AnalysisReportSchemaError, match="schema|number|float|type"):
                validate_analysis_report_payload(wrong_type)

        unexpected = deepcopy(payload)
        unexpected["msa"]["columns"][0]["unexpected"] = True
        with pytest.raises(AnalysisReportSchemaError, match="unknown|additional"):
            validate_analysis_report_payload(unexpected)

    non_finite = deepcopy(payload)
    non_finite["msa"]["columns"][0]["entropy_bits"] = float("nan")
    with pytest.raises(AnalysisReportSchemaError, match="finite|JSON"):
        validate_analysis_report_payload(non_finite)
