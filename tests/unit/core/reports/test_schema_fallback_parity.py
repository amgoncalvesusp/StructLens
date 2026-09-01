from __future__ import annotations

import re
from copy import deepcopy
from typing import Any

import pytest

import structlens.core.reports.schema as schema_module
from structlens.core.reports.schema import AnalysisReportSchemaError, load_analysis_report_schema
from structlens.core.reports.schema_interpreter import SchemaInterpreterError, validate_instance


def _valid_instance(node: object, root: dict[str, Any], *, path: str = "#") -> Any:
    """Build a small valid instance for the bundled wire schema.

    This deliberately uses only JSON-Schema semantics, not the fallback
    implementation, so the differential cases exercise both validators.
    """

    if not isinstance(node, dict):
        raise AssertionError(f"schema node at {path} is not an object")
    if "$ref" in node:
        reference = node["$ref"]
        assert isinstance(reference, str) and reference.startswith("#/$defs/")
        definition = reference.removeprefix("#/$defs/")
        return _valid_instance(root["$defs"][definition], root, path=reference)
    if "const" in node:
        return deepcopy(node["const"])
    if "enum" in node:
        values = node["enum"]
        assert isinstance(values, list) and values
        return deepcopy(values[0])
    if "anyOf" in node and not any(key in node for key in ("type", "properties", "required")):
        alternatives = node["anyOf"]
        assert isinstance(alternatives, list) and alternatives
        return _valid_instance(alternatives[0], root, path=f"{path}.anyOf[0]")

    schema_type = node.get("type")
    if isinstance(schema_type, list):
        schema_type = next((item for item in schema_type if item != "null"), schema_type[0])
    if schema_type == "object":
        properties = node.get("properties", {})
        assert isinstance(properties, dict)
        required = node.get("required", [])
        assert isinstance(required, list)
        return {field: _valid_instance(properties[field], root, path=f"{path}.{field}") for field in required}
    if schema_type == "array":
        item_schema = node.get("items")
        minimum = int(node.get("minItems", 0))
        if item_schema is None:
            return []
        return [_valid_instance(item_schema, root, path=f"{path}[0]") for _ in range(minimum or 1)]
    if schema_type == "string":
        pattern = node.get("pattern")
        if isinstance(pattern, str) and "64" in pattern:
            return "a" * 64
        minimum = int(node.get("minLength", 0))
        return "x" * max(1, minimum)
    if schema_type == "integer":
        return 1
    if schema_type == "number":
        return 1.0
    if schema_type == "boolean":
        return False
    if schema_type == "null":
        return None
    if schema_type is None:
        return {}
    raise AssertionError(f"unsupported test fixture schema type {schema_type!r} at {path}")


def _valid_payload() -> dict[str, Any]:
    schema = load_analysis_report_schema()
    payload = _valid_instance(schema, schema)
    assert isinstance(payload, dict)
    payload["report_id"] = "a" * 64
    return payload


def _set_path(payload: dict[str, Any], path: tuple[str | int, ...], value: Any) -> None:
    current: Any = payload
    for part in path[:-1]:
        current = current[part]
    current[path[-1]] = value


def _without_jsonschema(monkeypatch: pytest.MonkeyPatch) -> None:
    real_import = schema_module.importlib.import_module

    def missing_validator(name: str) -> object:
        if name == "jsonschema":
            raise ImportError("jsonschema is unavailable")
        return real_import(name)

    monkeypatch.setattr(schema_module.importlib, "import_module", missing_validator)


@pytest.mark.parametrize(
    ("path", "value"),
    (
        (("schema_version",), True),
        (("analysis", "sequence_identity"), True),
        (("analysis", "mapped_residue_count"), "bad"),
        (("input_quality", "reference", "settings"), {}),
        (("diagnostics",), [{}]),
        (("provenance",), {}),
        (("msa", "sequences"), [{}]),
        (("msa", "aligned_rows"), [["one"]]),
        (("interactions",), [{}]),
        (("sites",), [{}]),
        (("displacement_vectors",), [{}]),
        (("evidence_cards",), [{}]),
        (("pockets", "detections"), [{}]),
    ),
)
def test_invalid_complete_report_is_rejected_with_bundled_and_fallback_schema(
    monkeypatch: pytest.MonkeyPatch, path: tuple[str | int, ...], value: Any
) -> None:
    payload = _valid_payload()
    _set_path(payload, path, value)

    with pytest.raises(AnalysisReportSchemaError):
        schema_module.validate_analysis_report_payload(payload)

    _without_jsonschema(monkeypatch)
    with pytest.raises(AnalysisReportSchemaError, match=re.escape("report")):
        schema_module.validate_analysis_report_payload(payload)


def test_complete_report_is_accepted_with_bundled_and_fallback_schema(monkeypatch: pytest.MonkeyPatch) -> None:
    payload = _valid_payload()
    assert schema_module.validate_analysis_report_payload(payload) == payload
    _without_jsonschema(monkeypatch)
    assert schema_module.validate_analysis_report_payload(payload) == payload


@pytest.mark.parametrize(
    ("schema", "instance", "message"),
    (
        ({"type": "number"}, True, "expected number"),
        ({"type": "integer"}, False, "expected integer"),
        ({"type": "string", "minLength": 2}, "x", "at least"),
        ({"type": "string", "pattern": "^[a-z]+$"}, "1", "pattern"),
        ({"type": "array", "minItems": 1}, [], "at least"),
        ({"type": "array", "maxItems": 1}, [1, 2], "at most"),
        ({"type": "object", "required": ["value"]}, {}, "missing"),
        ({"type": "object", "properties": {}, "additionalProperties": False}, {"extra": 1}, "unknown"),
        ({"type": "string", "enum": ["a"]}, "b", "enum"),
        ({"const": "expected"}, "other", "const"),
        ({"anyOf": [{"type": "string"}, {"type": "null"}]}, 1, "anyOf"),
        ({"type": "number"}, float("inf"), "expected number"),
    ),
)
def test_dependency_free_interpreter_reports_path_bearing_wire_errors(
    schema: dict[str, Any], instance: Any, message: str
) -> None:
    with pytest.raises(SchemaInterpreterError, match=message):
        validate_instance(instance, schema, path="report.analysis.value")


def test_dependency_free_interpreter_accepts_nullable_numeric_and_closed_valid_shapes() -> None:
    validate_instance(None, {"type": ["number", "null"]})
    validate_instance(1.0, {"type": "integer"})
    validate_instance(
        {"value": ["x"]},
        {
            "type": "object",
            "additionalProperties": False,
            "required": ["value"],
            "properties": {"value": {"type": "array", "items": {"type": "string"}}},
        },
    )


@pytest.mark.parametrize(
    "schema",
    (
        None,
        {"unsupported": True},
        {"$defs": []},
        {"$defs": {"x": []}},
        {"$defs": {"x": {"$defs": {}}}},
        {"$defs": {"x": {"title": "nested metadata"}}},
        {"$defs": {"x": {"type": "bogus"}}},
        {"type": []},
        {"type": ["string", 1]},
        {"additionalProperties": {}},
        {"required": [""]},
        {"required": ["x", "x"]},
        {"properties": []},
        {"properties": {"x": 1}},
        {"items": 1},
        {"anyOf": []},
        {"anyOf": [1]},
        {"minItems": -1},
        {"minItems": True},
        {"minItems": 2, "maxItems": 1},
        {"pattern": 1},
        {"pattern": "["},
        {"enum": []},
        {"const": float("nan")},
        {"$ref": "other"},
        {"$defs": {"x": {"type": "string"}}, "$ref": "#/$defs/x", "type": "string"},
    ),
)
def test_dependency_free_interpreter_rejects_unsupported_schema_constructs(schema: object) -> None:
    with pytest.raises(SchemaInterpreterError):
        validate_instance({}, schema)  # type: ignore[arg-type]


def test_dependency_free_interpreter_rejects_unknown_and_cyclic_local_references() -> None:
    unknown = {"type": "object", "properties": {"value": {"$ref": "#/$defs/missing"}}, "$defs": {}}
    with pytest.raises(SchemaInterpreterError, match="unknown local"):
        validate_instance({"value": "x"}, unknown)
    with pytest.raises(SchemaInterpreterError, match="unknown local"):
        validate_instance({}, unknown)

    cyclic = {
        "type": "object",
        "properties": {"value": {"$ref": "#/$defs/node"}},
        "$defs": {"node": {"$ref": "#/$defs/node"}},
    }
    with pytest.raises(SchemaInterpreterError, match="cyclic"):
        validate_instance({"value": "x"}, cyclic)
