"""Small dependency-free interpreter for the report schema wire contract.

The bundled report schema intentionally uses a finite Draft 2020-12 subset.
Keeping this interpreter separate makes the optional :mod:`jsonschema`
dependency an optimization rather than a correctness requirement.
"""

from __future__ import annotations

import math
import re
from collections.abc import Mapping
from typing import Any


class SchemaInterpreterError(ValueError):
    """Raised when a schema or instance is outside the supported subset."""


_ALLOWED_TYPES = frozenset({"array", "boolean", "integer", "null", "number", "object", "string"})
_ALLOWED_METADATA = frozenset({"$id", "$schema", "title"})
_ALLOWED_KEYWORDS = frozenset(
    {
        "$defs",
        "$ref",
        "additionalProperties",
        "anyOf",
        "const",
        "enum",
        "items",
        "maxItems",
        "minItems",
        "minLength",
        "pattern",
        "properties",
        "required",
        "type",
    }
)


def validate_instance(instance: object, schema: Mapping[str, Any], *, path: str = "report") -> None:
    """Validate *instance* against the supported bundled-schema subset.

    Schema validation is performed on every call so a corrupted or replaced
    bundled resource fails closed instead of silently widening the contract.
    """

    if not isinstance(schema, Mapping):
        raise SchemaInterpreterError("schema must be an object")
    _validate_schema_document(schema, path="#")
    _validate_node(instance, schema, schema, path, refs=())


def _validate_schema_document(schema: Mapping[str, Any], *, path: str) -> None:
    _validate_schema_node(schema, path=path, root=True)
    definitions = schema.get("$defs", {})
    if not isinstance(definitions, Mapping):
        raise SchemaInterpreterError(f"{path}.$defs must be an object")
    for name, definition in definitions.items():
        if not isinstance(name, str) or not name:
            raise SchemaInterpreterError(f"{path}.$defs contains an invalid definition name")
        if not isinstance(definition, Mapping):
            raise SchemaInterpreterError(f"{path}.$defs.{name} must be an object")
        _validate_schema_node(definition, path=f"{path}.$defs.{name}", root=False)
    _validate_local_references(schema, definitions, path=path)


def _validate_local_references(node: object, definitions: Mapping[str, Any], *, path: str) -> None:
    if isinstance(node, Mapping):
        reference = node.get("$ref")
        if isinstance(reference, str) and reference.startswith("#/$defs/"):
            name = reference.removeprefix("#/$defs/")
            if name not in definitions:
                raise SchemaInterpreterError(f"{path}.$ref points to unknown local schema {reference}")
        for name, child in node.items():
            if name == "$defs":
                if isinstance(child, Mapping):
                    for definition_name, definition in child.items():
                        _validate_local_references(
                            definition,
                            definitions,
                            path=f"{path}.$defs.{definition_name}",
                        )
            elif isinstance(child, (Mapping, list)):
                _validate_local_references(child, definitions, path=f"{path}.{name}")
    elif isinstance(node, list):
        for index, child in enumerate(node):
            _validate_local_references(child, definitions, path=f"{path}[{index}]")


def _validate_schema_node(node: Mapping[str, Any], *, path: str, root: bool) -> None:
    allowed = _ALLOWED_KEYWORDS | (_ALLOWED_METADATA if root else frozenset())
    unknown = sorted(set(node).difference(allowed))
    if unknown:
        raise SchemaInterpreterError(f"{path} uses unsupported schema keyword(s): {', '.join(unknown)}")

    if "$defs" in node:
        if not root:
            raise SchemaInterpreterError(f"{path}.$defs is only supported at the schema root")
        definitions = node["$defs"]
        if not isinstance(definitions, Mapping):
            raise SchemaInterpreterError(f"{path}.$defs must be an object")
    if "$ref" in node:
        reference = node["$ref"]
        if not isinstance(reference, str) or not reference.startswith("#/$defs/"):
            raise SchemaInterpreterError(f"{path}.$ref must be a local #/$defs reference")
        name = reference.removeprefix("#/$defs/")
        if not name or "/" in name or ".." in name:
            raise SchemaInterpreterError(f"{path}.$ref escapes the local definitions")
        if len(node) != 1:
            raise SchemaInterpreterError(f"{path} cannot combine $ref with other schema keywords")

    if "type" in node:
        schema_type = node["type"]
        if isinstance(schema_type, str):
            types = (schema_type,)
        elif isinstance(schema_type, list) and schema_type and all(isinstance(item, str) for item in schema_type):
            types = tuple(schema_type)
        else:
            raise SchemaInterpreterError(f"{path}.type must be a type name or non-empty type-name array")
        if any(item not in _ALLOWED_TYPES for item in types):
            raise SchemaInterpreterError(f"{path}.type contains an unsupported type")

    if "additionalProperties" in node and not isinstance(node["additionalProperties"], bool):
        raise SchemaInterpreterError(f"{path}.additionalProperties must be a boolean")

    if "required" in node:
        required = node["required"]
        if not isinstance(required, list) or not all(isinstance(item, str) and item for item in required):
            raise SchemaInterpreterError(f"{path}.required must be an array of non-empty strings")
        if len(set(required)) != len(required):
            raise SchemaInterpreterError(f"{path}.required contains duplicate fields")

    if "properties" in node:
        properties = node["properties"]
        if not isinstance(properties, Mapping):
            raise SchemaInterpreterError(f"{path}.properties must be an object")
        for name, child in properties.items():
            if not isinstance(name, str) or not name:
                raise SchemaInterpreterError(f"{path}.properties contains an invalid field name")
            if not isinstance(child, Mapping):
                raise SchemaInterpreterError(f"{path}.properties.{name} must be an object")
            _validate_schema_node(child, path=f"{path}.properties.{name}", root=False)

    for keyword in ("items",):
        if keyword in node:
            child = node[keyword]
            if not isinstance(child, Mapping):
                raise SchemaInterpreterError(f"{path}.{keyword} must be an object")
            _validate_schema_node(child, path=f"{path}.{keyword}", root=False)

    if "anyOf" in node:
        alternatives = node["anyOf"]
        if not isinstance(alternatives, list) or not alternatives:
            raise SchemaInterpreterError(f"{path}.anyOf must be a non-empty schema array")
        for index, alternative in enumerate(alternatives):
            if not isinstance(alternative, Mapping):
                raise SchemaInterpreterError(f"{path}.anyOf[{index}] must be an object")
            _validate_schema_node(alternative, path=f"{path}.anyOf[{index}]", root=False)

    for keyword in ("minItems", "maxItems", "minLength"):
        if keyword in node:
            bound = node[keyword]
            if not isinstance(bound, int) or isinstance(bound, bool) or bound < 0:
                raise SchemaInterpreterError(f"{path}.{keyword} must be a non-negative integer")
    if "minItems" in node and "maxItems" in node and node["minItems"] > node["maxItems"]:
        raise SchemaInterpreterError(f"{path}.minItems cannot exceed maxItems")

    if "pattern" in node:
        pattern = node["pattern"]
        if not isinstance(pattern, str):
            raise SchemaInterpreterError(f"{path}.pattern must be a string")
        try:
            re.compile(pattern)
        except re.error as exc:
            raise SchemaInterpreterError(f"{path}.pattern is not a valid regular expression") from exc

    if "enum" in node and (not isinstance(node["enum"], list) or not node["enum"]):
        raise SchemaInterpreterError(f"{path}.enum must be a non-empty array")
    if "const" in node:
        _validate_literal(node["const"], f"{path}.const")


def _validate_literal(value: object, path: str) -> None:
    if isinstance(value, float) and not math.isfinite(value):
        raise SchemaInterpreterError(f"{path} must be finite")
    if isinstance(value, Mapping):
        for key, child in value.items():
            _validate_literal(child, f"{path}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _validate_literal(child, f"{path}[{index}]")


def _validate_node(
    instance: object,
    node: Mapping[str, Any],
    root: Mapping[str, Any],
    path: str,
    *,
    refs: tuple[str, ...],
) -> None:
    if "$ref" in node:
        reference = node["$ref"]
        assert isinstance(reference, str)
        if reference in refs:
            raise SchemaInterpreterError(f"{path}: cyclic schema reference {reference}")
        definitions = root.get("$defs")
        if not isinstance(definitions, Mapping):
            raise SchemaInterpreterError(f"{path}: schema has no local definitions")
        name = reference.removeprefix("#/$defs/")
        if name not in definitions:
            raise SchemaInterpreterError(f"{path}: unknown local schema reference {reference}")
        definition = definitions[name]
        assert isinstance(definition, Mapping)
        _validate_node(instance, definition, root, path, refs=(*refs, reference))
        return

    if "const" in node and not _json_equal(instance, node["const"]):
        raise SchemaInterpreterError(f"{path}: expected const {node['const']!r}")
    if "enum" in node and not any(_json_equal(instance, candidate) for candidate in node["enum"]):
        raise SchemaInterpreterError(f"{path}: value is not a supported enum value")

    if "type" in node and not _matches_type(instance, node["type"]):
        expected = node["type"]
        expected_text = ", ".join(expected) if isinstance(expected, list) else str(expected)
        raise SchemaInterpreterError(f"{path}: expected {expected_text}, got {_json_type(instance)}")

    if isinstance(instance, Mapping):
        _validate_object(instance, node, root, path, refs)
    elif isinstance(instance, list):
        _validate_array(instance, node, root, path, refs)
    elif isinstance(instance, str):
        _validate_string(instance, node, path)
    elif isinstance(instance, (int, float)) and not isinstance(instance, bool):
        if isinstance(instance, float) and not math.isfinite(instance):
            raise SchemaInterpreterError(f"{path}: numbers must be finite")

    if "anyOf" in node:
        errors: list[str] = []
        for alternative in node["anyOf"]:
            try:
                _validate_node(instance, alternative, root, path, refs=refs)
            except SchemaInterpreterError as exc:
                errors.append(str(exc))
            else:
                break
        else:
            detail = errors[0] if errors else "no alternative matched"
            raise SchemaInterpreterError(f"{path}: value does not match anyOf ({detail})")


def _validate_object(
    instance: Mapping[str, Any],
    node: Mapping[str, Any],
    root: Mapping[str, Any],
    path: str,
    refs: tuple[str, ...],
) -> None:
    required = node.get("required", ())
    for name in required:
        if name not in instance:
            raise SchemaInterpreterError(f"{path}: missing required field {name!r}")
    properties = node.get("properties", {})
    if not isinstance(properties, Mapping):
        return
    if node.get("additionalProperties") is False:
        unknown = sorted(set(instance).difference(properties))
        if unknown:
            raise SchemaInterpreterError(f"{path}: unknown field(s): {', '.join(unknown)}")
    for name, child in properties.items():
        if name in instance:
            assert isinstance(child, Mapping)
            _validate_node(instance[name], child, root, f"{path}.{name}", refs=refs)


def _validate_array(
    instance: list[Any],
    node: Mapping[str, Any],
    root: Mapping[str, Any],
    path: str,
    refs: tuple[str, ...],
) -> None:
    minimum = node.get("minItems")
    maximum = node.get("maxItems")
    if isinstance(minimum, int) and len(instance) < minimum:
        raise SchemaInterpreterError(f"{path}: must contain at least {minimum} item(s)")
    if isinstance(maximum, int) and len(instance) > maximum:
        raise SchemaInterpreterError(f"{path}: must contain at most {maximum} item(s)")
    child = node.get("items")
    if isinstance(child, Mapping):
        for index, value in enumerate(instance):
            _validate_node(value, child, root, f"{path}[{index}]", refs=refs)


def _validate_string(instance: str, node: Mapping[str, Any], path: str) -> None:
    minimum = node.get("minLength")
    if isinstance(minimum, int) and len(instance) < minimum:
        if minimum == 1:
            raise SchemaInterpreterError(f"{path}: must be non-empty (at least 1 character)")
        raise SchemaInterpreterError(f"{path}: must contain at least {minimum} character(s)")
    pattern = node.get("pattern")
    if isinstance(pattern, str) and re.search(pattern, instance) is None:
        pattern_name = "SHA-256 pattern" if "64" in pattern else "required pattern"
        raise SchemaInterpreterError(f"{path}: does not match required {pattern_name}")


def _matches_type(value: object, schema_type: object) -> bool:
    types = schema_type if isinstance(schema_type, list) else [schema_type]
    return any(
        (item == "null" and value is None)
        or (item == "boolean" and type(value) is bool)
        or (item == "object" and isinstance(value, Mapping))
        or (item == "array" and isinstance(value, list))
        or (item == "string" and isinstance(value, str))
        or (item == "number" and _is_number(value))
        or (item == "integer" and _is_integer(value))
        for item in types
    )


def _is_number(value: object) -> bool:
    return (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and (not isinstance(value, float) or math.isfinite(value))
    )


def _is_integer(value: object) -> bool:
    if not _is_number(value):
        return False
    if isinstance(value, int):
        return True
    return isinstance(value, float) and value.is_integer()


def _json_type(value: object) -> str:
    if value is None:
        return "null"
    if type(value) is bool:
        return "boolean"
    if isinstance(value, int):
        return "integer"
    if isinstance(value, float):
        return "number"
    if isinstance(value, str):
        return "string"
    if isinstance(value, list):
        return "array"
    if isinstance(value, Mapping):
        return "object"
    return type(value).__name__


def _json_equal(left: object, right: object) -> bool:
    if type(left) is not type(right):
        if (
            isinstance(left, (int, float))
            and not isinstance(left, bool)
            and isinstance(right, (int, float))
            and not isinstance(right, bool)
        ):
            return float(left) == float(right)
        return False
    return left == right


__all__ = ["SchemaInterpreterError", "validate_instance"]
