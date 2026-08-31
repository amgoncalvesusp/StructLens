"""Tests for deterministic, immutable method provenance."""

from dataclasses import FrozenInstanceError
from datetime import UTC, datetime

import pytest

import structlens.core.provenance as provenance_module
from structlens.core.models.results import AnalysisResult
from structlens.core.provenance import AuditEvent, MethodProvenance

HASH_A = "a" * 64


def _provenance() -> MethodProvenance:
    return MethodProvenance(
        method_id="structlens.pocket.detect",
        method_version="0.4.0",
        parameters={"grid_spacing": 0.5, "nested": {"enabled": True}, "steps": [1, 2]},
        units={"grid_spacing": "angstrom"},
        backend_versions={"scipy": "1.17.1"},
        input_hashes={"target": HASH_A},
        analyzed_representation="as_file",
    )


def test_method_provenance_deep_freezes_and_converts_to_json() -> None:
    provenance = _provenance()

    assert provenance.parameters["nested"]["enabled"] is True  # type: ignore[index]
    assert provenance.parameters["steps"] == (1, 2)
    assert provenance.to_json()["parameters"] == {
        "grid_spacing": 0.5,
        "nested": {"enabled": True},
        "steps": [1, 2],
    }
    with pytest.raises(TypeError):
        provenance.parameters["new"] = 1  # type: ignore[index]
    with pytest.raises(TypeError):
        provenance.parameters["nested"]["enabled"] = False  # type: ignore[index]
    with pytest.raises(FrozenInstanceError):
        provenance.method_id = "changed"  # type: ignore[misc]


def test_method_provenance_bounds_units_before_copy_and_rejects_oversized_strings() -> None:
    class OversizedUnits(dict[str, str]):
        def __len__(self) -> int:
            return provenance_module.MAX_UNIT_MAP_ITEMS + 1

        def items(self) -> object:
            raise AssertionError("oversized units must be rejected before iteration")

    with pytest.raises(ValueError, match="units.*item limit"):
        MethodProvenance("structlens.test", "1", units=OversizedUnits(), input_hashes={"target": HASH_A})

    huge = "x" * (provenance_module.MAX_UNIT_STRING_LENGTH + 1)
    for units in ({huge: "count"}, {"measure": huge}):
        with pytest.raises(ValueError, match="units.*string length"):
            MethodProvenance("structlens.test", "1", units=units, input_hashes={"target": HASH_A})


@pytest.mark.parametrize("value", [float("nan"), float("inf"), float("-inf")])
def test_method_provenance_rejects_non_finite_parameters(value: float) -> None:
    with pytest.raises(ValueError, match="finite"):
        MethodProvenance(
            "structlens.test",
            "1",
            {"grid_spacing": value},
            {"grid_spacing": "angstrom"},
            input_hashes={"target": HASH_A},
        )


def test_method_provenance_requires_units_for_physical_parameters() -> None:
    with pytest.raises(ValueError, match="unit"):
        MethodProvenance(
            "structlens.test",
            "1",
            {"grid_spacing": 0.5},
            {},
            input_hashes={"target": HASH_A},
        )


def test_method_provenance_rejects_inconsistent_source_hashes() -> None:
    with pytest.raises(ValueError, match="hash"):
        MethodProvenance(
            "structlens.test",
            "1",
            {},
            {},
            input_hashes={"target": "not-a-sha256"},
        )


def test_method_provenance_requires_at_least_one_source_hash() -> None:
    with pytest.raises(ValueError, match="at least one.*hash"):
        MethodProvenance("structlens.test", "1", {}, {})


def test_method_provenance_validates_nested_physical_paths_and_units() -> None:
    provenance = MethodProvenance(
        "structlens.pocket.measure",
        "1",
        {"pocket": {"probe_radius": 1.4}, "atom_count": 12},
        {"pocket.probe_radius": "angstrom"},
        input_hashes={"target": HASH_A},
    )

    assert provenance.units["pocket.probe_radius"] == "angstrom"

    with pytest.raises(ValueError, match="pocket.probe_radius.*unit"):
        MethodProvenance(
            "structlens.pocket.measure",
            "1",
            {"pocket": {"probe_radius": 1.4}, "atom_count": 12},
            {},
            input_hashes={"target": HASH_A},
        )


@pytest.mark.parametrize(
    ("parameter_name", "unit"),
    [("probe_radius", "angstrom^3"), ("pocket_volume", "angstrom")],
)
def test_method_provenance_rejects_dimensionally_incompatible_units(parameter_name: str, unit: str) -> None:
    with pytest.raises(ValueError, match="incompatible unit"):
        MethodProvenance(
            "structlens.pocket.measure",
            "1",
            {parameter_name: 1.4},
            {parameter_name: unit},
            input_hashes={"target": HASH_A},
        )


def test_method_provenance_rejects_unsupported_units() -> None:
    with pytest.raises(ValueError, match="unsupported unit"):
        MethodProvenance(
            "structlens.pocket.measure",
            "1",
            {"probe_radius": 1.4},
            {"probe_radius": "parsec"},
            input_hashes={"target": HASH_A},
        )


def test_method_provenance_rejects_physical_unit_on_dimensionless_count() -> None:
    with pytest.raises(ValueError, match="incompatible unit"):
        MethodProvenance(
            "structlens.pocket.measure",
            "1",
            {"atom_count": 12},
            {"atom_count": "angstrom"},
            input_hashes={"target": HASH_A},
        )


def test_method_provenance_from_json_recomputes_artifact_id() -> None:
    original = _provenance()
    restored = MethodProvenance.from_json(original.to_json())

    assert restored == original
    assert restored.canonical_bytes() == original.canonical_bytes()
    tampered = original.to_json()
    tampered["artifact_id"] = "0" * 64
    with pytest.raises(ValueError, match="artifact_id"):
        MethodProvenance.from_json(tampered)


def test_artifact_id_and_canonical_bytes_exclude_audit_time() -> None:
    first = _provenance()
    second = _provenance()
    first_event = AuditEvent(
        event_type="analysis.started",
        occurred_at=datetime(2026, 8, 30, 12, 0, tzinfo=UTC),
        user="analyst",
        workstation="one",
    )
    second_event = AuditEvent(
        event_type="analysis.started",
        occurred_at=datetime(2026, 8, 30, 12, 1, tzinfo=UTC),
        user="analyst",
        workstation="one",
    )

    assert first.artifact_id == second.artifact_id
    assert first.canonical_bytes() == second.canonical_bytes()
    assert first.with_audit_event(first_event)["scientific"] == second.with_audit_event(second_event)["scientific"]
    assert (
        first.with_audit_event(first_event)["audit_event"]["occurred_at"]
        != second.with_audit_event(second_event)["audit_event"]["occurred_at"]
    )


def test_analysis_result_requires_typed_record_in_explicit_field() -> None:
    typed = _provenance()
    with pytest.raises(TypeError, match="method_provenance"):
        AnalysisResult(
            reference_id="ref",
            target_id="target",
            correspondences=(),
            mutations=(),
            sequence_identity=1.0,
            sequence_coverage=1.0,
            alignment_decision="sequence",
            provenance=typed,  # type: ignore[arg-type]
        )

    result = AnalysisResult(
        reference_id="ref",
        target_id="target",
        correspondences=(),
        mutations=(),
        sequence_identity=1.0,
        sequence_coverage=1.0,
        alignment_decision="sequence",
        provenance={"backend": "legacy"},
        method_provenance=typed,
    )

    assert result.typed_provenance is typed
    assert result.provenance == {"backend": "legacy"}
