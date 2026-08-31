"""Tests for typed availability and diagnostic contracts."""

from dataclasses import FrozenInstanceError

import pytest

from structlens.core.evidence import Availability, Diagnostic, DiagnosticSeverity


def test_availability_is_a_closed_immutable_vocabulary() -> None:
    assert Availability.AVAILABLE.value == "available"
    assert Availability.NOT_DETECTED.value == "not_detected"
    assert Availability.NUMERICAL_FAILURE.value == "numerical_failure"
    with pytest.raises(ValueError, match="unknown"):
        Availability("unknown")


def test_diagnostic_normalizes_severity_and_is_json_ready() -> None:
    diagnostic = Diagnostic(
        code="coordinate.non_finite",
        severity="error",
        message="Coordinate is not finite.",
        source_id="target",
        atom_id="A:10:CA",
        residue_id="A:10",
        remediation="Review the source coordinates.",
    )

    assert diagnostic.severity is DiagnosticSeverity.ERROR
    assert diagnostic.to_json() == {
        "code": "coordinate.non_finite",
        "severity": "error",
        "message": "Coordinate is not finite.",
        "source_id": "target",
        "atom_id": "A:10:CA",
        "residue_id": "A:10",
        "remediation": "Review the source coordinates.",
    }
    with pytest.raises(FrozenInstanceError):
        diagnostic.code = "changed"  # type: ignore[misc]


def test_diagnostic_rejects_unknown_severity_and_empty_code() -> None:
    with pytest.raises(ValueError, match="severity"):
        Diagnostic("qc.warning", "notice", "message")
    with pytest.raises(ValueError, match="code"):
        Diagnostic("", DiagnosticSeverity.WARNING, "message")
