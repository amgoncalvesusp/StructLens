"""Typed availability and diagnostic contracts used by scientific services.

These values deliberately form a small closed vocabulary.  A consumer can map
the stable codes to localized user-facing text without making scientific code
depend on a particular GUI.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class Availability(str, Enum):
    """State of a scientific result or optional analysis section."""

    AVAILABLE = "available"
    NOT_APPLICABLE = "not_applicable"
    NOT_DETECTED = "not_detected"
    INVALID_INPUT = "invalid_input"
    DEPENDENCY_UNAVAILABLE = "dependency_unavailable"
    NUMERICAL_FAILURE = "numerical_failure"


class DiagnosticSeverity(str, Enum):
    """Severity of an actionable, stable diagnostic code."""

    INFO = "info"
    WARNING = "warning"
    ERROR = "error"


# Short alias for callers that prefer the generic term.
Severity = DiagnosticSeverity


@dataclass(frozen=True, slots=True)
class Diagnostic:
    """A deterministic diagnostic attached to a source or structural object."""

    code: str
    severity: DiagnosticSeverity
    message: str
    source_id: str | None = None
    atom_id: str | None = None
    residue_id: str | None = None
    remediation: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.code, str) or not self.code.strip() or any(char.isspace() for char in self.code):
            raise ValueError("diagnostic code must be a non-empty stable token without whitespace")
        if not isinstance(self.message, str) or not self.message.strip():
            raise ValueError("diagnostic message must be non-empty")
        severity = self.severity
        if not isinstance(severity, DiagnosticSeverity):
            try:
                severity = DiagnosticSeverity(severity)
            except (TypeError, ValueError) as exc:
                raise ValueError(f"unknown diagnostic severity: {self.severity!r}") from exc
            object.__setattr__(self, "severity", severity)
        for field_name in ("source_id", "atom_id", "residue_id", "remediation"):
            value = getattr(self, field_name)
            if value is not None and (not isinstance(value, str) or not value.strip()):
                raise ValueError(f"{field_name} must be a non-empty string when provided")

    def to_json(self) -> dict[str, str | None]:
        """Return a JSON-compatible copy without exposing internal state."""

        return {
            "code": self.code,
            "severity": self.severity.value,
            "message": self.message,
            "source_id": self.source_id,
            "atom_id": self.atom_id,
            "residue_id": self.residue_id,
            "remediation": self.remediation,
        }


__all__ = ["Availability", "Diagnostic", "DiagnosticSeverity", "Severity"]
