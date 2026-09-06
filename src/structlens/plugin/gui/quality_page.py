"""Immutable presentation records for the coordinate-quality page.

The presenter copies native QC evidence only.  It does not score structures,
recalculate diagnostics, or turn missing measurements into zeroes.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import Enum
from types import MappingProxyType

from structlens.core.evidence import Availability, DiagnosticSeverity
from structlens.core.quality import StructureQualityReport


class QualityFilter(str, Enum):
    """Stable filters exposed by the QC page."""

    ALL = "all"
    ERRORS = "errors"
    WARNINGS = "warnings"
    INFO = "info"


@dataclass(frozen=True, slots=True)
class QualityRow:
    """One native diagnostic with all context needed by the detail drawer."""

    code: str
    severity: DiagnosticSeverity
    message: str
    source_id: str | None = None
    atom_id: str | None = None
    residue_id: str | None = None
    remediation: str | None = None


@dataclass(frozen=True, slots=True)
class QualitySummary:
    error_count: int = 0
    warning_count: int = 0
    info_count: int = 0

    @property
    def total_count(self) -> int:
        return self.error_count + self.warning_count + self.info_count


@dataclass(frozen=True, slots=True)
class QualityPresentation:
    """Display projection for one :class:`StructureQualityReport`."""

    availability: Availability
    rows: tuple[QualityRow, ...] = ()
    summary: QualitySummary = field(default_factory=QualitySummary)
    counts: Mapping[str, int] = field(default_factory=dict)
    settings: object | None = None
    provenance: object | None = None
    filter_value: QualityFilter = QualityFilter.ALL

    def __post_init__(self) -> None:
        object.__setattr__(self, "rows", tuple(self.rows))
        object.__setattr__(self, "counts", MappingProxyType(dict(self.counts)))


class QualityPresenter:
    """Pure conversion of a typed quality report into GUI values."""

    @staticmethod
    def present(
        report: StructureQualityReport,
        *,
        filter_value: QualityFilter = QualityFilter.ALL,
        filter: QualityFilter | None = None,
    ) -> QualityPresentation:
        if not isinstance(report, StructureQualityReport):
            raise TypeError("report must be a StructureQualityReport")
        if filter is not None:
            filter_value = filter
        if not isinstance(filter_value, QualityFilter):
            filter_value = QualityFilter(filter_value)
        diagnostics = tuple(report.diagnostics)
        rows = tuple(
            QualityRow(
                code=item.code,
                severity=item.severity,
                message=item.message,
                source_id=item.source_id,
                atom_id=item.atom_id,
                residue_id=item.residue_id,
                remediation=item.remediation,
            )
            for item in diagnostics
            if _matches(item.severity, filter_value)
        )
        summary = QualitySummary(
            error_count=sum(item.severity is DiagnosticSeverity.ERROR for item in diagnostics),
            warning_count=sum(item.severity is DiagnosticSeverity.WARNING for item in diagnostics),
            info_count=sum(item.severity is DiagnosticSeverity.INFO for item in diagnostics),
        )
        return QualityPresentation(
            availability=report.availability,
            rows=rows,
            summary=summary,
            counts=report.counts,
            settings=report.settings,
            provenance=report.provenance,
            filter_value=filter_value,
        )


def _matches(severity: DiagnosticSeverity, filter_value: QualityFilter) -> bool:
    return filter_value is QualityFilter.ALL or (
        filter_value is QualityFilter.ERRORS and severity is DiagnosticSeverity.ERROR
    ) or (
        filter_value is QualityFilter.WARNINGS and severity is DiagnosticSeverity.WARNING
    ) or (filter_value is QualityFilter.INFO and severity is DiagnosticSeverity.INFO)


__all__ = ["QualityFilter", "QualityPresenter", "QualityPresentation", "QualityRow", "QualitySummary"]
