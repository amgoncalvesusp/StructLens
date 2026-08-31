"""Typed failure and diagnostic helpers for the report workflow."""

from __future__ import annotations

from collections.abc import Iterable

from structlens.application.dto import AnalysisReportRequest
from structlens.application.report_input import canonical_site_definition
from structlens.core.evidence import Availability, Diagnostic, DiagnosticSeverity
from structlens.core.provenance import MethodProvenance
from structlens.core.reports import AnalysisReport, InputQualityBundle, SectionAvailability


def invalid_input_report(
    request: AnalysisReportRequest,
    input_quality: InputQualityBundle,
    provenance: MethodProvenance,
    diagnostics: tuple[Diagnostic, ...],
) -> AnalysisReport:
    blocked_state = (
        Availability.NUMERICAL_FAILURE
        if Availability.NUMERICAL_FAILURE
        in {input_quality.reference.availability, input_quality.target.availability}
        else Availability.INVALID_INPUT
    )
    state = SectionAvailability(
        input_quality=blocked_state,
        analysis=blocked_state,
        msa=blocked_state,
        interactions=blocked_state,
        sites=blocked_state,
        distance_map=blocked_state,
        displacement_vectors=blocked_state,
        evidence_cards=blocked_state,
    )
    return AnalysisReport(
        request.reference_selection,
        request.target_selection,
        input_quality,
        site_definitions=tuple(canonical_site_definition(item) for item in request.site_definitions),
        availability=state,
        diagnostics=diagnostics,
        provenance=provenance,
    )


def analysis_failure_report(
    request: AnalysisReportRequest,
    input_quality: InputQualityBundle,
    provenance: MethodProvenance,
    diagnostics: tuple[Diagnostic, ...],
) -> AnalysisReport:
    state = SectionAvailability(
        input_quality=Availability.AVAILABLE,
        analysis=Availability.NUMERICAL_FAILURE,
        msa=Availability.NUMERICAL_FAILURE,
        interactions=Availability.NUMERICAL_FAILURE,
        sites=Availability.NUMERICAL_FAILURE,
        distance_map=Availability.NUMERICAL_FAILURE,
        displacement_vectors=Availability.NUMERICAL_FAILURE,
        evidence_cards=Availability.NUMERICAL_FAILURE,
    )
    return AnalysisReport(
        request.reference_selection,
        request.target_selection,
        input_quality,
        site_definitions=tuple(canonical_site_definition(item) for item in request.site_definitions),
        availability=state,
        diagnostics=diagnostics,
        provenance=provenance,
    )


def section_failure(name: str, error: Exception) -> Diagnostic:
    return Diagnostic(
        code=f"report.{name}.failed",
        severity=DiagnosticSeverity.ERROR,
        message=f"The {name.replace('_', ' ')} section could not be calculated ({type(error).__name__}).",
        remediation="Review the input QC and the section-specific method settings.",
    )


def merge_diagnostics(*groups: Iterable[Diagnostic]) -> tuple[Diagnostic, ...]:
    """Deduplicate diagnostics and sort them for stable report identity."""

    unique: dict[tuple[object, ...], Diagnostic] = {}
    for item in (diagnostic for group in groups for diagnostic in group):
        key = (
            item.code,
            item.severity,
            item.message,
            item.source_id,
            item.atom_id,
            item.residue_id,
            item.remediation,
        )
        unique.setdefault(key, item)
    return tuple(
        sorted(
            unique.values(),
            key=lambda item: (item.code, item.source_id or "", item.residue_id or "", item.atom_id or ""),
        )
    )


__all__ = [
    "analysis_failure_report",
    "invalid_input_report",
    "merge_diagnostics",
    "section_failure",
]
