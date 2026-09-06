"""RED contracts for the report-backed Quality GUI page."""

from __future__ import annotations

from dataclasses import is_dataclass

import pytest
from structlens.plugin.gui.quality_page import QualityFilter, QualityPresenter

from structlens.core.evidence import Availability, Diagnostic, DiagnosticSeverity
from structlens.core.quality import CoordinateQCSettings, StructureQualityReport


def _report() -> StructureQualityReport:
    return StructureQualityReport(
        availability=Availability.AVAILABLE,
        diagnostics=(
            Diagnostic(
                "quality.coordinate.nonfinite",
                DiagnosticSeverity.ERROR,
                "Atom coordinate is not finite.",
                source_id="reference",
                residue_id="A:42",
                remediation="Review the deposited coordinate record.",
            ),
            Diagnostic(
                "quality.coordinate.altloc",
                DiagnosticSeverity.WARNING,
                "Alternate conformer has lower occupancy.",
                source_id="reference",
                residue_id="A:43",
                remediation="Inspect the alternate conformer policy.",
            ),
            Diagnostic(
                "quality.coordinate.info",
                DiagnosticSeverity.INFO,
                "Coordinate record was retained.",
                source_id="reference",
            ),
        ),
        counts={"atoms": 120, "residues": 30},
        settings=CoordinateQCSettings(),
    )


def test_quality_presenter_is_immutable_and_retains_typed_severity_and_details() -> None:
    presentation = QualityPresenter.present(_report())

    assert is_dataclass(presentation)
    assert getattr(type(presentation), "__slots__", None)
    assert presentation.availability is Availability.AVAILABLE
    assert [row.code for row in presentation.rows] == [
        "quality.coordinate.nonfinite",
        "quality.coordinate.altloc",
        "quality.coordinate.info",
    ]
    assert presentation.rows[0].severity is DiagnosticSeverity.ERROR
    assert presentation.rows[0].residue_id == "A:42"
    assert presentation.rows[0].remediation == "Review the deposited coordinate record."
    assert presentation.summary.error_count == 1
    assert presentation.summary.warning_count == 1
    with pytest.raises((AttributeError, TypeError)):
        presentation.rows = ()


@pytest.mark.parametrize(
    ("filter_value", "expected"),
    (
        (QualityFilter.ALL, {"quality.coordinate.nonfinite", "quality.coordinate.altloc", "quality.coordinate.info"}),
        (QualityFilter.ERRORS, {"quality.coordinate.nonfinite"}),
        (QualityFilter.WARNINGS, {"quality.coordinate.altloc"}),
        (QualityFilter.INFO, {"quality.coordinate.info"}),
    ),
)
def test_quality_presenter_filters_without_reordering_or_losing_native_diagnostics(
    filter_value: QualityFilter, expected: set[str]
) -> None:
    presentation = QualityPresenter.present(_report(), filter_value=filter_value)

    assert {row.code for row in presentation.rows} == expected
    assert tuple(row.code for row in presentation.rows) == tuple(
        code
        for code in ("quality.coordinate.nonfinite", "quality.coordinate.altloc", "quality.coordinate.info")
        if code in expected
    )


def test_quality_presenter_does_not_invent_score_or_causal_language() -> None:
    presentation = QualityPresenter.present(_report())
    rendered = repr(presentation).casefold()

    assert "druggability" not in rendered
    assert "pathogenic" not in rendered
    assert "causal" not in rendered
    assert "coordinate" in rendered
