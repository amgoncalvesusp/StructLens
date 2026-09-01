"""RED contracts for the report-first GUI presentation boundary.

These tests intentionally stay independent from Qt.  The presenter is a
consumer of immutable report data: it may format values and attach method
text, but it must not run an analysis or manufacture scientific evidence.
"""

from __future__ import annotations

from dataclasses import asdict, fields, is_dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Any

import pytest

from structlens.application.dto import AnalysisReportRequest
from structlens.application.report_service import ReportService
from structlens.core.evidence import Availability
from structlens.core.models import AlignmentMode, AnalysisSettings, ResidueId
from structlens.core.msa import MSASettings
from structlens.core.parsing import InputSelection, StructureFormat, capture_snapshot
from structlens.core.sites import SiteDefinition, SiteDefinitionMode
from structlens.plugin.gui.presentation import present_report

_FIXTURES = Path(__file__).parents[2] / "fixtures" / "parsing"


def _request(tmp_path: Path) -> AnalysisReportRequest:
    source = (_FIXTURES / "enriched.pdb").read_bytes()
    reference_path = tmp_path / "reference.pdb"
    target_path = tmp_path / "target.pdb"
    reference_path.write_bytes(source)
    target_path.write_bytes(source)
    reference = capture_snapshot(reference_path)
    target = capture_snapshot(target_path)
    reference_selection = InputSelection(
        reference.content_id,
        reference.display_name,
        StructureFormat.PDB,
        "1",
        author_chain_ids=("A",),
        path=str(reference_path),
    )
    target_selection = InputSelection(
        target.content_id,
        target.display_name,
        StructureFormat.PDB,
        "1",
        author_chain_ids=("A",),
        path=str(target_path),
    )
    return AnalysisReportRequest(
        reference,
        target,
        reference_selection,
        target_selection,
        analysis_settings=AnalysisSettings(alignment_mode=AlignmentMode.SEQUENCE),
        msa_settings=MSASettings(),
        site_definitions=(
            SiteDefinition(
                "active",
                "Active site",
                SiteDefinitionMode.KEY_RESIDUES,
                (ResidueId("reference", "1", "A", "100", None, "ALA"),),
            ),
        ),
    )


def _field(value: object, name: str) -> object:
    if isinstance(value, dict) or isinstance(value, MappingProxyType):
        return value[name]
    return getattr(value, name)


def _section(presentation: object, name: str) -> object:
    sections = _field(presentation, "sections")
    return _field(sections, name)


def _status_value(section: object) -> str:
    value = _field(section, "status")
    return str(getattr(value, "value", value))


def _text(value: object) -> str:
    """Flatten a frozen presentation value for wording assertions."""

    if isinstance(value, str):
        return value
    if isinstance(value, MappingProxyType) or isinstance(value, dict):
        return " ".join(f"{key} {_text(item)}" for key, item in value.items())
    if isinstance(value, (tuple, list)):
        return " ".join(_text(item) for item in value)
    if is_dataclass(value):
        return _text(asdict(value))
    return str(value)


def test_present_report_projects_every_v03_section_from_one_report(tmp_path: Path) -> None:
    report = ReportService().analyze(_request(tmp_path))

    presentation = present_report(report)

    assert _field(presentation, "report_id") == report.report_id
    for name in (
        "summary",
        "quality",
        "msa",
        "mutations",
        "correspondences",
        "interactions",
        "sites",
        "distance_map",
        "vectors",
        "evidence_cards",
        "diagnostics",
        "export",
        "pymol",
    ):
        section = _section(presentation, name)
        assert section is not None, name

    assert _status_value(_section(presentation, "msa")) == Availability.AVAILABLE.value
    assert _status_value(_section(presentation, "sites")) == Availability.AVAILABLE.value
    assert _status_value(_section(presentation, "distance_map")) == Availability.AVAILABLE.value
    assert _status_value(_section(presentation, "vectors")) == Availability.AVAILABLE.value


def test_present_report_has_explicit_units_and_method_limits(tmp_path: Path) -> None:
    presentation = present_report(ReportService().analyze(_request(tmp_path)))

    units = _field(presentation, "units")
    assert units["length"] == "Å"
    assert units["area"] == "Å²"
    assert units["volume"] == "Å³"
    assert units["fraction"] == "fraction"

    explanations = _field(presentation, "explanations")
    assert "cα" in _text(explanations).casefold()
    assert "atomic envelope" in _text(explanations).casefold()
    assert "pocket free volume" in _text(explanations).casefold()
    assert "method" in _text(explanations).casefold() or "limit" in _text(explanations).casefold()


def test_present_report_is_deterministic_immutable_and_non_causal(tmp_path: Path) -> None:
    report = ReportService().analyze(_request(tmp_path))

    first = present_report(report)
    second = present_report(report)

    assert first == second
    assert is_dataclass(first)
    with pytest.raises((AttributeError, TypeError)):
        first.report_id = "tampered"

    rendered = _text(first).casefold()
    assert "druggability score" not in rendered
    assert "binding-affinity score" not in rendered
    assert "pathogenicity" not in rendered
    assert "causal claim" not in rendered


def test_none_evidence_is_unavailable_with_native_reason_not_zero(tmp_path: Path) -> None:
    class FailingMSA:
        def align(self, *_args: object, **_kwargs: object) -> object:
            raise RuntimeError("MSA backend is unavailable")

    report = ReportService(msa_engine=FailingMSA()).analyze(_request(tmp_path))
    assert report.availability.msa is Availability.NUMERICAL_FAILURE
    assert report.msa is None

    msa = _section(present_report(report), "msa")
    rendered = _text(msa)
    assert "Unavailable —" in rendered
    # The canonical report deliberately contains the typed public diagnostic,
    # not the backend exception text.  Presentation must preserve that native
    # reason without inventing or leaking hidden exception details.
    assert "report.msa.failed" in rendered
    assert "could not be calculated" in rendered
    assert "0.000" not in rendered
    assert _status_value(msa) == Availability.NUMERICAL_FAILURE.value


def test_optional_failure_does_not_hide_independent_report_sections(tmp_path: Path) -> None:
    class FailingMSA:
        def align(self, *_args: object, **_kwargs: object) -> object:
            raise RuntimeError("alignment backend failed")

    presentation = present_report(ReportService(msa_engine=FailingMSA()).analyze(_request(tmp_path)))

    assert _status_value(_section(presentation, "msa")) == Availability.NUMERICAL_FAILURE.value
    assert _status_value(_section(presentation, "interactions")) == Availability.AVAILABLE.value
    assert _status_value(_section(presentation, "distance_map")) == Availability.AVAILABLE.value
    assert _status_value(_section(presentation, "vectors")) == Availability.AVAILABLE.value
    assert _status_value(_section(presentation, "evidence_cards")) == Availability.AVAILABLE.value
    diagnostics = _text(_section(presentation, "diagnostics"))
    assert "report.msa.failed" in diagnostics


def test_presentation_values_are_frozen_slots_without_mutable_report_aliases(tmp_path: Path) -> None:
    presentation = present_report(ReportService().analyze(_request(tmp_path)))

    assert getattr(type(presentation), "__slots__", None)
    assert not hasattr(presentation, "payload")
    for field in fields(presentation):
        value: Any = getattr(presentation, field.name)
        if isinstance(value, dict):
            pytest.fail(f"presentation field {field.name} exposes a mutable dict")
