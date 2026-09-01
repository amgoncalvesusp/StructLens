from __future__ import annotations

from pathlib import Path
from typing import Any

from structlens.application.pocket_comparison_service import PocketComparisonReport
from structlens.application.pocket_service import PocketDetectionReport
from structlens.cli.commands import CliServices, run_pocket_compare, run_pockets
from structlens.core.evidence import Availability, Diagnostic, DiagnosticSeverity
from structlens.core.models import ResidueId
from structlens.core.pockets import AlphaSphere, PocketCandidate, PocketMatchingResult
from structlens.core.pockets.volume_models import PocketVolumeResult

_FIXTURE = Path("tests/fixtures/parsing/numbering.mmcif")


def _candidate(parsed: Any) -> PocketCandidate:
    residue = ResidueId("fixture", "1", "A", "1", None, "ALA")
    atom_ids = ("a", "b", "c", "d")
    sphere = AlphaSphere((0.0, 0.0, 0.0), 3.0, atom_ids, (residue,), atom_ids)
    return PocketCandidate(
        (sphere,),
        source_content_id=parsed.selection.content_id,
        selection_id=parsed.selection.selection_id,
    )


def _diagnostic(code: str) -> Diagnostic:
    return Diagnostic(code, DiagnosticSeverity.ERROR, f"{code} diagnostic")


class _FatalDetectionService:
    def analyze(self, parsed: Any) -> PocketDetectionReport:
        return PocketDetectionReport(Availability.NUMERICAL_FAILURE, diagnostics=(_diagnostic("detect.failed"),))

    def measure_volume(self, parsed: Any, candidate: Any) -> PocketVolumeResult:
        raise AssertionError("fatal detection must not request candidate volumes")


class _FatalVolumeService:
    def analyze(self, parsed: Any) -> PocketDetectionReport:
        candidate = _candidate(parsed)
        return PocketDetectionReport(Availability.AVAILABLE, candidates=(candidate,))

    def measure_volume(self, parsed: Any, candidate: PocketCandidate) -> PocketVolumeResult:
        return PocketVolumeResult(
            Availability.NUMERICAL_FAILURE,
            candidate_id=candidate.candidate_id,
            diagnostics=(_diagnostic("volume.failed"),),
        )


class _NoDetectionService:
    def analyze(self, parsed: Any) -> PocketDetectionReport:
        return PocketDetectionReport(Availability.NOT_DETECTED)

    def measure_volume(self, parsed: Any, candidate: Any) -> PocketVolumeResult:
        raise AssertionError("no-result detection must not request candidate volumes")


class _NeutralComparisonService:
    def compare(self, *args: Any, **kwargs: Any) -> PocketComparisonReport:
        return PocketComparisonReport(Availability.NOT_APPLICABLE, PocketMatchingResult(Availability.NOT_APPLICABLE))


class _FixedComparisonService:
    def __init__(self, report: PocketComparisonReport) -> None:
        self.report = report

    def compare(self, *args: Any, **kwargs: Any) -> PocketComparisonReport:
        return self.report


def test_pockets_retains_selection_candidates_volumes_and_lineage() -> None:
    result = run_pockets(_FIXTURE)
    assert "selection" in result.payload
    assert "pockets" in result.payload
    assert result.payload["pockets"]["detections"]
    detection = result.payload["pockets"]["detections"][0]
    if detection["candidates"]:
        candidate = detection["candidates"][0]
        assert candidate["source_content_id"] == result.payload["selection"]["content_id"]
        assert candidate["selection_id"] == result.payload["selection"]["selection_id"]
    assert "volumes" in result.payload["pockets"]


def test_pockets_no_detection_is_explicit_state() -> None:
    result = run_pockets(Path("tests/fixtures/parsing/numbering_altloc.pdb"))
    assert result.exit_code in {0, 2}
    assert result.status in {"available", "not_detected", "invalid_input", "numerical_failure"}
    assert "diagnostics" in result.payload["pockets"]


def test_pockets_propagates_fatal_detection_and_volume_states() -> None:
    fatal_detection = run_pockets(_FIXTURE, services=CliServices(pocket_service=_FatalDetectionService()))
    assert fatal_detection.status == Availability.NUMERICAL_FAILURE.value
    assert fatal_detection.exit_code != 0
    assert "detect.failed" in fatal_detection.summary

    fatal_volume = run_pockets(_FIXTURE, services=CliServices(pocket_service=_FatalVolumeService()))
    assert fatal_volume.status == Availability.NUMERICAL_FAILURE.value, fatal_volume.error
    assert fatal_volume.exit_code != 0
    assert "volume.failed" in fatal_volume.summary
    assert fatal_volume.payload["pockets"]["volumes"][0]["availability"] == Availability.NUMERICAL_FAILURE.value


def test_pockets_not_detected_remains_successful() -> None:
    result = run_pockets(_FIXTURE, services=CliServices(pocket_service=_NoDetectionService()))
    assert result.status == Availability.NOT_DETECTED.value
    assert result.exit_code == 0


def test_pocket_json_bytes_are_deterministic_with_typed_evidence() -> None:
    pockets_first = run_pockets(_FIXTURE, services=CliServices(pocket_service=_FatalVolumeService())).json_bytes()
    pockets_second = run_pockets(_FIXTURE, services=CliServices(pocket_service=_FatalVolumeService())).json_bytes()
    assert pockets_first == pockets_second

    compare_first = run_pocket_compare(
        _FIXTURE,
        _FIXTURE,
        services=CliServices(
            pocket_service=_FatalVolumeService(),
            comparison_service=_NeutralComparisonService(),
        ),
    ).json_bytes()
    compare_second = run_pocket_compare(
        _FIXTURE,
        _FIXTURE,
        services=CliServices(
            pocket_service=_FatalVolumeService(),
            comparison_service=_NeutralComparisonService(),
        ),
    ).json_bytes()
    assert compare_first == compare_second


def test_pocket_compare_retains_typed_match_and_both_inputs() -> None:
    result = run_pocket_compare(_FIXTURE, _FIXTURE)
    assert result.payload["reference_selection"]["content_id"] == result.payload["target_selection"]["content_id"]
    assert "pockets" in result.payload
    assert "matching" in result.payload["pockets"]
    assert "comparisons" in result.payload["pockets"]
    assert "concordance" in result.payload["pockets"]
    assert "analysis_report" in result.payload
    assert "comparison_report" in result.payload["pockets"]
    assert "units" in result.payload["pockets"]
    assert (
        result.payload["analysis_report"]["reference_selection"]["content_id"]
        == result.payload["reference_source"]["content_id"]
    )
    analysis_report = result.payload["analysis_report"]
    assert "input_quality" in analysis_report
    assert "availability" in analysis_report
    assert "diagnostics" in analysis_report
    assert "interactions" in analysis_report
    assert "displacement_vectors" in analysis_report
    assert "provenance" in analysis_report
    assert analysis_report["analysis"] is not None
    assert "mutations" in analysis_report["analysis"]
    assert "correspondences" in analysis_report["analysis"]
    assert "transform" in analysis_report["analysis"]
    comparison_report = result.payload["pockets"]["comparison_report"]
    assert comparison_report["availability"] in {item.value for item in Availability}
    for key in ("matching", "comparisons", "concordance"):
        assert comparison_report[key] == result.payload["pockets"][key]
    for diagnostic in comparison_report["diagnostics"]:
        assert diagnostic in result.payload["pockets"]["diagnostics"]
    for section, state in analysis_report["availability"].items():
        if analysis_report.get(section) is None:
            assert state != Availability.AVAILABLE.value


def test_pocket_compare_embeds_exact_typed_comparison_report() -> None:
    expected = PocketComparisonReport(Availability.NOT_APPLICABLE, PocketMatchingResult(Availability.NOT_APPLICABLE))
    result = run_pocket_compare(
        _FIXTURE,
        _FIXTURE,
        services=CliServices(
            pocket_service=_FatalDetectionService(),
            comparison_service=_FixedComparisonService(expected),
        ),
    )
    assert result.payload["pockets"]["comparison_report"] == expected.to_json()


def test_pocket_compare_invalid_selection_fails_without_silent_chain_choice() -> None:
    result = run_pocket_compare(_FIXTURE, _FIXTURE, reference_author_chains=("missing",))
    assert result.exit_code != 0
    assert result.error is not None


def test_pocket_commands_reject_boundary_and_manual_errors() -> None:
    fixture = _FIXTURE
    assert run_pockets(Path("missing-structure.pdb")).exit_code == 2
    assert run_pocket_compare(fixture, fixture, mode="manual").exit_code == 2


def test_pocket_compare_propagates_fatal_detection_and_volume_states() -> None:
    fatal_detection = run_pocket_compare(
        _FIXTURE,
        _FIXTURE,
        services=CliServices(pocket_service=_FatalDetectionService()),
    )
    assert fatal_detection.status == Availability.NUMERICAL_FAILURE.value
    assert fatal_detection.exit_code != 0
    assert any(item["code"] == "detect.failed" for item in fatal_detection.payload["pockets"]["diagnostics"])

    fatal_volume = run_pocket_compare(
        _FIXTURE,
        _FIXTURE,
        services=CliServices(
            pocket_service=_FatalVolumeService(),
            comparison_service=_NeutralComparisonService(),
        ),
    )
    assert fatal_volume.status == Availability.NUMERICAL_FAILURE.value, fatal_volume.error
    assert fatal_volume.exit_code != 0
    assert any(item["code"] == "volume.failed" for item in fatal_volume.payload["pockets"]["diagnostics"])


def test_application_type_errors_are_not_swallowed() -> None:
    class BrokenService:
        def analyze(self, parsed: Any) -> PocketDetectionReport:
            raise TypeError("programming fault")

    try:
        run_pockets(_FIXTURE, services=CliServices(pocket_service=BrokenService()))
    except TypeError as error:
        assert str(error) == "programming fault"
    else:
        raise AssertionError("application TypeError must remain visible")
