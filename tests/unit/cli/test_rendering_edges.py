from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from structlens.application.pocket_comparison_service import PocketComparisonReport, PocketComparisonService
from structlens.application.pocket_service import PocketDetectionReport
from structlens.cli import rendering
from structlens.cli.commands import _interaction_differences, run_compare
from structlens.core.evidence import Availability
from structlens.core.models import ResidueCorrespondence, ResidueId
from structlens.core.pockets import (
    AlphaSphere,
    PocketCandidate,
    PocketMatch,
    PocketMatchingResult,
    PocketMatchingSettings,
    PocketMatchState,
)
from structlens.core.pockets.volume_models import PocketVolumeComparison, PocketVolumeResult


def _candidate(structure: str = "structure") -> PocketCandidate:
    residue = ResidueId(structure, "1", "A", "1", None, "ALA")
    atoms = ("atom-a", "atom-b", "atom-c", "atom-d")
    sphere = AlphaSphere((0.0, 0.0, 0.0), 3.0, atoms, (residue,), atoms)
    return PocketCandidate((sphere,), source_content_id="a" * 64, selection_id="b" * 64)


def test_rendering_handles_empty_analysis_and_matched_pocket() -> None:
    result = run_compare(
        Path("tests/fixtures/parsing/numbering.mmcif"),
        Path("tests/fixtures/parsing/numbering.mmcif"),
        mode="sequence",
    )
    assert result.report is not None
    unavailable = replace(
        result.report,
        analysis=None,
        availability=replace(result.report.availability, analysis=Availability.NOT_APPLICABLE),
    )
    assert "Analysis: unavailable" in rendering.render_compare(unavailable)
    candidate = _candidate("reference")
    target = _candidate("target")
    comparison = PocketComparisonService(PocketMatchingSettings(minimum_lining_jaccard=0.0)).compare(
        (candidate,),
        (target,),
        (
            ResidueCorrespondence(
                0,
                ResidueId("reference", "1", "A", "1", None, "ALA"),
                ResidueId("target", "1", "A", "1", None, "ALA"),
                "A",
                "A",
                "conserved",
            ),
        ),
    )
    text = rendering.render_pocket_compare(comparison)
    assert "matched" in text
    assert "Score:" in text

    measured = replace(
        comparison.comparisons[0],
        volume=PocketVolumeComparison(
            Availability.AVAILABLE,
            delta_angstrom3=1.25,
            reference_volume_angstrom3=10.0,
            target_volume_angstrom3=11.25,
        ),
        surface_delta_angstrom2=2.5,
        ca_displacement_angstrom=0.75,
        local_displacement_angstrom=1.5,
    )
    rendered_measures = rendering.render_pocket_compare(replace(comparison, comparisons=(measured,)))
    assert "Volume delta: 1.250" in rendered_measures
    assert "Surface delta: 2.500" in rendered_measures
    assert "Cα displacement: 0.750" in rendered_measures
    assert "Local displacement: 1.500" in rendered_measures


def test_rendering_lists_candidates_and_stable_json() -> None:
    candidate = _candidate()
    report = PocketDetectionReport(Availability.AVAILABLE, candidates=(candidate,))
    text = rendering.render_pockets(report)
    assert candidate.candidate_id in text
    assert rendering.json_bytes({"z": 1, "a": 2}) == b'{"a":2,"z":1}'


def test_interaction_projection_accepts_legacy_typed_sequence() -> None:
    result = run_compare(
        Path("tests/fixtures/parsing/numbering.mmcif"),
        Path("tests/fixtures/parsing/numbering.mmcif"),
        mode="sequence",
    )
    assert result.report is not None
    changed = replace(result.report, interactions=())
    assert _interaction_differences(changed) == ()


def test_rendering_exposes_ambiguous_and_both_unmatched_typed_states() -> None:
    reference = _candidate("reference")
    target = _candidate("target")
    matching = PocketMatchingResult(
        Availability.AVAILABLE,
        matches=(
            PocketMatch(reference, target, PocketMatchState.AMBIGUOUS, score=0.4, alternative_score=0.35),
            PocketMatch(reference, None, PocketMatchState.UNMATCHED_REFERENCE),
            PocketMatch(None, target, PocketMatchState.UNMATCHED_TARGET),
        ),
    )
    report = PocketComparisonReport(Availability.AVAILABLE, matching)
    reference_volume = PocketVolumeResult(Availability.NOT_APPLICABLE, candidate_id=reference.candidate_id)
    target_volume = PocketVolumeResult(Availability.NOT_APPLICABLE, candidate_id=target.candidate_id)
    text = rendering.render_pocket_compare(
        report,
        reference_volumes=(reference_volume,),
        target_volumes=(target_volume,),
    )
    assert "ambiguous" in text
    assert "unmatched_reference" in text
    assert "unmatched_target" in text
    assert "Reference volumes: not_applicable" in text
    assert "Target volumes: not_applicable" in text
    assert reference.candidate_id in text
    assert target.candidate_id in text
