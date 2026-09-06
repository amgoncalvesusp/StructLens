"""RED contracts for the typed Sites & Pockets GUI page."""

from __future__ import annotations

from types import MappingProxyType

from structlens.plugin.gui.pockets_page import PocketPresenter

from structlens.core.evidence import Availability, Diagnostic, DiagnosticSeverity
from structlens.core.models import ResidueId
from structlens.core.pockets import AlphaSphere, PocketCandidate
from structlens.core.pockets.matching import PocketMatch, PocketMatchingResult, PocketMatchState
from structlens.core.pockets.volume_models import PocketVolumeResult, PocketVolumeSensitivity
from structlens.core.reports.pockets import PocketDetectionSnapshot, PocketReportSnapshot, PocketVolumeSnapshot


def _candidate(structure_id: str, offset: float = 0.0) -> PocketCandidate:
    residue = ResidueId(structure_id, "1", "A", "42", None, "ALA")
    spheres = tuple(
        AlphaSphere(
            (offset + x, y, z),
            3.0,
            (f"{structure_id}-a{i}", f"{structure_id}-b{i}", f"{structure_id}-c{i}", f"{structure_id}-d{i}"),
            (residue,),
            (f"{structure_id}-s{i}", f"{structure_id}-t{i}", f"{structure_id}-u{i}", f"{structure_id}-v{i}"),
        )
        for i, (x, y, z) in enumerate(((0.0, 0.0, 0.0), (0.8, 0.0, 0.0)))
    )
    return PocketCandidate(spheres)


def _snapshot() -> tuple[PocketReportSnapshot, PocketCandidate, PocketCandidate, PocketVolumeResult]:
    reference = _candidate("reference")
    target = _candidate("target", 0.5)
    matched = PocketMatch(reference, target, PocketMatchState.MATCHED, 0.8, 0.5, 0.9)
    ambiguous = PocketMatch(reference, target, PocketMatchState.AMBIGUOUS, 0.7, 0.8, 0.7)
    unmatched_reference = PocketMatch(reference, None, PocketMatchState.UNMATCHED_REFERENCE)
    unmatched_target = PocketMatch(None, target, PocketMatchState.UNMATCHED_TARGET)
    matching = PocketMatchingResult(
        Availability.AVAILABLE,
        matches=(matched, ambiguous, unmatched_reference, unmatched_target),
    )
    volume = PocketVolumeResult(
        Availability.AVAILABLE,
        coarse_voxel_count=100,
        fine_voxel_count=800,
        coarse_volume_angstrom3=110.0,
        fine_volume_angstrom3=104.0,
        coarse_grid_shape=(10, 10, 10),
        fine_grid_shape=(20, 20, 20),
        sensitivity=PocketVolumeSensitivity(6.0, 6.0 / 104.0),
        candidate_id=reference.candidate_id,
        units=MappingProxyType({"length": "angstrom", "volume": "angstrom^3"}),
    )
    snapshot = PocketReportSnapshot(
        detections=(
            PocketDetectionSnapshot("reference", Availability.AVAILABLE, candidates=(reference,)),
            PocketDetectionSnapshot("target", Availability.AVAILABLE, candidates=(target,)),
        ),
        volumes=(PocketVolumeSnapshot("reference", volume),),
        matching=matching,
        diagnostics=(
            Diagnostic(
                "pocket.presenter.note",
                DiagnosticSeverity.INFO,
                "Volume is a dual-resolution free-volume estimate.",
            ),
        ),
    )
    return snapshot, reference, target, volume


def test_pocket_presenter_lists_candidates_and_explicit_match_states() -> None:
    snapshot, reference, target, _ = _snapshot()

    presentation = PocketPresenter.present(snapshot)

    assert [row.state.value for row in presentation.match_rows] == [
        "matched",
        "ambiguous",
        "unmatched_reference",
        "unmatched_target",
    ]
    assert presentation.candidate_rows
    assert {row.candidate_id for row in presentation.candidate_rows} == {
        reference.candidate_id,
        target.candidate_id,
    }
    assert all(row.state.value != "unavailable" for row in presentation.match_rows)


def test_pocket_presenter_detail_drawer_retains_volume_units_and_coarse_fine_sensitivity() -> None:
    snapshot, reference, _, volume = _snapshot()

    presentation = PocketPresenter.present(snapshot)
    detail = presentation.detail_for(reference.candidate_id)

    assert detail is not None
    assert detail.candidate_id == reference.candidate_id
    assert detail.units["volume"] == "Å³"
    assert detail.coarse_volume_angstrom3 == volume.coarse_volume_angstrom3
    assert detail.fine_volume_angstrom3 == volume.fine_volume_angstrom3
    assert detail.absolute_sensitivity_angstrom3 == volume.absolute_sensitivity_angstrom3
    assert detail.relative_sensitivity == volume.relative_sensitivity
    assert detail.method_explanation
    assert "druggability" not in repr(detail).casefold()
    assert "binding-affinity" not in repr(detail).casefold()


def test_pocket_presenter_actions_are_disabled_when_no_capability_is_available() -> None:
    empty = PocketReportSnapshot(
        diagnostics=(
            Diagnostic(
                "pocket.detect.unavailable",
                DiagnosticSeverity.WARNING,
                "Pocket detector is unavailable for this input.",
            ),
        ),
        availability=Availability.DEPENDENCY_UNAVAILABLE,
    )

    presentation = PocketPresenter.present(empty)

    assert presentation.can_detect is False
    assert presentation.can_measure is False
    assert presentation.can_compare is False
    assert presentation.disabled_reasons
    assert "druggability" not in repr(presentation).casefold()
