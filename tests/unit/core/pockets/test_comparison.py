from __future__ import annotations

from collections.abc import Sequence
from dataclasses import FrozenInstanceError

import pytest

import structlens.core.pockets.comparison as comparison_module
from structlens.core.difference_maps import ResidueDisplacementVector
from structlens.core.evidence import Availability, Diagnostic, DiagnosticSeverity
from structlens.core.interactions import (
    InteractionChange,
    InteractionDifference,
    InteractionRecord,
    InteractionType,
    ReferenceInteractionKey,
)
from structlens.core.models import MutationEvent, MutationKind, ResidueCorrespondence, ResidueId
from structlens.core.pockets import AlphaSphere, PocketCandidate
from structlens.core.pockets.comparison import (
    compare_pocket_candidates,
    compare_pocket_match,
)
from structlens.core.pockets.matching import PocketMatch, PocketMatchState
from structlens.core.pockets.volume import compare_pocket_volumes, measure_pocket_volume
from structlens.core.pockets.volume_models import (
    PocketVolumeComparison,
    PocketVolumeResult,
    PocketVolumeSensitivity,
)


def _residue(number: int, structure: str, name: str = "ALA") -> ResidueId:
    return ResidueId(structure, "1", "A", str(number), None, name)


def _candidate(structure: str, lining: tuple[int, ...]) -> PocketCandidate:
    atom_ids = tuple(f"{structure}-atom-{offset}" for offset in range(4))
    return PocketCandidate(
        (
            AlphaSphere(
                (0.0, 0.0, 0.0), 3.0, atom_ids, tuple(_residue(number, structure) for number in lining), atom_ids
            ),
        )
    )


def test_comparison_reports_target_minus_reference_and_typed_evidence_changes() -> None:
    reference = _candidate("ref", (10, 20))
    target = _candidate("target", (10, 30))
    match = PocketMatch(
        reference_candidate=reference,
        target_candidate=target,
        state=PocketMatchState.MATCHED,
        lining_jaccard=1 / 3,
        centroid_distance_angstrom=1.25,
        score=0.5,
    )
    mutation = MutationEvent(
        0,
        MutationKind.SUBSTITUTION,
        _residue(20, "ref"),
        _residue(20, "target"),
        "A",
        "V",
        "20",
        "20",
        "A20V",
        0,
        64,
        "hydrophobic",
    )
    # Use a valid lost difference; the interaction contract requires a
    # reference record for a lost interaction.
    interaction = InteractionDifference(
        ReferenceInteractionKey(InteractionType.HBOND_GEOMETRIC, "A:10"),
        InteractionChange.LOST,
        reference_record=__import__("structlens.core.interactions", fromlist=["InteractionRecord"]).InteractionRecord(
            "ref", InteractionType.HBOND_GEOMETRIC, _residue(10, "ref"), None, "N", None, 3.0
        ),
    )
    unrelated_interaction = InteractionDifference(
        ReferenceInteractionKey(InteractionType.HBOND_GEOMETRIC, "A:99"),
        InteractionChange.LOST,
        reference_record=__import__("structlens.core.interactions", fromlist=["InteractionRecord"]).InteractionRecord(
            "ref", InteractionType.HBOND_GEOMETRIC, _residue(99, "ref"), None, "N", None, 3.0
        ),
    )
    result = compare_pocket_match(
        match,
        correspondences=tuple(
            ResidueCorrespondence(i, _residue(i + 10, "ref"), _residue(i + 10, "target"), "A", "A", "conserved")
            for i in (0, 10, 20)
        ),
        reference_volume=PocketVolumeComparison(Availability.AVAILABLE, 10.0, 0.5, 20.0, 10.0),
        target_volume=PocketVolumeComparison(Availability.AVAILABLE, 10.0, 0.5, 20.0, 10.0),
        reference_surface_area_angstrom2=40.0,
        target_surface_area_angstrom2=50.0,
        mutations=(mutation,),
        interaction_differences=(unrelated_interaction, interaction),
        local_displacements={_residue(10, "ref"): 2.0},
        qc_diagnostics=(Diagnostic("qc.warning", DiagnosticSeverity.WARNING, "local coordinate warning"),),
    )

    assert result.volume_delta_angstrom3 == -10.0
    assert result.surface_delta_angstrom2 == 10.0
    assert result.lining_residue_losses == (_residue(20, "ref"),)
    assert result.lining_residue_gains == (_residue(30, "target"),)
    assert result.associated_mutations == (mutation,)
    assert result.interaction_changes == (interaction,)
    serialized_interaction = result.to_json()["interaction_changes"][0]
    assert serialized_interaction["reference_record"]["distance_angstrom"] == 3.0
    assert serialized_interaction["key"]["reference_position_a"] == "A:10"
    assert result.local_displacement_angstrom == 2.0
    assert result.qc_diagnostics[0].code == "qc.warning"
    assert all("cause" not in str(result.to_json()).lower() for _ in [0])


def test_comparison_does_not_invent_relative_delta_for_zero_reference() -> None:
    candidate = _candidate("ref", (10,))
    match = PocketMatch(candidate, candidate, PocketMatchState.MATCHED, 1.0, 0.0, 1.0)
    result = compare_pocket_match(
        match,
        correspondences=(ResidueCorrespondence(0, _residue(10, "ref"), _residue(10, "target"), "A", "A", "conserved"),),
        reference_volume=PocketVolumeComparison(Availability.AVAILABLE, 4.0, None, 0.0, 0.0),
        target_volume=PocketVolumeComparison(Availability.AVAILABLE, 5.0, None, 0.0, 5.0),
    )
    assert result.volume_delta_angstrom3 == 5.0
    assert result.relative_volume_delta_fraction is None
    assert result.availability is Availability.AVAILABLE


def test_comparison_uses_authoritative_mapping_for_insertions_deletions_and_displacements() -> None:
    reference = _candidate("ref", (10, 20))
    target = _candidate("target", (10, 30))
    match = PocketMatch(reference, target, PocketMatchState.MATCHED)
    rows = (
        ResidueCorrespondence(
            0, _residue(10, "ref"), _residue(10, "target"), "A", "A", "conserved", ca_displacement_angstrom=1.0
        ),
        ResidueCorrespondence(1, _residue(20, "ref"), None, "A", None, "deletion", ca_displacement_angstrom=3.0),
        ResidueCorrespondence(2, None, _residue(30, "target"), None, "A", "insertion"),
    )
    result = compare_pocket_match(match, correspondences=rows)
    assert result.lining_residue_conserved == (_residue(10, "ref"),)
    assert result.lining_residue_losses == (_residue(20, "ref"),)
    assert result.lining_residue_gains == (_residue(30, "target"),)
    assert result.local_displacement_angstrom == 1.0


def test_comparison_propagates_qc_states_and_accepts_displacement_vectors() -> None:
    residue = _residue(10, "ref")
    target_residue = _residue(10, "target")
    match = PocketMatch(_candidate("ref", (10,)), _candidate("target", (10,)), PocketMatchState.AMBIGUOUS)
    vector = ResidueDisplacementVector(
        "A:10", residue, target_residue, (0.0, 0.0, 0.0), (0.0, 0.0, 2.0), (0.0, 0.0, 2.0), 2.0
    )
    warning = Diagnostic("qc.warning", DiagnosticSeverity.WARNING, "incomplete local QC")
    result = compare_pocket_match(
        match,
        correspondences=(ResidueCorrespondence(0, residue, target_residue, "A", "A", "conserved"),),
        local_displacements=(vector,),
        reference_qc=(warning,),
    )
    assert result.qc_availability is Availability.AVAILABLE
    assert result.local_displacement_angstrom == 2.0
    assert result.availability is Availability.AVAILABLE
    numerical = compare_pocket_match(match, qc=Availability.NUMERICAL_FAILURE)
    assert numerical.qc_availability is Availability.NUMERICAL_FAILURE
    with pytest.raises(TypeError):
        compare_pocket_match(match, local_displacements={"A:10": 1.0})  # type: ignore[dict-item]


def test_comparison_validates_contracts_and_volume_compatibility() -> None:
    candidate = _candidate("ref", (10,))
    match = PocketMatch(candidate, candidate, PocketMatchState.MATCHED)
    with pytest.raises(TypeError):
        compare_pocket_match(object())  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        compare_pocket_match(match, mutations=(object(),))  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        compare_pocket_match(match, interaction_differences=(object(),))  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        compare_pocket_match(match, reference_surface_area_angstrom2=-1.0, target_surface_area_angstrom2=1.0)
    with pytest.raises(TypeError):
        compare_pocket_match(match, qc_diagnostics=(object(),))  # type: ignore[arg-type]
    volume = PocketVolumeResult(
        Availability.AVAILABLE,
        1,
        1,
        1.0,
        1.0,
        (1, 1, 1),
        (0.0, 0.0, 0.0),
        sensitivity=PocketVolumeSensitivity(0.0, 0.0),
        compatibility_signature=(("grid", 1),),
    )
    incompatible = PocketVolumeResult(
        Availability.AVAILABLE,
        1,
        1,
        1.0,
        2.0,
        (1, 1, 1),
        (0.0, 0.0, 0.0),
        sensitivity=PocketVolumeSensitivity(1.0, 1.0),
        compatibility_signature=(("grid", 2),),
    )
    result = compare_pocket_match(match, reference_volume=volume, target_volume=incompatible)
    assert result.availability is Availability.NOT_APPLICABLE
    assert result.volume_delta_angstrom3 is None
    with pytest.raises(TypeError):
        compare_pocket_match(match, reference_volume=object())  # type: ignore[arg-type]


def test_comparison_wrapper_and_serialization_aliases_are_deterministic() -> None:
    reference = _candidate("ref", (10,))
    target = _candidate("target", (10,))
    rows = (ResidueCorrespondence(0, _residue(10, "ref"), _residue(10, "target"), "A", "A", "conserved"),)
    first = compare_pocket_candidates(reference, target, rows, reference_volume=2.0, target_volume=3.0)
    second = compare_pocket_candidates(reference, target, rows, reference_volume=2.0, target_volume=3.0)
    assert first == second
    assert first.volume_delta_angstrom3 == 1.0
    assert first.relative_volume_delta_fraction == 0.5
    with pytest.raises(TypeError):
        compare_pocket_candidates(object(), target)  # type: ignore[arg-type]
    with pytest.raises(FrozenInstanceError):
        first.surface_delta_angstrom2 = 2.0  # type: ignore[misc]


def test_comparison_validates_mapping_qc_and_displacement_input_shapes() -> None:
    candidate = _candidate("ref", (10,))
    match = PocketMatch(candidate, candidate, PocketMatchState.MATCHED)
    row = ResidueCorrespondence(0, _residue(10, "ref"), _residue(10, "target"), "A", "A", "conserved")
    with pytest.raises(TypeError):
        compare_pocket_match(match, correspondences={object(): object()})  # type: ignore[dict-item]
    duplicate_target = (
        row,
        ResidueCorrespondence(1, _residue(11, "ref"), _residue(10, "target"), "A", "A", "conserved"),
    )
    with pytest.raises(ValueError):
        compare_pocket_match(match, correspondences=duplicate_target)
    with pytest.raises(TypeError):
        compare_pocket_match(match, correspondences=(object(),))  # type: ignore[arg-type]
    for displacement in (
        {object(): 1.0},
        {_residue(10, "ref"): -1.0},
        ((object(), 1.0),),
        ((_residue(10, "ref"), float("nan")),),
        (object(),),
    ):
        with pytest.raises((TypeError, ValueError)):
            compare_pocket_match(match, local_displacements=displacement)  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        compare_pocket_match(match, qc=(object(),))  # type: ignore[arg-type]

    class BadQC:
        availability = "bad"
        diagnostics = ()

    with pytest.raises(ValueError):
        compare_pocket_match(match, qc=BadQC())
    assert compare_pocket_match(match, qc=Availability.NOT_APPLICABLE).qc_availability is Availability.NOT_APPLICABLE


def test_comparison_rejects_duplicate_reference_in_authoritative_mapping() -> None:
    candidate = _candidate("ref", (10,))
    match = PocketMatch(candidate, candidate, PocketMatchState.MATCHED)
    row = ResidueCorrespondence(0, _residue(10, "ref"), _residue(10, "target"), "A", "A", "conserved")
    duplicate = (
        row,
        ResidueCorrespondence(1, _residue(10, "ref"), _residue(11, "target"), "A", "A", "conserved"),
    )
    with pytest.raises(ValueError, match="duplicate reference"):
        compare_pocket_match(match, correspondences=duplicate)


def test_comparison_is_unavailable_without_authoritative_correspondence_or_for_unmatched() -> None:
    reference = _candidate("ref", (10,))
    target = _candidate("target", (10,))
    matched_without_map = compare_pocket_match(PocketMatch(reference, target, PocketMatchState.MATCHED))
    assert matched_without_map.availability is Availability.NOT_APPLICABLE
    assert matched_without_map.lining_residue_conserved == ()
    assert matched_without_map.lining_residue_gains == ()
    assert matched_without_map.lining_residue_losses == ()

    unmatched = compare_pocket_match(PocketMatch(reference, None, PocketMatchState.UNMATCHED_REFERENCE))
    assert unmatched.lining_residue_losses == ()
    assert unmatched.local_displacement_angstrom is None


def test_comparison_does_not_measure_displacement_for_deleted_target() -> None:
    reference = _candidate("ref", (10,))
    target = _candidate("target", (11,))
    match = PocketMatch(reference, target, PocketMatchState.MATCHED)
    rows = (
        ResidueCorrespondence(
            0,
            _residue(10, "ref"),
            None,
            "A",
            None,
            "deletion",
            ca_displacement_angstrom=9.0,
        ),
    )
    result = compare_pocket_match(match, correspondences=rows)
    assert result.local_displacement_angstrom is None
    assert (
        compare_pocket_match(
            match,
            correspondences=rows,
            local_displacements={_residue(10, "ref"): 9.0},
        ).local_displacement_angstrom
        is None
    )


def test_comparison_keeps_ca_and_sidechain_displacement_channels_separate() -> None:
    reference = _candidate("ref", (10,))
    target = _candidate("target", (10,))
    match = PocketMatch(reference, target, PocketMatchState.MATCHED)
    row = ResidueCorrespondence(
        0,
        _residue(10, "ref"),
        _residue(10, "target"),
        "A",
        "A",
        "conserved",
        ca_displacement_angstrom=2.0,
        sidechain_rmsd_angstrom=4.0,
    )
    result = compare_pocket_match(match, correspondences=(row,))
    assert result.ca_displacement_angstrom == 2.0
    assert result.sidechain_displacement_angstrom == 4.0


def test_comparison_canonicalizes_interaction_order_with_record_tiebreakers() -> None:
    reference = _candidate("ref", (10,))
    target = _candidate("target", (10,))
    match = PocketMatch(reference, target, PocketMatchState.MATCHED)
    row = ResidueCorrespondence(0, _residue(10, "ref"), _residue(10, "target"), "A", "A", "conserved")
    key = ReferenceInteractionKey(InteractionType.HBOND_GEOMETRIC, "A:10")
    near = InteractionDifference(
        key,
        InteractionChange.LOST,
        reference_record=InteractionRecord(
            "ref", InteractionType.HBOND_GEOMETRIC, _residue(10, "ref"), None, "N", None, 2.0
        ),
    )
    far = InteractionDifference(
        key,
        InteractionChange.LOST,
        reference_record=InteractionRecord(
            "ref", InteractionType.HBOND_GEOMETRIC, _residue(10, "ref"), None, "O", None, 3.0
        ),
    )

    first = compare_pocket_match(match, correspondences=(row,), interaction_differences=(far, near))
    second = compare_pocket_match(match, correspondences=(row,), interaction_differences=(near, far))

    assert first.interaction_changes == second.interaction_changes
    assert first.to_json()["interaction_changes"] == second.to_json()["interaction_changes"]


def test_comparison_validates_unilateral_surface_and_preserves_available_side() -> None:
    reference = _candidate("ref", (10,))
    target = _candidate("target", (10,))
    row = ResidueCorrespondence(0, _residue(10, "ref"), _residue(10, "target"), "A", "A", "conserved")
    result = compare_pocket_match(
        PocketMatch(reference, target, PocketMatchState.MATCHED),
        correspondences=(row,),
        reference_surface_area_angstrom2=12.0,
    )
    assert result.reference_surface_area_angstrom2 == 12.0
    assert result.target_surface_area_angstrom2 is None
    assert result.surface_delta_angstrom2 is None
    assert any(item.code == "pocket.surface.incomplete" for item in result.diagnostics)


def test_comparison_boundary_helpers_cover_native_missing_and_invalid_states() -> None:
    candidate = _candidate("ref", (10,))
    target = _candidate("target", (10,))
    reference = _residue(10, "ref")
    target_residue = _residue(10, "target")
    row = ResidueCorrespondence(0, reference, target_residue, "A", "A", "conserved")
    with pytest.raises(ValueError):
        comparison_module._finite("bad", "value")
    with pytest.raises(ValueError):
        comparison_module._finite(float("nan"), "value")
    assert comparison_module._correspondence_maps({reference: target_residue})[0][reference] == target_residue
    with pytest.raises(ValueError):
        comparison_module._correspondence_maps({reference: target_residue, _residue(11, "ref"): target_residue})
    assert comparison_module._volume_from_input(None, "reference") is None
    assert comparison_module._volume_from_input(1.0, "reference") == 1.0
    with pytest.raises(TypeError):
        comparison_module._volume_from_input(object(), "reference")  # type: ignore[arg-type]
    assert comparison_module._surface_delta(None, None)[0] is None
    assert comparison_module._surface_delta(1.0, None)[2] == 1.0
    assert comparison_module._normalise_displacements(((reference, 2.0),)) == ((reference, 2.0),)
    with pytest.raises(TypeError):
        comparison_module._normalise_displacements(object())

    class Oversized(Sequence[object]):
        def __len__(self) -> int:
            return comparison_module.MAX_COMPARISON_ITEMS + 1

        def __getitem__(self, index: int) -> object:
            raise AssertionError("oversized helper input must not be materialized")

    with pytest.raises(ValueError):
        comparison_module._correspondence_maps(Oversized())  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        comparison_module._normalise_displacements(Oversized())
    assert comparison_module._volume_delta(None, None).availability is Availability.NOT_APPLICABLE
    assert comparison_module._qc_values(None) == (None, ())
    assert comparison_module._qc_values(Availability.NOT_APPLICABLE)[0] is Availability.NOT_APPLICABLE

    class TypedQC:
        availability = "available"
        diagnostics = ()

    assert comparison_module._qc_values(TypedQC())[0] is Availability.AVAILABLE
    with pytest.raises(TypeError):
        comparison_module._qc_values(object())
    with pytest.raises(TypeError):
        comparison_module.PocketSurfaceMeasurement("a", "b", "c", 1.0, "m", object())
    with pytest.raises(ValueError):
        comparison_module.PocketSurfaceMeasurement("a", "b", "c", 1.0, "m", {"x": "y"}, {"surface_area": "angstrom"})
    unavailable = PocketVolumeResult(
        Availability.INVALID_INPUT, diagnostics=(Diagnostic("invalid", DiagnosticSeverity.ERROR, "invalid"),)
    )
    result = compare_pocket_match(
        PocketMatch(candidate, target, PocketMatchState.MATCHED),
        correspondences=(row,),
        reference_volume=unavailable,
        target_volume=unavailable,
    )
    assert result.volume is not None
    assert result.volume.availability is Availability.INVALID_INPUT
    propagated = compare_pocket_match(
        PocketMatch(candidate, target, PocketMatchState.MATCHED),
        correspondences=(row,),
        reference_volume=unavailable,
    )
    assert propagated.volume is not None
    assert propagated.volume.availability is Availability.INVALID_INPUT

    available_volume = PocketVolumeResult(
        Availability.AVAILABLE,
        1,
        1,
        1.0,
        1.0,
        (1, 1, 1),
        (0.0, 0.0, 0.0),
        sensitivity=PocketVolumeSensitivity(0.0, 0.0),
    )
    incomplete = compare_pocket_match(
        PocketMatch(candidate, target, PocketMatchState.MATCHED),
        correspondences=(row,),
        reference_volume=available_volume,
    )
    assert incomplete.volume is not None
    assert incomplete.volume.availability is Availability.NOT_APPLICABLE


def test_comparison_qc_error_diagnostics_are_invalid_input_and_object_is_rejected() -> None:
    candidate = _candidate("ref", (10,))
    match = PocketMatch(candidate, candidate, PocketMatchState.MATCHED)
    error = Diagnostic("qc.error", DiagnosticSeverity.ERROR, "invalid local coordinates")
    result = compare_pocket_match(match, qc=(error,))
    assert result.qc_availability is Availability.INVALID_INPUT
    assert result.availability is Availability.INVALID_INPUT
    with pytest.raises(TypeError):
        compare_pocket_match(match, qc=object())


def test_comparison_qc_diagnostic_argument_preserves_error_as_invalid_input() -> None:
    reference = _candidate("ref", (10,))
    target = _candidate("target", (10,))
    match = PocketMatch(reference, target, PocketMatchState.MATCHED)
    row = ResidueCorrespondence(0, _residue(10, "ref"), _residue(10, "target"), "A", "A", "conserved")
    error = Diagnostic("qc.error", DiagnosticSeverity.ERROR, "invalid local coordinates")

    result = compare_pocket_match(match, correspondences=(row,), qc_diagnostics=(error,))

    assert result.qc_availability is Availability.INVALID_INPUT
    assert result.availability is Availability.INVALID_INPUT
    assert result.qc_diagnostics == (error,)


def test_compare_real_candidate_volumes_with_equal_settings_and_reject_incompatible_comparison() -> None:
    reference = _candidate("ref", (10,))
    target = _candidate("target", (10,))
    settings = __import__(
        "structlens.core.pockets.volume_models", fromlist=["PocketVolumeSettings"]
    ).PocketVolumeSettings(
        coarse_grid_spacing_angstrom=2.0,
        fine_grid_spacing_angstrom=1.0,
    )
    reference_result = measure_pocket_volume(reference, settings=settings)
    target_result = measure_pocket_volume(target, settings=settings)
    assert reference_result.compatibility_signature == target_result.compatibility_signature
    volume = compare_pocket_volumes(reference_result, target_result)
    assert volume.availability is Availability.AVAILABLE
    incompatible = PocketVolumeComparison(
        Availability.AVAILABLE, 1.0, 0.1, 10.0, 11.0, compatibility_signature=(("grid", 2),)
    )
    comparison = compare_pocket_match(
        PocketMatch(reference, target, PocketMatchState.MATCHED),
        correspondences=(ResidueCorrespondence(0, _residue(10, "ref"), _residue(10, "target"), "A", "A", "conserved"),),
        reference_volume=volume,
        target_volume=incompatible,
    )
    assert comparison.volume_delta_angstrom3 is None
    assert comparison.availability is Availability.NOT_APPLICABLE


def test_comparison_of_volume_comparisons_preserves_sensitivity_and_provenance() -> None:
    reference = _candidate("ref", (10,))
    target = _candidate("target", (10,))
    row = ResidueCorrespondence(0, _residue(10, "ref"), _residue(10, "target"), "A", "A", "conserved")
    signature = (("method", "pocket_free_volume"), ("grid", 1))
    reference_measure = PocketVolumeComparison(
        Availability.AVAILABLE,
        2.0,
        0.2,
        10.0,
        12.0,
        compatibility_signature=signature,
        reference_sensitivity=PocketVolumeSensitivity(1.0, 0.1),
    )
    target_measure = PocketVolumeComparison(
        Availability.AVAILABLE,
        3.0,
        0.3,
        20.0,
        23.0,
        compatibility_signature=signature,
        target_sensitivity=PocketVolumeSensitivity(2.0, 0.1),
    )

    result = compare_pocket_match(
        PocketMatch(reference, target, PocketMatchState.MATCHED),
        correspondences=(row,),
        reference_volume=reference_measure,
        target_volume=target_measure,
    )

    assert result.volume is not None
    assert result.volume.reference_sensitivity == reference_measure.reference_sensitivity
    assert result.volume.target_sensitivity == target_measure.target_sensitivity
