from __future__ import annotations

from collections.abc import Sequence
from dataclasses import FrozenInstanceError

import pytest

import structlens.core.provenance as provenance_module
from structlens.core.evidence import Availability, Diagnostic, DiagnosticSeverity
from structlens.core.models import ResidueCorrespondence, ResidueId, StructuralTransform
from structlens.core.pockets import AlphaSphere, PocketCandidate
from structlens.core.pockets.matching import (
    PocketMatch,
    PocketMatchingSettings,
    PocketMatchState,
    match_pocket_candidates,
)


def _residue(number: int, structure: str) -> ResidueId:
    return ResidueId(structure, "1", "A", str(number), None, "ALA")


def _candidate(
    index: int, structure: str, lining: tuple[int, ...], center: tuple[float, float, float]
) -> PocketCandidate:
    atom_ids = tuple(f"{structure}-{index}-atom-{offset}" for offset in range(4))
    return PocketCandidate(
        (
            AlphaSphere(
                center_xyz=center,
                radius_angstrom=3.0,
                touching_atom_ids=atom_ids,
                lining_residues=tuple(_residue(number, structure) for number in lining),
                source_simplex_atom_ids=atom_ids,
            ),
        )
    )


def _correspondences(numbers: tuple[int, ...]) -> tuple[ResidueCorrespondence, ...]:
    return tuple(
        ResidueCorrespondence(
            index,
            _residue(number, "ref"),
            _residue(number, "target"),
            "A",
            "A",
            "conserved",
        )
        for index, number in enumerate(numbers)
    )


def test_pocket_match_bounds_units_before_copy_and_rejects_oversized_strings() -> None:
    candidate = _candidate(1, "ref", (10,), (0.0, 0.0, 0.0))

    class OversizedUnits(dict[str, str]):
        def __len__(self) -> int:
            return provenance_module.MAX_UNIT_MAP_ITEMS + 1

        def items(self) -> object:
            raise AssertionError("oversized units must be rejected before iteration")

    with pytest.raises(ValueError, match="units.*item limit"):
        PocketMatch(candidate, None, PocketMatchState.UNMATCHED_REFERENCE, units=OversizedUnits())

    huge = "x" * (provenance_module.MAX_UNIT_STRING_LENGTH + 1)
    for units in ({huge: "count"}, {"measure": huge}):
        with pytest.raises(ValueError, match="units.*string length"):
            PocketMatch(candidate, None, PocketMatchState.UNMATCHED_REFERENCE, units=units)


def test_matching_uses_authoritative_lining_overlap_and_transformed_centroid() -> None:
    reference = _candidate(1, "ref", (10, 20, 30), (0.0, 0.0, 0.0))
    target = _candidate(1, "target", (10, 20, 40), (5.0, 0.0, 0.0))

    result = match_pocket_candidates(
        (reference,),
        (target,),
        _correspondences((10, 20, 30, 40)),
        transform=StructuralTransform(translation=(-5.0, 0.0, 0.0)),
        settings=PocketMatchingSettings(max_centroid_distance_angstrom=0.1),
    )

    match = result.matches[0]
    assert match.state is PocketMatchState.MATCHED
    assert match.lining_jaccard == 2 / 4
    assert match.centroid_distance_angstrom == 0.0


def test_matching_retains_ambiguity_and_unmatched_states_independent_of_input_order() -> None:
    reference = _candidate(1, "ref", (10,), (0.0, 0.0, 0.0))
    target_a = _candidate(1, "target", (10,), (0.0, 0.0, 0.0))
    target_b = _candidate(2, "target", (10,), (0.0, 0.0, 0.0))
    target_extra = _candidate(3, "target", (90,), (20.0, 0.0, 0.0))

    result = match_pocket_candidates(
        (reference,),
        (target_extra, target_b, target_a),
        _correspondences((10,)),
        settings=PocketMatchingSettings(ambiguity_margin=0.0, max_centroid_distance_angstrom=4.0),
    )

    assert result.matches[0].state is PocketMatchState.AMBIGUOUS
    assert any(item.state is PocketMatchState.UNMATCHED_TARGET for item in result.matches)
    assert (
        result.matches
        == match_pocket_candidates(
            (reference,),
            (target_a, target_b, target_extra),
            _correspondences((10,)),
            settings=PocketMatchingSettings(ambiguity_margin=0.0, max_centroid_distance_angstrom=4.0),
        ).matches
    )


def test_matching_handles_insertions_deletions_and_empty_collections_without_fabricating_pairs() -> None:
    reference = _candidate(1, "ref", (10,), (0.0, 0.0, 0.0))
    target = _candidate(1, "target", (11,), (0.0, 0.0, 0.0))
    insertion = ResidueCorrespondence(1, None, _residue(11, "target"), None, "A", "insertion")
    deletion = ResidueCorrespondence(2, _residue(10, "ref"), None, "A", None, "deletion")

    result = match_pocket_candidates((reference,), (target,), (insertion, deletion))
    assert {item.state for item in result.matches} == {
        PocketMatchState.UNMATCHED_REFERENCE,
        PocketMatchState.UNMATCHED_TARGET,
    }
    assert match_pocket_candidates((), (), ()).availability is Availability.NOT_APPLICABLE


@pytest.mark.parametrize(
    "kwargs",
    (
        {"minimum_lining_jaccard": "bad"},
        {"minimum_lining_jaccard": -0.1},
        {"minimum_lining_jaccard": 1.1},
        {"maximum_centroid_distance_angstrom": float("nan")},
        {"ambiguity_margin": -0.1},
        {"lining_weight": 0.0, "centroid_weight": 0.0},
        {"maximum_candidate_count": True},
        {"maximum_candidate_count": 0},
    ),
)
def test_matching_settings_validate_bounds(kwargs: dict[str, object]) -> None:
    with pytest.raises(ValueError):
        PocketMatchingSettings(**kwargs)


def test_matching_contracts_validate_inputs_and_keep_serialization_immutable() -> None:
    settings = PocketMatchingSettings(
        min_lining_jaccard=0.0,
        max_centroid_distance_angstrom=0.0,
        lining_overlap_weight=1.0,
        centroid_distance_weight=0.0,
    )
    candidate = _candidate(1, "ref", (10,), (0.0, 0.0, 0.0))
    target = _candidate(2, "target", (10,), (0.0, 0.0, 0.0))
    match = match_pocket_candidates(
        (candidate,), (target,), {_residue(10, "ref"): _residue(10, "target")}, settings=settings
    ).matches[0]
    assert match.status is match.kind is PocketMatchState.MATCHED
    assert match.lining_residue_jaccard == match.lining_overlap_fraction == 1.0
    assert match.transformed_centroid_distance_angstrom == 0.0
    assert match.to_json()["status"] == "matched"
    with pytest.raises(FrozenInstanceError):
        match.score = 2.0  # type: ignore[misc]
    with pytest.raises(TypeError):
        match_pocket_candidates((object(),), (), ())  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        match_pocket_candidates((candidate,), (), (object(),))  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        match_pocket_candidates((candidate,), (), {object(): _residue(10, "target")})  # type: ignore[dict-item]
    duplicate = (
        ResidueCorrespondence(0, _residue(10, "ref"), _residue(10, "target"), "A", "A", "conserved"),
        ResidueCorrespondence(1, _residue(10, "ref"), _residue(11, "target"), "A", "A", "conserved"),
    )
    with pytest.raises(ValueError):
        match_pocket_candidates((candidate,), (), duplicate)
    assert Diagnostic("match.info", DiagnosticSeverity.INFO, "ok").to_json()["code"] == "match.info"


def test_matching_uses_bounded_greedy_path_for_large_collections() -> None:
    references = tuple(_candidate(index, "ref", (10,), (float(index), 0.0, 0.0)) for index in range(9))
    targets = tuple(_candidate(index, "target", (10,), (float(index), 0.0, 0.0)) for index in range(9))
    correspondences = _correspondences((10,))
    result = match_pocket_candidates(references, targets, correspondences)
    assert len(result.matches) == 9
    assert all(item.state in {PocketMatchState.MATCHED, PocketMatchState.AMBIGUOUS} for item in result.matches)


def test_matching_reports_candidate_resource_limit() -> None:
    candidate = _candidate(1, "ref", (10,), (0.0, 0.0, 0.0))
    result = match_pocket_candidates(
        (candidate, candidate), (), settings=PocketMatchingSettings(maximum_candidate_count=1)
    )
    assert result.availability is Availability.NUMERICAL_FAILURE


def test_matching_checks_sequence_length_before_materializing_oversized_candidates() -> None:
    class OversizedSequence(Sequence[PocketCandidate]):
        def __len__(self) -> int:
            return 2

        def __getitem__(self, index: int) -> PocketCandidate:
            raise AssertionError("oversized sequence must not be materialized")

    result = match_pocket_candidates(
        OversizedSequence(),
        (),
        settings=PocketMatchingSettings(maximum_candidate_count=1),
    )

    assert result.availability is Availability.NUMERICAL_FAILURE


def test_matching_match_and_result_contract_validation_paths() -> None:
    candidate = _candidate(1, "ref", (10,), (0.0, 0.0, 0.0))
    target = _candidate(1, "target", (10,), (0.0, 0.0, 0.0))
    valid = match_pocket_candidates((candidate,), (target,), _correspondences((10,)))
    match = valid.matches[0]
    assert (match.reference, match.target) == (candidate, target)
    assert valid.status is valid.availability
    assert valid.candidate_matches[0] == valid[0]
    assert len(valid) == len(tuple(valid)) == 1
    settings = PocketMatchingSettings()
    assert settings.min_lining_jaccard == settings.minimum_lining_jaccard
    assert settings.max_centroid_distance_angstrom == settings.maximum_centroid_distance_angstrom
    assert settings.centroid_distance_weight == settings.centroid_weight
    assert settings.lining_overlap_weight == settings.lining_weight
    for kwargs in (
        {"reference_candidate": object(), "target_candidate": target, "state": PocketMatchState.MATCHED},
        {"reference_candidate": candidate, "target_candidate": object(), "state": PocketMatchState.MATCHED},
        {"reference_candidate": candidate, "target_candidate": None, "state": "bad"},
        {"reference_candidate": None, "target_candidate": target, "state": PocketMatchState.UNMATCHED_REFERENCE},
        {"reference_candidate": candidate, "target_candidate": None, "state": PocketMatchState.UNMATCHED_TARGET},
        {"reference_candidate": None, "target_candidate": None, "state": PocketMatchState.MATCHED},
        {
            "reference_candidate": candidate,
            "target_candidate": target,
            "state": PocketMatchState.MATCHED,
            "lining_jaccard": 2.0,
        },
        {
            "reference_candidate": candidate,
            "target_candidate": target,
            "state": PocketMatchState.MATCHED,
            "diagnostics": (object(),),
        },
        {
            "reference_candidate": candidate,
            "target_candidate": target,
            "state": PocketMatchState.MATCHED,
            "units": {"distance": ""},
        },
    ):
        with pytest.raises((TypeError, ValueError)):
            from structlens.core.pockets.matching import PocketMatch

            PocketMatch(**kwargs)  # type: ignore[arg-type]
    from structlens.core.pockets.matching import PocketMatchingResult

    for kwargs in (
        {"availability": "bad"},
        {"availability": Availability.AVAILABLE, "matches": (object(),)},
        {"availability": Availability.AVAILABLE, "diagnostics": (object(),)},
        {"availability": Availability.AVAILABLE, "settings": object()},
        {"availability": Availability.AVAILABLE, "transform": object()},
    ):
        with pytest.raises((TypeError, ValueError)):
            PocketMatchingResult(**kwargs)  # type: ignore[arg-type]


def test_matching_assignment_skips_used_targets_and_marks_greedy_conflicts() -> None:
    reference_a = _candidate(1, "ref", (10,), (0.0, 0.0, 0.0))
    reference_b = _candidate(2, "ref", (10,), (0.0, 0.0, 0.0))
    target = _candidate(1, "target", (10,), (0.0, 0.0, 0.0))
    rows = _correspondences((10,))
    result = match_pocket_candidates((reference_a, reference_b), (target,), rows)
    assert {item.state for item in result.matches} == {
        PocketMatchState.AMBIGUOUS,
        PocketMatchState.UNMATCHED_REFERENCE,
    }
    refs = tuple(_candidate(index, "ref", (10,), (float(index), 0.0, 0.0)) for index in range(9))
    greedy = match_pocket_candidates(refs, (target,), rows)
    assert sum(item.state in {PocketMatchState.MATCHED, PocketMatchState.AMBIGUOUS} for item in greedy.matches) == 1


def test_matching_rejects_duplicate_target_and_invalid_runtime_arguments() -> None:
    candidate = _candidate(1, "ref", (10,), (0.0, 0.0, 0.0))
    target = _candidate(1, "target", (10,), (0.0, 0.0, 0.0))
    duplicate = (
        ResidueCorrespondence(0, _residue(10, "ref"), _residue(10, "target"), "A", "A", "conserved"),
        ResidueCorrespondence(1, _residue(11, "ref"), _residue(10, "target"), "A", "A", "conserved"),
    )
    with pytest.raises(ValueError):
        match_pocket_candidates((candidate,), (target,), duplicate)
    with pytest.raises(TypeError):
        match_pocket_candidates((candidate,), (target,), settings=object())  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        match_pocket_candidates((candidate,), (target,), transform=object())  # type: ignore[arg-type]


def test_matching_uses_global_optimum_for_nine_by_nine_assignment(monkeypatch: pytest.MonkeyPatch) -> None:
    references = tuple(_candidate(index, "ref", (10,), (float(index), 0.0, 0.0)) for index in range(9))
    targets = tuple(_candidate(index, "target", (10,), (float(index), 0.0, 0.0)) for index in range(9))

    def weighted_edge(reference: PocketCandidate, target: PocketCandidate, *args: object) -> tuple[float, float, float]:
        row = int(reference.centroid_xyz[0])
        column = int(target.centroid_xyz[0])
        matrix = {(0, 0): 0.90, (0, 1): 0.80, (1, 0): 0.85, (1, 1): 0.10}
        score = matrix.get((row, column), 0.50 if row == column else 0.10)
        return score, 0.0, score

    monkeypatch.setattr("structlens.core.pockets.matching._edge", weighted_edge)
    result = match_pocket_candidates(references, targets, _correspondences((10,)))
    paired = {
        (int(item.reference_candidate.centroid_xyz[0]), int(item.target_candidate.centroid_xyz[0]))
        for item in result.matches
        if item.reference_candidate is not None and item.target_candidate is not None
    }
    assert (0, 1) in paired
    assert (1, 0) in paired
    assert result.optimal_score == pytest.approx(5.15)


def test_matching_ambiguity_compares_global_totals_and_is_deterministic(monkeypatch: pytest.MonkeyPatch) -> None:
    references = tuple(_candidate(index, "ref", (10,), (float(index), 0.0, 0.0)) for index in range(2))
    targets = tuple(_candidate(index, "target", (10,), (float(index), 0.0, 0.0)) for index in range(2))

    scores = [[0.90, 0.89], [0.89, 0.10]]

    def weighted_edge(reference: PocketCandidate, target: PocketCandidate, *args: object) -> tuple[float, float, float]:
        row = int(reference.centroid_xyz[0])
        column = int(target.centroid_xyz[0])
        score = scores[row][column]
        return score, 0.0, score

    monkeypatch.setattr("structlens.core.pockets.matching._edge", weighted_edge)
    result = match_pocket_candidates(
        references,
        targets,
        _correspondences((10,)),
        settings=PocketMatchingSettings(ambiguity_margin=0.02),
    )
    assert result.optimal_score == pytest.approx(1.78)
    assert result.second_best_score == pytest.approx(1.00)
    assert all(item.state is PocketMatchState.MATCHED for item in result.matches)

    scores[:] = [[0.50, 0.50], [0.50, 0.50]]
    tie = match_pocket_candidates(
        references,
        targets,
        _correspondences((10,)),
        settings=PocketMatchingSettings(ambiguity_margin=0.0),
    )
    assert all(item.state is PocketMatchState.AMBIGUOUS for item in tie.matches)
    assert (
        tie.matches
        == match_pocket_candidates(
            tuple(reversed(references)),
            tuple(reversed(targets)),
            _correspondences((10,)),
            settings=PocketMatchingSettings(ambiguity_margin=0.0),
        ).matches
    )


def test_matching_rejects_duplicate_references_even_when_the_target_is_equal() -> None:
    candidate = _candidate(1, "ref", (10,), (0.0, 0.0, 0.0))
    target = _candidate(1, "target", (10,), (0.0, 0.0, 0.0))
    duplicate = (
        ResidueCorrespondence(0, _residue(10, "ref"), _residue(10, "target"), "A", "A", "conserved"),
        ResidueCorrespondence(1, _residue(10, "ref"), _residue(10, "target"), "A", "A", "conserved"),
    )
    with pytest.raises(ValueError, match="duplicate reference"):
        match_pocket_candidates((candidate,), (target,), duplicate)


def test_matching_global_alternative_includes_explicit_unmatched_assignment() -> None:
    reference = _candidate(1, "ref", (10,), (0.0, 0.0, 0.0))
    target = _candidate(1, "target", (10,), (0.0, 0.0, 0.0))

    result = match_pocket_candidates(
        (reference,),
        (target,),
        _correspondences((10,)),
        settings=PocketMatchingSettings(ambiguity_margin=1.0),
    )

    assert result.matches[0].state is PocketMatchState.AMBIGUOUS
    assert result.matches[0].alternative_assignment_score == 0.0


def test_matching_does_not_trade_a_real_score_for_tie_break_epsilon(monkeypatch: pytest.MonkeyPatch) -> None:
    references = tuple(_candidate(index, "ref", (10,), (float(index), 0.0, 0.0)) for index in range(2))
    targets = (_candidate(1, "target", (10,), (0.0, 0.0, 0.0)),)

    def weighted_edge(reference: PocketCandidate, target: PocketCandidate, *args: object) -> tuple[float, float, float]:
        return (1.0 - 1.0e-13, 0.0, 1.0 - 1.0e-13) if reference.centroid_xyz[0] == 0.0 else (1.0, 0.0, 1.0)

    monkeypatch.setattr("structlens.core.pockets.matching._edge", weighted_edge)
    result = match_pocket_candidates(references, targets, _correspondences((10,)))

    matched = next(item for item in result.matches if item.target_candidate is not None)
    assert matched.reference_candidate is not None
    assert matched.reference_candidate.centroid_xyz[0] == 1.0
