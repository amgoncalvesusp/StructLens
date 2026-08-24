from structlens.core.metrics.sequence_metrics import SequenceAlignmentMetrics
from structlens.core.models import (
    ComparisonMode,
    MultiStructurePosition,
    PairwiseMatrix,
    ReferenceVsManyAnalysis,
    TargetAnalysis,
)


def _target(target_id: str) -> TargetAnalysis:
    return TargetAnalysis(
        target_id=target_id,
        correspondence=(),
        mutations=(),
        sequence_metrics=SequenceAlignmentMetrics(1.0, 1.0, 1.0, 1.0, 1.0, 0),
    )


def test_reference_vs_many_preserves_independent_target_results() -> None:
    analysis = ReferenceVsManyAnalysis("ref", {"b": _target("b"), "a": _target("a")})
    assert analysis.comparison_mode is ComparisonMode.REFERENCE_VS_MANY
    assert analysis.target_ids == ("b", "a")
    assert analysis.targets["a"].target_id == "a"


def test_pairwise_matrix_mirrors_one_stored_value() -> None:
    matrix = PairwiseMatrix.from_pairs(
        "sequence_identity", ("A", "B", "C"), {("B", "A"): 0.75}, unit="fraction"
    )
    assert matrix.value("A", "B") == 0.75
    assert matrix.value("B", "A") == 0.75
    assert matrix.value("A", "C") is None
    assert tuple(matrix.values) == (("A", "B"),)


def test_multi_structure_position_validates_coverage() -> None:
    position = MultiStructurePosition(0, None, {"A": None}, 0.5, 1.2)
    assert position.coverage == 0.5
