from __future__ import annotations

from structlens.application.project_state import ProjectState
from structlens.core.metrics.sequence_metrics import SequenceAlignmentMetrics
from structlens.core.models import (
    AlignmentMode,
    AnalysisSettings,
    ComparisonMode,
    ReferenceVsManyAnalysis,
    ResidueId,
    TargetAnalysis,
)


def test_project_state_json_round_trip_preserves_settings_and_keys() -> None:
    state = ProjectState(
        reference_source="reference.pdb",
        target_sources=("target.pdb",),
        settings=AnalysisSettings(alignment_mode=AlignmentMode.SEQUENCE),
        key_residues=(ResidueId("ref", "1", "A", "130", "A", "SER"),),
        source_objects={"reference": "reference_obj", "target": "target_obj"},
    )

    restored = ProjectState.from_json(state.to_json())

    assert restored == state
    assert restored.settings.alignment_mode is AlignmentMode.SEQUENCE
    assert restored.key_residues[0].insertion_code == "A"
    assert restored.source_objects["reference"] == "reference_obj"


def test_project_state_can_record_sha256_source_hash(tmp_path) -> None:
    source = tmp_path / "reference.pdb"
    source.write_text("ATOM\n", encoding="utf-8")
    state = ProjectState(reference_source=str(source)).with_source_hashes()

    assert len(state.source_hashes[str(source)]) == 64


def test_project_state_round_trips_reference_vs_many() -> None:
    target = TargetAnalysis(
        target_id="target",
        correspondence=(),
        mutations=(),
        sequence_metrics=SequenceAlignmentMetrics(0.8, 0.9, 1.0, 1.0, 1.0, 10),
    )
    state = ProjectState(
        comparison_mode=ComparisonMode.REFERENCE_VS_MANY,
        reference_vs_many=ReferenceVsManyAnalysis("reference", {"target": target}),
    )
    restored = ProjectState.from_json(state.to_json())
    assert restored.comparison_mode is ComparisonMode.REFERENCE_VS_MANY
    assert restored.reference_vs_many is not None
    assert restored.reference_vs_many.targets["target"].sequence_metrics.identity == 0.8
