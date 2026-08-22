from __future__ import annotations

from structlens.application.project_state import ProjectState
from structlens.core.models import AlignmentMode, AnalysisSettings, ResidueId


def test_project_state_json_round_trip_preserves_settings_and_keys() -> None:
    state = ProjectState(
        reference_source="reference.pdb",
        target_sources=("target.pdb",),
        settings=AnalysisSettings(alignment_mode=AlignmentMode.SEQUENCE),
        key_residues=(ResidueId("ref", "1", "A", "130", "A", "SER"),),
    )

    restored = ProjectState.from_json(state.to_json())

    assert restored == state
    assert restored.settings.alignment_mode is AlignmentMode.SEQUENCE
    assert restored.key_residues[0].insertion_code == "A"


def test_project_state_can_record_sha256_source_hash(tmp_path) -> None:
    source = tmp_path / "reference.pdb"
    source.write_text("ATOM\n", encoding="utf-8")
    state = ProjectState(reference_source=str(source)).with_source_hashes()

    assert len(state.source_hashes[str(source)]) == 64
