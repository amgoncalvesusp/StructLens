from __future__ import annotations

from pathlib import Path

import pytest

import structlens.application.project_state as project_state_module
from structlens.application.project_state import ProjectState
from structlens.core.errors import ProjectSchemaError
from structlens.core.metrics.sequence_metrics import SequenceAlignmentMetrics
from structlens.core.models import (
    AlignmentMode,
    AnalysisResult,
    AnalysisSettings,
    ComparisonMode,
    ReferenceVsManyAnalysis,
    ResidueId,
    TargetAnalysis,
)
from structlens.core.parsing import AltlocPolicy, AssemblyScope, InputSelection, StructureFormat
from structlens.core.provenance import MethodProvenance


def test_project_state_round_trips_full_typed_provenance() -> None:
    provenance = MethodProvenance(
        "structlens.compare",
        "0.4.0",
        {"grid_spacing": 0.5, "enabled": True},
        {"grid_spacing": "angstrom"},
        {"scipy": "1.17.1"},
        {"target": "a" * 64},
        "as_file",
    )
    analysis = AnalysisResult(
        reference_id="ref",
        target_id="target",
        correspondences=(),
        mutations=(),
        sequence_identity=1.0,
        sequence_coverage=1.0,
        alignment_decision="sequence",
        method_provenance=provenance,
    )
    restored = ProjectState.from_json(ProjectState(analysis_results=(analysis,)).to_json())

    assert restored.analysis_results[0].method_provenance == provenance
    assert restored.analysis_results[0].method_provenance is not None
    assert restored.analysis_results[0].method_provenance.parameters["grid_spacing"] == 0.5


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


def test_project_state_can_record_sha256_source_hash(tmp_path: Path) -> None:
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


def test_project_state_round_trips_reference_vs_many_typed_provenance() -> None:
    provenance = MethodProvenance(
        "structlens.compare",
        "0.4.0",
        {"grid_spacing": 0.5},
        {"grid_spacing": "angstrom"},
        input_hashes={"target": "b" * 64},
    )
    target = TargetAnalysis(
        target_id="target",
        correspondence=(),
        mutations=(),
        sequence_metrics=SequenceAlignmentMetrics(0.8, 0.9, 1.0, 1.0, 1.0, 10),
        provenance={"backend": "legacy"},
        method_provenance=provenance,
    )
    state = ProjectState(
        comparison_mode=ComparisonMode.REFERENCE_VS_MANY,
        reference_vs_many=ReferenceVsManyAnalysis("reference", {"target": target}),
    )

    restored = ProjectState.from_json(state.to_json())

    assert restored.reference_vs_many is not None
    restored_target = restored.reference_vs_many.targets["target"]
    assert restored_target.provenance == {"backend": "legacy"}
    assert restored_target.method_provenance == provenance


def test_legacy_project_migrates_to_v04_with_selection_and_pocket_settings() -> None:
    selection = InputSelection(
        "a" * 64,
        "reference.cif",
        StructureFormat.MMCIF,
        "2",
        author_chain_ids=("A",),
        label_chain_ids=("L",),
        altloc_policy=AltlocPolicy.FIRST,
        assembly_scope=AssemblyScope.BIOLOGICAL_ASSEMBLY,
        path="/private/reference.cif",
    )
    legacy = ProjectState(
        schema_version="0.3",
        source_hashes={"/private/reference.cif": "b" * 64},
        msa_settings={"algorithm": "muscle5"},
    )
    migrated = legacy.migrate()
    migrated = ProjectState(
        schema_version="4.0",
        source_hashes=dict(migrated.source_hashes),
        msa_settings=dict(migrated.msa_settings),
        reference_selection=selection,
        target_selection=selection,
        pocket_settings={"minimum_cluster_size": 2},
        report_references=("report.json",),
    )
    restored = ProjectState.from_json(migrated.to_json())
    assert restored.schema_version == "4.0"
    assert restored.reference_selection == selection
    assert restored.reference_selection is not None
    assert restored.reference_selection.selection_id == selection.selection_id
    assert restored.source_hashes["/private/reference.cif"] == "b" * 64
    assert restored.pocket_settings["minimum_cluster_size"] == 2


def test_project_state_rejects_future_schema() -> None:
    with pytest.raises(ProjectSchemaError, match="Unsupported"):
        ProjectState.from_dict({"schema_version": "5.0"})


def test_project_state_loads_handwritten_v02_as_v04_without_dropping_hashes() -> None:
    restored = ProjectState.from_dict(
        {
            "schema_version": "0.2",
            "reference_source": "reference.pdb",
            "source_hashes": {"reference.pdb": "a" * 64},
            "settings": {"alignment_mode": "sequence"},
        }
    )
    assert restored.schema_version == "4.0"
    assert restored.settings.alignment_mode is AlignmentMode.SEQUENCE
    assert restored.source_hashes["reference.pdb"] == "a" * 64


def test_project_migration_rejects_nonfinite_settings_and_oversized_collections() -> None:
    with pytest.raises(ProjectSchemaError, match="finite"):
        ProjectState.from_dict({"schema_version": "0.2", "settings": {"minimum_sequence_identity": float("nan")}})

    with pytest.raises(ProjectSchemaError, match="too many|limit"):
        ProjectState.from_dict({"schema_version": "0.3", "target_sources": ["target.pdb"] * 100_001})


def test_project_load_bounds_file_size(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    path = tmp_path / "project.json"
    path.write_text("{}", encoding="utf-8")
    monkeypatch.setattr(project_state_module, "_MAX_PROJECT_BYTES", 1)
    with pytest.raises(ProjectSchemaError, match="size"):
        ProjectState.load(path)
