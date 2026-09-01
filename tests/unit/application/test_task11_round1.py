from __future__ import annotations

import gzip
import hashlib
import importlib
import json
from copy import deepcopy
from dataclasses import replace
from pathlib import Path
from typing import Any, cast

import pytest

import structlens.application.report_snapshot_io as _snapshot_io_module
import structlens.core.reports.safe_io as _safe_io_module
from structlens.application.report_serialization import (
    deserialize_report,
    report_to_dict,
    serialize_report,
    store_snapshot,
)
from structlens.application.report_snapshot_io import verify_report_inputs, verify_source_path
from structlens.core.evidence import (
    Availability,
    Diagnostic,
    DiagnosticSeverity,
    PocketEvidenceChannel,
    PocketEvidenceSections,
)
from structlens.core.models import ResidueId
from structlens.core.parsing import InputSelection, SourceSnapshot, StructureFormat
from structlens.core.pockets.comparison_models import PocketComparison
from structlens.core.pockets.matching import PocketMatch, PocketMatchingResult, PocketMatchState
from structlens.core.pockets.models import (
    AlphaSphere,
    AlphaSphereDetectionResult,
    PocketCandidate,
    PocketDetectionSettings,
)
from structlens.core.pockets.volume_models import (
    PocketVolumeComparison,
    PocketVolumeResult,
    PocketVolumeSensitivity,
    PocketVolumeSettings,
)
from structlens.core.provenance import MethodProvenance
from structlens.core.quality import StructureQualityReport
from structlens.core.reports import (
    AnalysisReport,
    InputQualityBundle,
    PocketDetectionSnapshot,
    PocketLiningSnapshot,
    PocketReportSnapshot,
    PocketVolumeSnapshot,
    SectionAvailability,
)
from structlens.core.reports.safe_io import atomic_write_bytes, read_bounded_bytes
from structlens.core.reports.schema import load_analysis_report_schema, validate_analysis_report_payload

safe_io: Any = _safe_io_module
snapshot_io: Any = _snapshot_io_module


def _snapshot(name: str, content: bytes, *, gzip: bool = False) -> SourceSnapshot:
    digest = hashlib.sha256(content).hexdigest()
    return SourceSnapshot(content, digest, digest, digest, name, "pdb", gzip)


def _report(
    reference: SourceSnapshot, target: SourceSnapshot, *, pockets: PocketReportSnapshot | None = None
) -> AnalysisReport:
    ref = InputSelection(reference.content_id, reference.display_name, StructureFormat.PDB, "1")
    tgt = InputSelection(target.content_id, target.display_name, StructureFormat.PDB, "1")
    return AnalysisReport(
        ref,
        tgt,
        InputQualityBundle(
            StructureQualityReport(Availability.NOT_APPLICABLE), StructureQualityReport(Availability.NOT_APPLICABLE)
        ),
        availability=SectionAvailability(pockets=Availability.AVAILABLE if pockets else Availability.NOT_APPLICABLE),
        provenance=MethodProvenance(
            "test.report",
            "1",
            input_hashes={
                "reference_content": ref.content_id,
                "target_content": tgt.content_id,
                "reference_raw": reference.raw_sha256,
                "target_raw": target.raw_sha256,
            },
        ),
        pockets=pockets,
    )


def _pocket_report(reference: SourceSnapshot, target: SourceSnapshot | None = None) -> PocketReportSnapshot:
    residue = ResidueId("reference", "1", "A", "10", None, "ALA")
    sphere = AlphaSphere((1.0, 2.0, 3.0), 3.0, ("a", "b", "c", "d"), (residue,), ("a", "b", "c", "d"))
    candidate = PocketCandidate(
        (sphere,),
        reference.content_id,
        InputSelection(reference.content_id, reference.display_name, StructureFormat.PDB, "1").selection_id,
    )
    detection = PocketDetectionSnapshot(
        "reference",
        Availability.AVAILABLE,
        candidates=(candidate,),
        detection=AlphaSphereDetectionResult((sphere,), ()),
        settings=PocketDetectionSettings(),
        counts={"candidate_count": 1},
        provenance=MethodProvenance("pocket.detect", "1", input_hashes={"source": reference.content_id}),
    )
    detections = (detection,)
    target_candidate: PocketCandidate | None = None
    if target is not None:
        target_residue = ResidueId("target", "1", "A", "10", None, "SER")
        target_sphere = AlphaSphere((4.0, 5.0, 6.0), 3.0, ("a", "b", "c", "d"), (target_residue,), ("a", "b", "c", "d"))
        target_candidate = PocketCandidate(
            (target_sphere,),
            target.content_id,
            InputSelection(target.content_id, target.display_name, StructureFormat.PDB, "1").selection_id,
        )
        detections += (
            PocketDetectionSnapshot(
                "target",
                Availability.AVAILABLE,
                candidates=(target_candidate,),
                detection=AlphaSphereDetectionResult((target_sphere,), ()),
                settings=PocketDetectionSettings(),
                counts={"candidate_count": 1},
                provenance=MethodProvenance("pocket.detect", "1", input_hashes={"source": target.content_id}),
            ),
        )
    volume = PocketVolumeResult(
        Availability.AVAILABLE,
        coarse_voxel_count=10,
        fine_voxel_count=20,
        coarse_volume_angstrom3=10.0,
        fine_volume_angstrom3=8.0,
        sensitivity=PocketVolumeSensitivity(2.0, 0.25),
        candidate_id=candidate.candidate_id,
        settings=PocketVolumeSettings(),
        provenance=MethodProvenance("pocket.volume", "1", input_hashes={"source": reference.content_id}),
    )
    match = (
        PocketMatch(
            candidate,
            target_candidate,
            PocketMatchState.MATCHED,
            lining_jaccard=0.8,
            centroid_distance_angstrom=0.5,
            score=0.75,
            alternative_score=0.2,
            ambiguity_margin=0.55,
        )
        if target_candidate is not None
        else None
    )
    comparison = (
        PocketComparison(
            Availability.AVAILABLE,
            match,
            volume=PocketVolumeComparison(
                Availability.AVAILABLE,
                delta_angstrom3=1.0,
                relative_delta_fraction=0.125,
                reference_volume_angstrom3=8.0,
                target_volume_angstrom3=9.0,
                reference_sensitivity=PocketVolumeSensitivity(2.0, 0.25),
                target_sensitivity=PocketVolumeSensitivity(2.0, 0.22),
            ),
            surface_delta_angstrom2=0.5,
            relative_surface_delta_fraction=0.1,
            reference_surface_area_angstrom2=5.0,
            target_surface_area_angstrom2=5.5,
            reference_surface_method="surface-v1",
            target_surface_method="surface-v1",
            reference_surface_provenance={"artifact": "ref"},
            target_surface_provenance={"artifact": "tgt"},
            reference_surface_units={"surface_area": "angstrom^2"},
            target_surface_units={"surface_area": "angstrom^2"},
            lining_residue_conserved=(residue,),
            local_displacement_angstrom=1.0,
            ca_displacement_angstrom=0.8,
            sidechain_displacement_angstrom=1.2,
            local_displacements=((residue, 1.0),),
        )
        if match is not None
        else None
    )
    unavailable_comparison = PocketComparison(Availability.NOT_APPLICABLE, match) if match is not None else None
    return PocketReportSnapshot(
        detections=detections,
        volumes=(
            PocketVolumeSnapshot("reference", volume),
            PocketVolumeSnapshot("target", None, Availability.NOT_APPLICABLE),
        ),
        matching=(
            PocketMatchingResult(Availability.AVAILABLE, matches=(match,), optimal_score=0.75, second_best_score=0.2)
            if match is not None
            else None
        ),
        comparisons=tuple(item for item in (comparison, unavailable_comparison) if item is not None),
        lining_residues=(
            PocketLiningSnapshot("reference", candidate.candidate_id, (residue,), units={"residue": "count"}),
        ),
        concordance=PocketEvidenceSections(
            **{
                name: PocketEvidenceChannel(Availability.NOT_APPLICABLE)
                for name in (
                    "geometry",
                    "ligand_support",
                    "parameter_persistence",
                    "volume_sensitivity",
                    "match_ambiguity",
                    "qc",
                    "interactions",
                )
            }
        ),
        diagnostics=(Diagnostic("pocket.test", DiagnosticSeverity.WARNING, " =formula"),),
        units={"volume": "angstrom^3"},
    )


def test_roundtrip_rejects_unknown_nested_fields() -> None:
    reference = _snapshot("reference.pdb", b"reference")
    target = _snapshot("target.pdb", b"target")
    payload = report_to_dict(_report(reference, target))
    payload["input_quality"]["reference"]["unexpected"] = "loss"
    with pytest.raises(ValueError, match="canonical|unknown|schema"):
        deserialize_report(payload, snapshots=(reference, target))


def test_roundtrip_rejects_unknown_nested_fields_without_jsonschema(monkeypatch: pytest.MonkeyPatch) -> None:
    reference = _snapshot("reference.pdb", b"reference")
    target = _snapshot("target.pdb", b"target")
    payload = report_to_dict(_report(reference, target))
    payload["input_quality"]["target"]["unexpected"] = True
    real_import = importlib.import_module

    def missing_jsonschema(name: str) -> object:
        if name == "jsonschema":
            raise ImportError("optional dependency absent")
        return real_import(name)

    monkeypatch.setattr("structlens.core.reports.schema.importlib.import_module", missing_jsonschema)
    with pytest.raises(ValueError, match="canonical|unknown|schema"):
        deserialize_report(payload, snapshots=(reference, target))


def test_canonical_serialization_requires_complete_raw_evidence_and_gzip_identity() -> None:
    reference = _snapshot("reference.pdb", b"reference")
    target = _snapshot("target.pdb", b"target")
    report = _report(reference, target)
    with pytest.raises(ValueError, match="evidence|raw|snapshot"):
        serialize_report(report)
    alternate = SourceSnapshot(
        b"reference",
        "1" * 64,
        reference.content_id,
        reference.content_id,
        "reference.pdb.gz",
        "pdb",
        True,
    )
    with pytest.raises(ValueError, match="raw|snapshot|source"):
        serialize_report(report, snapshots=(alternate, target))


def test_snapshot_store_is_idempotent_and_refuses_conflicts_or_links(tmp_path: Path) -> None:
    snapshot = _snapshot("reference.pdb", b"reference")
    evidence = tmp_path / "evidence"
    store_snapshot(snapshot, evidence)
    store_snapshot(snapshot, evidence)
    target = evidence / f"{snapshot.content_id}.snapshot"
    target.write_bytes(b"tampered")
    with pytest.raises(ValueError, match="hash|tamper|snapshot"):
        store_snapshot(snapshot, evidence)
    target.unlink()
    target.write_bytes(snapshot.decompressed_bytes)
    hardlink = evidence / "hardlink.snapshot"
    hardlink.hardlink_to(target)
    with pytest.raises(ValueError, match="link|snapshot"):
        store_snapshot(snapshot, evidence)


def test_safe_artifact_io_is_bounded_atomic_and_link_resistant(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    path = tmp_path / "artifact.bin"
    with pytest.raises(TypeError, match="bytes"):
        atomic_write_bytes(path, "not bytes", max_bytes=100)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="size|supported"):
        atomic_write_bytes(path, b"123", max_bytes=2)
    atomic_write_bytes(path, b"same", max_bytes=100, idempotent=True)
    atomic_write_bytes(path, b"same", max_bytes=100, idempotent=True)
    with pytest.raises(FileExistsError, match="conflicting"):
        atomic_write_bytes(path, b"different", max_bytes=100, idempotent=True)
    hardlink = tmp_path / "artifact-hardlink.bin"
    hardlink.hardlink_to(path)
    with pytest.raises(ValueError, match="link"):
        read_bounded_bytes(hardlink, max_bytes=100)
    with pytest.raises(ValueError, match="link"):
        atomic_write_bytes(hardlink, b"new", max_bytes=100, idempotent=True)

    original_replace = safe_io.os.replace

    def fail_replace(*_args: object, **_kwargs: object) -> None:
        raise OSError("injected replace failure")

    monkeypatch.setattr(safe_io.os, "replace", fail_replace)
    with pytest.raises(OSError, match="replace"):
        atomic_write_bytes(tmp_path / "partial.bin", b"partial", max_bytes=100)
    assert not tuple(tmp_path.glob(".partial.bin.*"))
    monkeypatch.setattr(safe_io.os, "replace", original_replace)


def test_safe_io_handles_missing_growth_text_and_cleanup_edges(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    missing = tmp_path / "missing.bin"
    with pytest.raises(FileNotFoundError):
        read_bounded_bytes(missing, max_bytes=100)
    growing = tmp_path / "growing.bin"
    growing.write_bytes(b"ok")
    with pytest.raises(ValueError, match="size|supported"):
        read_bounded_bytes(growing, max_bytes=1)
    monkeypatch.setattr(safe_io.os, "read", lambda _descriptor, _size: b"x" * 3)
    with pytest.raises(ValueError, match="size|supported"):
        read_bounded_bytes(growing, max_bytes=2)
    monkeypatch.undo()

    overwrite = tmp_path / "overwrite.bin"
    atomic_write_bytes(overwrite, b"old", max_bytes=100)
    with pytest.raises(FileExistsError, match="changed"):
        atomic_write_bytes(overwrite, b"new", max_bytes=100)
    assert overwrite.read_bytes() == b"old"
    with pytest.raises(TypeError, match="text"):
        safe_io.atomic_write_text(tmp_path / "text", b"not text", max_bytes=100)

    original_replace = safe_io.os.replace
    original_unlink = safe_io.os.unlink

    def fail_replace(*_args: object, **_kwargs: object) -> None:
        raise OSError("injected replace failure")

    def missing_unlink(*_args: object, **_kwargs: object) -> None:
        raise FileNotFoundError("already cleaned")

    monkeypatch.setattr(safe_io.os, "replace", fail_replace)
    monkeypatch.setattr(safe_io.os, "unlink", missing_unlink)
    with pytest.raises(OSError, match="replace"):
        atomic_write_bytes(tmp_path / "cleanup.bin", b"cleanup", max_bytes=100)
    monkeypatch.setattr(safe_io.os, "replace", original_replace)
    monkeypatch.setattr(safe_io.os, "unlink", original_unlink)

    appeared = tmp_path / "appeared.bin"
    original_regular = safe_io._regular
    calls = 0

    def appear_between_stat_and_replace(path: Path, *, label: str) -> object:
        nonlocal calls
        calls += 1
        if calls == 2:
            path.write_bytes(b"appeared")
        return original_regular(path, label=label)

    monkeypatch.setattr(safe_io, "_regular", appear_between_stat_and_replace)
    with pytest.raises(FileExistsError, match="appeared"):
        atomic_write_bytes(appeared, b"new", max_bytes=100)
    monkeypatch.setattr(safe_io, "_regular", original_regular)


def test_source_binding_rejects_incomplete_duplicate_or_linked_evidence(tmp_path: Path) -> None:
    reference = _snapshot("reference.pdb", b"reference")
    target = _snapshot("target.pdb", b"target")
    report = _report(reference, target)
    with pytest.raises(ValueError, match="at most"):
        verify_report_inputs(report, (reference, target, reference), (), None)
    with pytest.raises(TypeError, match="SourceSnapshot"):
        verify_report_inputs(report, (reference, object()), (), None)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="multiple|conflicting"):
        verify_report_inputs(
            report,
            (
                reference,
                SourceSnapshot(
                    reference.decompressed_bytes,
                    "1" * 64,
                    reference.content_id,
                    reference.content_id,
                    reference.display_name,
                    reference.logical_format,
                    reference.is_gzip,
                ),
            ),
            (),
            None,
        )
    source = tmp_path / "reference.pdb"
    source.write_bytes(reference.decompressed_bytes)
    hardlink = tmp_path / "reference-hardlink.pdb"
    hardlink.hardlink_to(source)
    with pytest.raises(ValueError, match="link"):
        verify_source_path(hardlink, report.reference_selection)


def test_snapshot_store_and_load_fail_closed_at_all_artifact_boundaries(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    reference = _snapshot("reference.pdb", b"reference")
    target = _snapshot("target.pdb", b"target")
    report = _report(reference, target)
    evidence = tmp_path / "evidence"
    with pytest.raises(ValueError, match="digest"):
        snapshot_io.snapshot_path("invalid", evidence)
    with pytest.raises(TypeError, match="SourceSnapshot"):
        snapshot_io.store_snapshot(object(), evidence)
    monkeypatch.setattr(snapshot_io, "MAX_SNAPSHOT_BYTES", 1)
    with pytest.raises(ValueError, match="size|supported"):
        snapshot_io.store_snapshot(reference, evidence)
    monkeypatch.setattr(snapshot_io, "MAX_SNAPSHOT_BYTES", 100 * 1024 * 1024)

    snapshot_io.store_snapshot(reference, evidence)
    with pytest.raises(ValueError, match="conflicts"):
        snapshot_io.store_snapshot(
            SourceSnapshot(
                reference.decompressed_bytes,
                "1" * 64,
                reference.content_id,
                reference.content_id,
                reference.display_name,
                reference.logical_format,
                reference.is_gzip,
            ),
            evidence,
        )
    with pytest.raises(ValueError, match="missing"):
        snapshot_io.load_snapshot(target.content_id, evidence)

    metadata_path = evidence / f"{reference.content_id}.json"
    original_metadata = metadata_path.read_bytes()
    metadata_path.write_bytes(b"[]")
    with pytest.raises(ValueError, match="object|metadata"):
        snapshot_io.load_snapshot(reference.content_id, evidence)
    metadata_path.write_bytes(original_metadata)
    metadata_path.write_text('{"is_gzip": "no"}', encoding="utf-8")
    with pytest.raises(ValueError, match="metadata"):
        snapshot_io.load_snapshot(reference.content_id, evidence)
    metadata_path.write_bytes(original_metadata)
    metadata = json.loads(original_metadata)
    metadata["content_id"] = "0" * 64
    metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
    with pytest.raises(ValueError, match="content ID|metadata"):
        snapshot_io.load_snapshot(reference.content_id, evidence)
    metadata_path.write_bytes(original_metadata)

    original_write = snapshot_io.atomic_write_bytes
    calls = 0

    def fail_metadata(*args: object, **kwargs: object) -> object:
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("injected metadata write failure")
        return original_write(*args, **kwargs)

    fresh = tmp_path / "fresh-evidence"
    monkeypatch.setattr(snapshot_io, "atomic_write_bytes", fail_metadata)
    with pytest.raises(OSError, match="metadata"):
        snapshot_io.store_snapshot(target, fresh)
    assert not tuple(fresh.glob("*"))
    monkeypatch.setattr(snapshot_io, "atomic_write_bytes", original_write)

    gzip_dir = tmp_path / "gzip-evidence"
    snapshot_io.store_snapshot(reference, gzip_dir)
    with pytest.raises(TypeError, match="typed"):
        snapshot_io.verify_snapshot(object(), report.reference_selection)
    wrong_format = SourceSnapshot(
        reference.decompressed_bytes,
        reference.raw_sha256,
        reference.content_id,
        reference.content_id,
        reference.display_name,
        "mmcif",
        False,
    )
    with pytest.raises(ValueError, match="selection"):
        snapshot_io.verify_snapshot(wrong_format, report.reference_selection)

    gzip_source = tmp_path / "reference.pdb.gz"
    gzip_source.write_bytes(gzip.compress(reference.decompressed_bytes))
    with pytest.raises(ValueError, match="conflicts|source"):
        verify_report_inputs(report, (reference, target), (gzip_source,), None)
    with pytest.raises(ValueError, match="conflicts"):
        verify_report_inputs(report, (wrong_format, target), (), gzip_dir)
    with pytest.raises(ValueError, match="content"):
        verify_report_inputs(report, (_snapshot("other.pdb", b"other"), target), (), None)
    with pytest.raises(ValueError, match="provenance|hash"):
        verify_report_inputs(
            replace(
                report,
                provenance=MethodProvenance(
                    "test.report",
                    "1",
                    input_hashes={
                        "reference_content": reference.content_id,
                        "target_content": target.content_id,
                    },
                ),
            ),
            (reference, target),
            (),
            None,
        )


def test_report_schema_fallback_closes_nested_shapes_and_bounds() -> None:
    reference = _snapshot("reference.pdb", b"reference")
    target = _snapshot("target.pdb", b"target")
    base = report_to_dict(_report(reference, target))
    real_import = importlib.import_module

    def missing_jsonschema(name: str) -> object:
        if name == "jsonschema":
            raise ImportError("optional dependency absent")
        return real_import(name)

    import structlens.core.reports.schema as _schema_module

    schema_module: Any = _schema_module

    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(schema_module.importlib, "import_module", missing_jsonschema)
        invalid_cases: list[tuple[str, object]] = [
            ("comparison_mode", "invalid"),
            ("reference_selection", None),
            ("input_quality", None),
            ("analysis", []),
            ("displacement_vectors", None),
            ("interactions", 1),
        ]
        for field, value in invalid_cases:
            candidate = deepcopy(base)
            candidate[field] = value
            with pytest.raises(ValueError, match="supported|object|array"):
                validate_analysis_report_payload(candidate)
        candidate = deepcopy(base)
        candidate["reference_selection"]["unexpected"] = True
        with pytest.raises(ValueError, match="unknown"):
            validate_analysis_report_payload(candidate)
        candidate = deepcopy(base)
        del candidate["target_selection"]["selection_id"]
        with pytest.raises(ValueError, match="missing"):
            validate_analysis_report_payload(candidate)
        candidate = deepcopy(base)
        candidate["provenance"] = "not an object"
        with pytest.raises(ValueError, match="object"):
            validate_analysis_report_payload(candidate)
        candidate = deepcopy(base)
        candidate["diagnostics"] = ["diagnostic"] * 100_001
        with pytest.raises(ValueError, match="items"):
            validate_analysis_report_payload(candidate)
        candidate = deepcopy(base)
        candidate["schema_version"] = "4.0"
        candidate["provenance"] = {"oversized": "x" * 1_000_001}
        with pytest.raises(ValueError, match="oversized"):
            validate_analysis_report_payload(candidate)
        for field, value in (
            ("content_id", "invalid"),
            ("selection_id", "invalid"),
            ("format", "invalid"),
            ("model_id", ""),
            ("author_chain_ids", "A"),
            ("label_chain_ids", "A"),
            ("chain_locators", "A"),
        ):
            malformed = deepcopy(base)
            malformed["reference_selection"][field] = value
            with pytest.raises(ValueError, match="SHA|supported|non-empty|array"):
                validate_analysis_report_payload(malformed)
        malformed = deepcopy(base)
        malformed["reference_selection"] = {"content_id": base["reference_selection"]["content_id"]}
        with pytest.raises(ValueError, match="missing"):
            validate_analysis_report_payload(malformed)

    class OversizedBytes(bytes):
        def __len__(self) -> int:
            return 100 * 1024 * 1024 + 1

    class OversizedString(str):
        def encode(self, _encoding: str = "utf-8", _errors: str = "strict") -> bytes:
            return OversizedBytes()

    with pytest.raises(ValueError, match="size|supported"):
        validate_analysis_report_payload(OversizedString("{}"))
    with pytest.raises(ValueError, match="size|supported"):
        validate_analysis_report_payload(OversizedBytes())
    with pytest.raises(ValueError, match="unsupported"):
        validate_analysis_report_payload({"schema_version": "4.0", "nested": object()})
    nested: object = "leaf"
    for _ in range(66):
        nested = [nested]
    with pytest.raises(ValueError, match="bounds"):
        validate_analysis_report_payload({"schema_version": "4.0", "nested": nested})

    def fail_files(_package: str) -> object:
        raise OSError("schema resource unavailable")

    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(schema_module.resources, "files", fail_files)
        with pytest.raises(ValueError, match="schema"):
            load_analysis_report_schema()

    class ArrayResource:
        def joinpath(self, _name: str) -> ArrayResource:
            return self

        def read_text(self, *, encoding: str) -> str:
            return "[]"

    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(schema_module.resources, "files", lambda _package: ArrayResource())
        with pytest.raises(ValueError, match="object"):
            load_analysis_report_schema()


def test_typed_pockets_survive_canonical_roundtrip() -> None:
    reference = _snapshot("reference.pdb", b"reference")
    target = _snapshot("target.pdb", b"target")
    report = _report(reference, target, pockets=_pocket_report(reference, target))
    restored = deserialize_report(
        serialize_report(report, snapshots=(reference, target)), snapshots=(reference, target)
    )
    assert restored.pockets is not None
    assert isinstance(restored.pockets.detections[0], PocketDetectionSnapshot)
    expected_pockets = report.pockets
    assert expected_pockets is not None
    assert (
        restored.pockets.detections[0].candidates[0].candidate_id
        == expected_pockets.detections[0].candidates[0].candidate_id
    )
    assert restored.pockets.volumes[0].result is not None


def test_pocket_wire_reconstruction_rejects_lossy_or_invalid_nested_values() -> None:
    reference = _snapshot("reference.pdb", b"reference")
    target = _snapshot("target.pdb", b"target")
    payload = report_to_dict(_report(reference, target, pockets=_pocket_report(reference, target)))
    payload["pockets"]["detections"][0]["candidates"][0]["candidate_id"] = "0" * 64
    with pytest.raises(ValueError, match="candidate_id|canonical"):
        deserialize_report(payload, snapshots=(reference, target))

    payload = report_to_dict(_report(reference, target, pockets=_pocket_report(reference, target)))
    payload["pockets"]["volumes"][0]["result"]["grid"]["shape"] = [1, 2]
    with pytest.raises(ValueError, match="shape"):
        deserialize_report(payload, snapshots=(reference, target))

    payload = report_to_dict(_report(reference, target, pockets=_pocket_report(reference, target)))
    payload["pockets"]["comparisons"][0]["associated_mutations"] = ["p.A10G"]
    with pytest.raises(ValueError, match="mutation|reconstruct"):
        deserialize_report(payload, snapshots=(reference, target))


def test_typed_pocket_snapshots_reject_invalid_values_and_lineage() -> None:
    reference = _snapshot("reference.pdb", b"reference")
    valid = _pocket_report(reference)
    detection = valid.detections[0]
    volume = valid.volumes[0]
    volume_result = volume.result
    assert volume_result is not None
    lining = valid.lining_residues[0]
    candidate = detection.candidates[0]

    invalid_constructors: tuple[tuple[str, Any], ...] = (
        ("role", lambda: PocketDetectionSnapshot("other", Availability.AVAILABLE)),
        (
            "candidate sequence",
            lambda: PocketDetectionSnapshot("reference", Availability.AVAILABLE, cast(Any, "bad")),
        ),
        (
            "candidate type",
            lambda: PocketDetectionSnapshot("reference", Availability.AVAILABLE, candidates=cast(Any, (object(),))),
        ),
        (
            "detection type",
            lambda: PocketDetectionSnapshot("reference", Availability.AVAILABLE, detection=cast(Any, object())),
        ),
        (
            "detection availability",
            lambda: PocketDetectionSnapshot(
                "reference",
                Availability.NOT_APPLICABLE,
                detection=detection.detection,
            ),
        ),
        (
            "settings type",
            lambda: PocketDetectionSnapshot("reference", Availability.AVAILABLE, settings=cast(Any, object())),
        ),
        (
            "diagnostics type",
            lambda: PocketDetectionSnapshot("reference", Availability.AVAILABLE, diagnostics=cast(Any, (object(),))),
        ),
        (
            "counts type",
            lambda: PocketDetectionSnapshot("reference", Availability.AVAILABLE, counts=cast(Any, [])),
        ),
        (
            "count key",
            lambda: PocketDetectionSnapshot("reference", Availability.AVAILABLE, counts={"": 1}),
        ),
        (
            "count value",
            lambda: PocketDetectionSnapshot("reference", Availability.AVAILABLE, counts={"count": -1}),
        ),
        (
            "provenance type",
            lambda: PocketDetectionSnapshot("reference", Availability.AVAILABLE, provenance=cast(Any, object())),
        ),
        ("volume role", lambda: PocketVolumeSnapshot("other")),
        ("volume result type", lambda: PocketVolumeSnapshot("reference", cast(Any, object()))),
        (
            "volume availability",
            lambda: PocketVolumeSnapshot("reference", volume.result, Availability.NOT_APPLICABLE),
        ),
        ("lining id", lambda: PocketLiningSnapshot("reference", "")),
        (
            "lining residues",
            lambda: PocketLiningSnapshot("reference", candidate.candidate_id, residues=cast(Any, (object(),))),
        ),
        (
            "lining diagnostics",
            lambda: PocketLiningSnapshot("reference", candidate.candidate_id, diagnostics=cast(Any, (object(),))),
        ),
    )
    for _name, constructor in invalid_constructors:
        with pytest.raises((TypeError, ValueError), match="."):
            constructor()

    with pytest.raises(TypeError, match="detections"):
        PocketReportSnapshot(detections=cast(Any, (object(),)))
    with pytest.raises(TypeError, match="volumes"):
        PocketReportSnapshot(volumes=cast(Any, (object(),)))
    with pytest.raises(TypeError, match="comparisons"):
        PocketReportSnapshot(comparisons=cast(Any, (object(),)))
    with pytest.raises(TypeError, match="lining"):
        PocketReportSnapshot(lining_residues=cast(Any, (object(),)))
    with pytest.raises(TypeError, match="diagnostics"):
        PocketReportSnapshot(diagnostics=cast(Any, (object(),)))
    with pytest.raises(TypeError, match="matching"):
        PocketReportSnapshot(matching=cast(Any, object()))
    with pytest.raises(TypeError, match="concordance"):
        PocketReportSnapshot(concordance=cast(Any, object()))
    with pytest.raises(TypeError, match="provenance"):
        PocketReportSnapshot(provenance=cast(Any, object()))
    with pytest.raises(ValueError, match="unique"):
        PocketReportSnapshot(detections=(detection, detection))
    with pytest.raises(ValueError, match="candidate_id"):
        PocketReportSnapshot(
            detections=(detection,),
            volumes=(PocketVolumeSnapshot("target", replace(volume_result, candidate_id="unrelated")),),
        )
    with pytest.raises(ValueError, match="lining"):
        PocketReportSnapshot(
            detections=(detection,),
            lining_residues=(PocketLiningSnapshot("target", "unrelated"),),
        )
    with pytest.raises(ValueError, match="comparison"):
        unrelated_sphere = replace(candidate.alpha_spheres[0], center_xyz=(9.0, 9.0, 9.0))
        unrelated = PocketCandidate((unrelated_sphere,), "0" * 64, "1" * 64)
        PocketReportSnapshot(
            detections=(detection,),
            comparisons=(
                PocketComparison(
                    Availability.AVAILABLE,
                    PocketMatch(unrelated, unrelated, PocketMatchState.MATCHED),
                ),
            ),
        )
    assert PocketReportSnapshot().availability is Availability.NOT_APPLICABLE
    assert volume.to_json()["availability"] == Availability.AVAILABLE.value
    assert lining.to_json()["candidate_id"] == candidate.candidate_id
