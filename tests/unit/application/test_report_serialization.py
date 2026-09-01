from __future__ import annotations

import hashlib
import json
from copy import deepcopy
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, cast

import numpy as np
import pytest

from structlens.application.report_serialization import (
    deserialize_report,
    load_report,
    load_snapshot,
    report_to_dict,
    save_report,
    serialize_report,
    store_snapshot,
)
from structlens.core.difference_maps import DistanceDifferenceMatrix, ResidueDisplacementVector
from structlens.core.evidence import (
    Availability,
    Diagnostic,
    DiagnosticSeverity,
    EvidenceCard,
    EvidenceQuality,
    InteractionEvidence,
    PocketEvidenceSections,
    SequenceEvidence,
    SiteEvidence,
    StructureEvidence,
)
from structlens.core.interactions import (
    InteractionChange,
    InteractionDifference,
    InteractionRecord,
    InteractionType,
    ReferenceInteractionKey,
)
from structlens.core.models import AnalysisResult, CorrespondenceStatus, ResidueCorrespondence, ResidueId
from structlens.core.msa import (
    AnalysisSequence,
    MSAColumn,
    MSAResidueCell,
    MultipleSequenceAlignment,
    SequenceResidueRef,
)
from structlens.core.parsing import InputSelection, SourceSnapshot, StructureFormat
from structlens.core.pockets.matching import PocketMatchingResult
from structlens.core.provenance import AuditEvent, MethodProvenance
from structlens.core.quality import StructureQualityReport
from structlens.core.reports import (
    AnalysisReport,
    AnalysisSnapshot,
    DistanceMapSnapshot,
    InputQualityBundle,
    PocketReportSnapshot,
    SectionAvailability,
)
from structlens.core.sites import SiteDefinition, SiteDefinitionMode, SiteMetrics


def _snapshot(name: str, content: bytes) -> SourceSnapshot:
    digest = hashlib.sha256(content).hexdigest()
    return SourceSnapshot(content, digest, digest, digest, name, "pdb", False)


def _report(reference: SourceSnapshot, target: SourceSnapshot) -> AnalysisReport:
    ref = InputSelection(reference.content_id, reference.display_name, StructureFormat.PDB, "1")
    tgt = InputSelection(target.content_id, target.display_name, StructureFormat.PDB, "1")
    return AnalysisReport(
        ref,
        tgt,
        InputQualityBundle(
            StructureQualityReport(Availability.NOT_APPLICABLE), StructureQualityReport(Availability.NOT_APPLICABLE)
        ),
        availability=SectionAvailability(),
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
    )


def test_canonical_report_round_trip_is_deterministic_and_excludes_audit_time(tmp_path: Path) -> None:
    reference = _snapshot("reference.pdb", b"reference")
    target = _snapshot("target.pdb", b"target")
    report = _report(reference, target)
    first = serialize_report(report, snapshots=(reference, target))
    second = serialize_report(
        report, snapshots=(reference, target), audit_event=AuditEvent("export", datetime(2030, 1, 1, tzinfo=UTC))
    )
    assert first == second
    assert json.loads(first)["report_id"] == report.report_id

    report_path = tmp_path / "report.json"
    save_report(report, report_path, snapshots=(reference, target))
    restored = load_report(report_path, snapshots=(reference, target))
    assert restored.report_id == report.report_id
    assert restored.reference_selection.selection_id == report.reference_selection.selection_id
    assert restored.target_selection.selection_id == report.target_selection.selection_id


def test_snapshot_store_is_content_addressed_and_export_rejects_mismatch(tmp_path: Path) -> None:
    source = tmp_path / "reference.pdb"
    source.write_bytes(b"original")
    snapshot = _snapshot(source.name, source.read_bytes())
    store_snapshot(snapshot, tmp_path / "evidence")
    assert load_snapshot(snapshot.content_id, tmp_path / "evidence").content_id == snapshot.content_id
    report = _report(snapshot, _snapshot("target.pdb", b"target"))
    source.write_bytes(b"replacement")
    with pytest.raises(ValueError, match="hash|snapshot"):
        save_report(report, tmp_path / "report.json", source_paths=(source,), snapshot_dir=tmp_path / "evidence")


def test_snapshot_store_rejects_tampered_metadata(tmp_path: Path) -> None:
    reference = _snapshot("reference.pdb", b"reference")
    target = _snapshot("target.pdb", b"target")
    evidence = tmp_path / "evidence"
    store_snapshot(reference, evidence)
    store_snapshot(target, evidence)
    metadata_path = evidence / f"{reference.content_id}.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata["decompressed_sha256"] = "0" * 64
    metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
    with pytest.raises(ValueError, match="metadata|hash"):
        serialize_report(_report(reference, target), snapshot_dir=evidence)


def test_report_serializer_requires_typed_inputs_and_idempotent_snapshot_store(tmp_path: Path) -> None:
    reference = _snapshot("reference.pdb", b"reference")
    target = _snapshot("target.pdb", b"target")
    report = _report(reference, target)
    with pytest.raises(TypeError, match="AnalysisReport"):
        serialize_report(cast(Any, object()), snapshots=(reference, target))
    with pytest.raises(TypeError, match="AuditEvent"):
        serialize_report(report, snapshots=(reference, target), audit_event=cast(Any, object()))

    evidence = tmp_path / "evidence"
    first = serialize_report(report, snapshots=(reference, target), snapshot_dir=evidence)
    second = serialize_report(report, snapshot_dir=evidence)
    assert first == second
    with pytest.raises(ValueError, match="report_id"):
        payload = json.loads(first)
        payload["report_id"] = "0" * 64
        deserialize_report(payload, snapshots=(reference, target))

    malformed = json.loads(first)
    malformed["reference_selection"]["selection_id"] = "0" * 64
    with pytest.raises(ValueError, match="selection_id"):
        deserialize_report(malformed, snapshots=(reference, target))


def test_rich_typed_report_round_trips_all_serialized_sections() -> None:
    reference = _snapshot("reference.pdb", b"reference")
    target = _snapshot("target.pdb", b"target")
    reference_residue = ResidueId("reference", "1", "A", "10", None, "ALA")
    target_residue = ResidueId("target", "1", "A", "10", None, "SER")
    correspondence = ResidueCorrespondence(
        0,
        reference_residue,
        target_residue,
        "A",
        "S",
        CorrespondenceStatus.SUBSTITUTION,
        ca_displacement_angstrom=1.0,
    )
    analysis = AnalysisResult(
        "reference",
        "target",
        (correspondence,),
        (),
        0.5,
        1.0,
        "sequence",
        strict_rmsd_angstrom=1.0,
        mapped_residue_count=1,
    )
    reference_record = InteractionRecord(
        "reference", InteractionType.HYDROPHOBIC, reference_residue, None, "CA", None, 4.0
    )
    target_record = InteractionRecord("target", InteractionType.HYDROPHOBIC, target_residue, None, "CA", None, 4.1)
    difference = InteractionDifference(
        ReferenceInteractionKey(InteractionType.HYDROPHOBIC, "A:10"),
        InteractionChange.CONSERVED,
        reference_record,
        target_record,
    )
    sequence_ref = SequenceResidueRef(0, "A", reference_residue)
    msa = MultipleSequenceAlignment(
        (AnalysisSequence("reference", "A", "A", (sequence_ref,), "structure"),),
        (("reference", "A"),),
        (
            MSAColumn(
                0,
                "A:10",
                reference_residue,
                (MSAResidueCell("reference", 0, sequence_ref, "A"),),
                1,
                0.0,
                0.0,
                1.0,
                0.0,
            ),
        ),
        "reference",
        "fallback",
    )
    site = SiteMetrics("active", "target", 1, 1.0, global_frame_backbone_rmsd_angstrom=0.2)
    vector = ResidueDisplacementVector(
        "A:10", reference_residue, target_residue, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 0.0, 1.0), 1.0
    )
    card = EvidenceCard(
        sequence_ref,
        "target",
        sequence=SequenceEvidence("A", "S", 0, 0.5, 0.5, 0.0, 0.0, 0.0, (sequence_ref,)),
        structure=StructureEvidence(ca_displacement_angstrom=1.0, available=True),
        interactions=InteractionEvidence((difference,), (reference_record,), (target_record,)),
        site=SiteEvidence((site,)),
        quality=EvidenceQuality("available", ("sequence", "structure", "sites"), source_count=3),
        pocket_sections=PocketEvidenceSections(
            Availability.AVAILABLE,
            Availability.NOT_APPLICABLE,
            Availability.NOT_APPLICABLE,
            Availability.NOT_APPLICABLE,
            Availability.NOT_APPLICABLE,
            Availability.NOT_APPLICABLE,
            Availability.NOT_APPLICABLE,
        ),
    )
    report = AnalysisReport(
        InputSelection(reference.content_id, reference.display_name, StructureFormat.PDB, "1"),
        InputSelection(target.content_id, target.display_name, StructureFormat.PDB, "1"),
        InputQualityBundle(
            StructureQualityReport(Availability.AVAILABLE), StructureQualityReport(Availability.AVAILABLE)
        ),
        analysis=AnalysisSnapshot.from_result(analysis),
        msa=msa,
        interactions=InteractionEvidence((difference,), (reference_record,), (target_record,)),
        sites=SiteEvidence((site,)),
        distance_map=DistanceMapSnapshot.from_matrix(
            DistanceDifferenceMatrix(
                ("A:10",), np.zeros((1, 1)), np.zeros((1, 1)), np.zeros((1, 1)), np.ones((1, 1), dtype=bool)
            )
        ),
        displacement_vectors=(vector,),
        evidence_cards=(card,),
        site_definitions=(
            SiteDefinition("active", mode=SiteDefinitionMode.KEY_RESIDUES, reference_residues=(reference_residue,)),
        ),
        availability=SectionAvailability(
            input_quality=Availability.AVAILABLE,
            analysis=Availability.AVAILABLE,
            msa=Availability.AVAILABLE,
            interactions=Availability.AVAILABLE,
            sites=Availability.AVAILABLE,
            distance_map=Availability.AVAILABLE,
            displacement_vectors=Availability.AVAILABLE,
            evidence_cards=Availability.AVAILABLE,
        ),
        provenance=MethodProvenance(
            "test.report",
            "1",
            input_hashes={
                "reference_content": reference.content_id,
                "target_content": target.content_id,
                "reference_raw": reference.raw_sha256,
                "target_raw": target.raw_sha256,
            },
        ),
    )
    restored = deserialize_report(
        serialize_report(report, snapshots=(reference, target)), snapshots=(reference, target)
    )
    assert restored.report_id == report.report_id
    assert restored.msa is not None and restored.msa.columns[0].cells[0].character == "A"
    assert restored.interactions is not None and restored.sites is not None
    assert restored.evidence_cards[0].quality.source_count == 3
    assert restored.evidence_cards[0].pocket_sections is not None

    canonical = report_to_dict(report)
    malformed = deepcopy(canonical)
    malformed["evidence_cards"][0]["pocket_concordance"]["unexpected"] = {}
    with pytest.raises(ValueError, match="unknown"):
        deserialize_report(malformed, snapshots=(reference, target))
    malformed = deepcopy(canonical)
    malformed["evidence_cards"][0]["pocket_concordance"].pop("qc")
    with pytest.raises(ValueError, match="seven"):
        deserialize_report(malformed, snapshots=(reference, target))
    malformed = deepcopy(canonical)
    malformed["displacement_vectors"][0]["reference_residue"] = None
    with pytest.raises(ValueError, match="reference_residue"):
        deserialize_report(malformed, snapshots=(reference, target))
    malformed = deepcopy(canonical)
    malformed["displacement_vectors"][0]["start_xyz"] = [0.0, 0.0]
    with pytest.raises(ValueError, match="start_xyz"):
        deserialize_report(malformed, snapshots=(reference, target))
    malformed = deepcopy(canonical)
    malformed["evidence_cards"][0]["interactions"] = []
    with pytest.raises(ValueError, match="interactions"):
        deserialize_report(malformed, snapshots=(reference, target))
    malformed = deepcopy(canonical)
    malformed["evidence_cards"][0]["site"] = []
    with pytest.raises(ValueError, match="site"):
        deserialize_report(malformed, snapshots=(reference, target))
    malformed = deepcopy(canonical)
    malformed["input_quality"]["reference"] = "not an object"
    with pytest.raises(ValueError, match="reference"):
        deserialize_report(malformed, snapshots=(reference, target))
    malformed = deepcopy(canonical)
    malformed["analysis"]["correspondences"] = "not an array"
    with pytest.raises(ValueError, match="correspondences"):
        deserialize_report(malformed, snapshots=(reference, target))


def test_pocket_report_infers_availability_from_matching_and_rejects_empty_available() -> None:
    inferred = PocketReportSnapshot(matching=PocketMatchingResult(Availability.AVAILABLE))
    assert inferred.availability is Availability.AVAILABLE
    with pytest.raises(ValueError, match="evidence"):
        PocketReportSnapshot(availability=Availability.AVAILABLE)


def test_non_available_typed_pocket_payload_preserves_diagnostic_missingness() -> None:
    reference = _snapshot("reference.pdb", b"reference")
    target = _snapshot("target.pdb", b"target")
    unavailable = PocketReportSnapshot(
        availability=Availability.INVALID_INPUT,
        diagnostics=(Diagnostic("POCKET_INPUT", DiagnosticSeverity.ERROR, "invalid pocket input"),),
    )
    report = AnalysisReport(
        InputSelection(reference.content_id, reference.display_name, StructureFormat.PDB, "1"),
        InputSelection(target.content_id, target.display_name, StructureFormat.PDB, "1"),
        InputQualityBundle(
            StructureQualityReport(Availability.NOT_APPLICABLE), StructureQualityReport(Availability.NOT_APPLICABLE)
        ),
        availability=SectionAvailability(pockets=Availability.INVALID_INPUT),
        pockets=unavailable,
    )
    assert report.pockets is unavailable
    assert report.pockets.diagnostics[0].code == "POCKET_INPUT"
