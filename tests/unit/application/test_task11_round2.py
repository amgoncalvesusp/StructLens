from __future__ import annotations

import importlib
import json
from dataclasses import replace
from pathlib import Path
from typing import Any

import pytest
from openpyxl import load_workbook

import structlens.application.report_exports as report_exports
import structlens.core.reports.schema as schema_module
from structlens.application.export_service import export_report_csv, export_report_xlsx
from structlens.application.project_state import ProjectState
from structlens.application.report_serialization import (
    deserialize_report,
    report_to_dict,
    save_report,
    serialize_report,
)
from structlens.core.evidence import Availability, PocketEvidenceChannel, PocketEvidenceSections
from structlens.core.models import MutationEvent, MutationKind, ResidueId
from structlens.core.reports.safe_io import read_bounded_bytes
from structlens.core.reports.schema import validate_analysis_report_payload

fixtures = importlib.import_module("tests.unit.application.test_task11_round1")


def _mutation() -> MutationEvent:
    return MutationEvent(
        4,
        MutationKind.SUBSTITUTION,
        ResidueId("reference", "1", "A", "10", None, "ALA"),
        ResidueId("target", "1", "A", "10", None, "VAL"),
        "A",
        "V",
        "A10",
        "V10",
        "p.A10V",
        -1,
        64,
        "hydrophobic",
    )


def _mutation_report() -> tuple[Any, Any, Any, MutationEvent]:
    reference = fixtures._snapshot("reference.pdb", b"reference")
    target = fixtures._snapshot("target.pdb", b"target")
    report = fixtures._report(reference, target, pockets=fixtures._pocket_report(reference, target))
    event = _mutation()
    pockets = report.pockets
    assert pockets is not None
    comparison = replace(pockets.comparisons[0], associated_mutations=(event,))
    report = replace(report, pockets=replace(pockets, comparisons=(comparison, *pockets.comparisons[1:])))
    return report, reference, target, event


def test_associated_mutations_round_trip_as_complete_typed_events() -> None:
    report, reference, target, event = _mutation_report()
    restored = deserialize_report(
        serialize_report(report, snapshots=(reference, target)), snapshots=(reference, target)
    )
    assert restored.pockets is not None
    assert restored.pockets.comparisons[0].associated_mutations == (event,)


def test_mutation_wire_rejects_unknown_fields_and_noncanonical_numbers() -> None:
    report, reference, target, _event = _mutation_report()
    payload = report_to_dict(report)
    mutation = payload["pockets"]["comparisons"][0]["associated_mutations"][0]
    mutation["unexpected"] = True
    with pytest.raises(ValueError, match="canonical|unknown"):
        deserialize_report(payload, snapshots=(reference, target))

    payload = report_to_dict(report)
    payload["pockets"]["comparisons"][0]["associated_mutations"][0]["blosum62_score"] = -1.0
    with pytest.raises(ValueError, match="canonical|round-trip"):
        deserialize_report(payload, snapshots=(reference, target))


def test_pocket_candidate_lineage_compares_identity_and_source_tuple() -> None:
    report, reference, target, _event = _mutation_report()
    pockets = report.pockets
    assert pockets is not None
    target_candidate = pockets.detections[1].candidates[0]
    forged = replace(target_candidate, source_content_id="0" * 64)
    forged_match = replace(pockets.matching.matches[0], target_candidate=forged)
    forged_matching = replace(pockets.matching, matches=(forged_match,))
    with pytest.raises(ValueError, match="lineage"):
        replace(report, pockets=replace(pockets, matching=forged_matching))


def test_pocket_export_rows_include_auditable_match_and_comparison_fields() -> None:
    report, reference, target, event = _mutation_report()
    rows = report_exports._rows(report)
    match_rows = {row["metric"]: row for row in rows if row["section"] == "pocket_matches"}
    assert (
        match_rows["reference_candidate_id"]["value"]
        == report.pockets.matching.matches[0].reference_candidate.candidate_id
    )
    assert (
        match_rows["target_candidate_id"]["value"] == report.pockets.matching.matches[0].target_candidate.candidate_id
    )
    assert match_rows["score"]["value"] == 0.75
    assert match_rows["alternative_assignment_score"]["value"] == 0.2

    comparison = report.pockets.comparisons[0]
    comparison_role = (
        f"{comparison.match.reference_candidate.candidate_id}->{comparison.match.target_candidate.candidate_id}"
    )
    comparison_rows = {
        row["metric"]: row
        for row in rows
        if row["section"] == "pocket_comparisons" and row["role"] == comparison_role and row["status"] == "available"
    }
    assert comparison_rows["reference_volume"]["value"] == comparison.reference_volume_angstrom3
    assert comparison_rows["target_volume"]["value"] == comparison.target_volume_angstrom3
    assert comparison_rows["volume_delta"]["value"] == comparison.volume_delta_angstrom3
    assert json.loads(comparison_rows["associated_mutations"]["value"]) == [
        {
            "alignment_index": event.alignment_index,
            "kind": event.kind.value,
            "reference": {
                "structure_id": event.reference.structure_id,
                "model_id": event.reference.model_id,
                "chain_id": event.reference.chain_id,
                "auth_seq_id": event.reference.auth_seq_id,
                "insertion_code": event.reference.insertion_code,
                "residue_name": event.reference.residue_name,
            },
            "target": {
                "structure_id": event.target.structure_id,
                "model_id": event.target.model_id,
                "chain_id": event.target.chain_id,
                "auth_seq_id": event.target.auth_seq_id,
                "insertion_code": event.target.insertion_code,
                "residue_name": event.target.residue_name,
            },
            "reference_aa": event.reference_aa,
            "target_aa": event.target_aa,
            "reference_label": event.reference_label,
            "target_label": event.target_label,
            "canonical_notation": event.canonical_notation,
            "blosum62_score": event.blosum62_score,
            "grantham_distance": event.grantham_distance,
            "physicochemical_class": event.physicochemical_class,
        }
    ]
    lining_rows = [row for row in rows if row["section"] == "lining_residues" and row["metric"] == "residue"]
    assert lining_rows and all(row["role"] == "reference" for row in lining_rows)
    assert all(row["provenance"] == comparison.match.reference_candidate.candidate_id for row in lining_rows)


def test_concordance_export_keeps_native_units_and_unavailable_status() -> None:
    reference = fixtures._snapshot("reference.pdb", b"reference")
    target = fixtures._snapshot("target.pdb", b"target")
    report = fixtures._report(reference, target, pockets=fixtures._pocket_report(reference, target))
    assert report.pockets is not None
    sections = PocketEvidenceSections(
        geometry=PocketEvidenceChannel(
            Availability.INVALID_INPUT,
            measure=1.5,
            units={"distance": "angstrom"},
        ),
        ligand_support=Availability.NOT_APPLICABLE,
        parameter_persistence=Availability.NOT_APPLICABLE,
        volume_sensitivity=Availability.NOT_APPLICABLE,
        match_ambiguity=Availability.NOT_APPLICABLE,
        qc=Availability.NOT_APPLICABLE,
        interactions=Availability.NOT_APPLICABLE,
    )
    rich = replace(report, pockets=replace(report.pockets, concordance=sections))
    geometry = next(row for row in report_exports._rows(rich) if row["metric"] == "geometry")
    assert geometry["units"] == "angstrom"
    assert geometry["status"] == "invalid_input"
    assert geometry["reason"] == "invalid input"


def test_project_state_save_allows_atomic_update_of_existing_file(tmp_path: Path) -> None:
    path = tmp_path / "project.json"
    first = ProjectState(reference_source="reference.pdb")
    first.save(path)
    second = replace(first, reference_source="updated.pdb")
    second.save(path)
    assert ProjectState.from_json(path.read_text(encoding="utf-8")).reference_source == "updated.pdb"


def test_safe_reads_reject_linked_ancestor(tmp_path: Path) -> None:
    real = tmp_path / "real"
    real.mkdir()
    payload = real / "payload.bin"
    payload.write_bytes(b"payload")
    link = tmp_path / "linked"
    try:
        link.symlink_to(real, target_is_directory=True)
    except (OSError, NotImplementedError):
        pytest.skip("test host cannot create directory symlinks")
    with pytest.raises(ValueError, match="link"):
        read_bounded_bytes(link / "payload.bin", max_bytes=100, label="payload")


def test_schema_fallback_rejects_unknown_nested_object(monkeypatch: pytest.MonkeyPatch) -> None:
    report = fixtures._report(
        fixtures._snapshot("reference.pdb", b"reference"), fixtures._snapshot("target.pdb", b"target")
    )
    payload = report_to_dict(report)
    payload["availability"]["unexpected"] = "bad"

    real_import = schema_module.importlib.import_module

    def missing_validator(name: str) -> object:
        if name == "jsonschema":
            raise ImportError("optional dependency absent")
        return real_import(name)

    monkeypatch.setattr(schema_module.importlib, "import_module", missing_validator)
    with pytest.raises(ValueError, match="unknown"):
        validate_analysis_report_payload(payload)


def test_exports_are_tidy_formula_safe_and_byte_deterministic(tmp_path: Path) -> None:
    reference = fixtures._snapshot("reference.pdb", b"reference")
    target = fixtures._snapshot("target.pdb", b"target")
    report = fixtures._report(reference, target, pockets=fixtures._pocket_report(reference, target))
    first = tmp_path / "first.xlsx"
    second = tmp_path / "second.xlsx"
    export_report_xlsx(report, first, snapshots=(reference, target))
    export_report_xlsx(report, second, snapshots=(reference, target))
    assert first.read_bytes() == second.read_bytes()
    workbook = load_workbook(first, data_only=False, read_only=True)
    assert (
        "free volume" in str(workbook["Pockets"]["A1"].value).lower()
        or "free volume" in " ".join(str(cell.value) for cell in workbook["Pockets"][1]).lower()
    )
    assert workbook["Diagnostics"]["D2"].value == "' =formula"
    csv_path = tmp_path / "report.csv"
    export_report_csv(report, csv_path, snapshots=(reference, target))
    text = csv_path.read_text(encoding="utf-8")
    assert "pocket" in text
    assert "' =formula" in text


def test_project_legacy_state_remains_unverified(tmp_path: Path) -> None:
    project = ProjectState.from_dict({"schema_version": "3.0", "settings": {}})
    assert getattr(project, "report_verification", "legacy_unverified") == "legacy_unverified"
    path = tmp_path / "project.json"
    project.save(path)
    assert json.loads(path.read_text(encoding="utf-8"))["v04"]["report_verification"] == "legacy_unverified"


def test_project_report_verification_checks_id_and_snapshots(tmp_path: Path) -> None:
    reference = fixtures._snapshot("reference.pdb", b"reference")
    target = fixtures._snapshot("target.pdb", b"target")
    report = fixtures._report(reference, target)
    report_path = tmp_path / "report.json"
    save_report(report, report_path, snapshots=(reference, target))
    state = ProjectState(report_path=str(report_path), report_hash=report.report_id)
    assert state.verify_report(snapshots=(reference, target)) is True
    with pytest.raises(Exception, match="ID|hash|verification"):
        ProjectState(report_path=str(report_path), report_hash="0" * 64).verify_report(snapshots=(reference, target))
