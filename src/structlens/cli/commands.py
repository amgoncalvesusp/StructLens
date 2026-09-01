"""Application-level orchestration for the canonical StructLens CLI."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from structlens.application.dto import AnalysisReportRequest
from structlens.application.pocket_comparison_service import PocketComparisonReport, PocketComparisonService
from structlens.application.pocket_service import PocketDetectionReport, StructurePocketService
from structlens.application.quality_service import StructureQualityService
from structlens.application.report_exports import (
    export_report_csv,
    export_report_json,
    export_report_tsv,
    export_report_xlsx,
)
from structlens.application.report_service import ReportService
from structlens.core.errors import StructLensError
from structlens.core.evidence import Availability, InteractionEvidence
from structlens.core.models import AlignmentMode, AnalysisSettings, ResidueCorrespondence
from structlens.core.parsing import ParsedStructure, SourceSnapshot
from structlens.core.pockets import PocketVolumeResult
from structlens.core.reports.safe_io import atomic_write_bytes

from .inputs import CapturedInput, capture_input
from .rendering import (
    json_bytes,
    render_compare,
    render_pocket_compare,
    render_pockets,
    render_quality,
    selection_payload,
    source_payload,
)


@dataclass(frozen=True, slots=True)
class CliServices:
    """Injectable application services; scientific logic stays outside the CLI."""

    report_service: ReportService = field(default_factory=ReportService)
    quality_service: StructureQualityService = field(default_factory=StructureQualityService)
    pocket_service: StructurePocketService = field(default_factory=StructurePocketService)
    comparison_service: PocketComparisonService = field(default_factory=PocketComparisonService)


@dataclass(frozen=True, slots=True)
class WorkflowResult:
    status: str
    payload: dict[str, Any] = field(default_factory=dict)
    summary: str = ""
    exit_code: int = 0
    error: str | None = None
    report: Any | None = None
    snapshots: tuple[SourceSnapshot, ...] = ()

    def json_bytes(self) -> bytes:
        return json_bytes(self.payload)


def run_compare(
    reference_path: str | Path,
    target_path: str | Path,
    *,
    mode: str | AlignmentMode = AlignmentMode.AUTO,
    reference_model: str = "1",
    target_model: str = "1",
    reference_author_chains: tuple[str, ...] = (),
    target_author_chains: tuple[str, ...] = (),
    reference_label_chains: tuple[str, ...] = (),
    target_label_chains: tuple[str, ...] = (),
    reference_altloc_policy: str = "highest_occupancy",
    target_altloc_policy: str = "highest_occupancy",
    reference_assembly_scope: str = "asymmetric_unit",
    target_assembly_scope: str = "asymmetric_unit",
    json_path: str | Path | None = None,
    report_path: str | Path | None = None,
    csv_path: str | Path | None = None,
    tsv_path: str | Path | None = None,
    xlsx_path: str | Path | None = None,
    services: CliServices | None = None,
) -> WorkflowResult:
    try:
        alignment_mode = _alignment_mode(mode)
        if alignment_mode is AlignmentMode.MANUAL:
            raise ValueError("manual mode requires explicit residue-pair input; no manual pair input is configured")
        if json_path is not None and report_path is not None:
            raise ValueError("use either json_path or report_path, not both")
        _validate_outputs(json_path, report_path, csv_path, tsv_path, xlsx_path)
        reference = capture_input(
            reference_path,
            model_id=reference_model,
            author_chain_ids=reference_author_chains,
            label_chain_ids=reference_label_chains,
            altloc_policy=reference_altloc_policy,
            assembly_scope=reference_assembly_scope,
        )
        target = capture_input(
            target_path,
            model_id=target_model,
            author_chain_ids=target_author_chains,
            label_chain_ids=target_label_chains,
            altloc_policy=target_altloc_policy,
            assembly_scope=target_assembly_scope,
        )
        _require_pairwise_chains(reference.parsed, target.parsed)
        request = AnalysisReportRequest(
            reference_snapshot=reference.snapshot,
            target_snapshot=target.snapshot,
            reference_selection=reference.selection,
            target_selection=target.selection,
            analysis_settings=AnalysisSettings(alignment_mode=alignment_mode),
        )
        active = services or CliServices()
        report = active.report_service.analyze(request)
        snapshots = (reference.snapshot, target.snapshot)
        if json_path is not None or report_path is not None:
            export_report_json(report, report_path or json_path, snapshots=snapshots)  # type: ignore[arg-type]
        if csv_path is not None:
            export_report_csv(report, csv_path, snapshots=snapshots)
        if tsv_path is not None:
            export_report_tsv(report, tsv_path, snapshots=snapshots)
        if xlsx_path is not None:
            export_report_xlsx(report, xlsx_path, snapshots=snapshots)
        return WorkflowResult(
            status=report.availability.analysis.value,
            # Keep stdout/legacy invocation usable even when an optional MSA
            # section contains a report-schema edge value such as ``-0.0``.
            # Canonical file outputs still go through ``export_report_*`` and
            # therefore retain the strict Task 11 schema gate.
            payload=report.to_json(),
            summary=render_compare(report),
            exit_code=0 if report.availability.analysis in _SUCCESS_STATES else 2,
            report=report,
            snapshots=snapshots,
        )
    except (StructLensError, OSError, ValueError) as error:
        return _failure(error)


def run_qc(
    path: str | Path,
    *,
    model: str = "1",
    author_chains: tuple[str, ...] = (),
    label_chains: tuple[str, ...] = (),
    altloc_policy: str = "highest_occupancy",
    assembly_scope: str = "asymmetric_unit",
    json_path: str | Path | None = None,
    services: CliServices | None = None,
) -> WorkflowResult:
    try:
        captured = capture_input(
            path,
            model_id=model,
            author_chain_ids=author_chains,
            label_chain_ids=label_chains,
            altloc_policy=altloc_policy,
            assembly_scope=assembly_scope,
        )
        quality = (services or CliServices()).quality_service.analyze(captured.parsed)
        payload = {
            "selection": selection_payload(captured.selection),
            "source": source_payload(captured.snapshot),
            "quality": quality.to_json(),
        }
        if json_path is not None:
            _write_json(json_path, payload)
        return WorkflowResult(
            quality.availability.value,
            payload,
            render_quality(quality),
            exit_code=_exit_code(quality.availability),
            snapshots=(captured.snapshot,),
        )
    except (StructLensError, OSError, ValueError) as error:
        return _failure(error)


def run_pockets(
    path: str | Path,
    *,
    model: str = "1",
    author_chains: tuple[str, ...] = (),
    label_chains: tuple[str, ...] = (),
    altloc_policy: str = "highest_occupancy",
    assembly_scope: str = "asymmetric_unit",
    json_path: str | Path | None = None,
    services: CliServices | None = None,
) -> WorkflowResult:
    try:
        captured = capture_input(
            path,
            model_id=model,
            author_chain_ids=author_chains,
            label_chain_ids=label_chains,
            altloc_policy=altloc_policy,
            assembly_scope=assembly_scope,
        )
        active = services or CliServices()
        detection = active.pocket_service.analyze(captured.parsed)
        volumes = _measure_volumes(active.pocket_service, captured.parsed, detection)
        workflow_status = _aggregate_availability(detection.availability, *(item.availability for item in volumes))
        payload = _single_pocket_payload(captured, detection, volumes, workflow_status)
        if json_path is not None:
            _write_json(json_path, payload)
        return WorkflowResult(
            workflow_status.value,
            payload,
            render_pockets(detection, volumes=volumes, status=workflow_status),
            exit_code=_exit_code(workflow_status),
            snapshots=(captured.snapshot,),
        )
    except (StructLensError, OSError, ValueError) as error:
        return _failure(error)


def run_pocket_compare(
    reference_path: str | Path,
    target_path: str | Path,
    *,
    mode: str | AlignmentMode = AlignmentMode.AUTO,
    reference_model: str = "1",
    target_model: str = "1",
    reference_author_chains: tuple[str, ...] = (),
    target_author_chains: tuple[str, ...] = (),
    reference_label_chains: tuple[str, ...] = (),
    target_label_chains: tuple[str, ...] = (),
    reference_altloc_policy: str = "highest_occupancy",
    target_altloc_policy: str = "highest_occupancy",
    reference_assembly_scope: str = "asymmetric_unit",
    target_assembly_scope: str = "asymmetric_unit",
    json_path: str | Path | None = None,
    services: CliServices | None = None,
) -> WorkflowResult:
    try:
        alignment_mode = _alignment_mode(mode)
        if alignment_mode is AlignmentMode.MANUAL:
            raise ValueError("manual mode requires explicit residue-pair input; no manual pair input is configured")
        reference = capture_input(
            reference_path,
            model_id=reference_model,
            author_chain_ids=reference_author_chains,
            label_chain_ids=reference_label_chains,
            altloc_policy=reference_altloc_policy,
            assembly_scope=reference_assembly_scope,
        )
        target = capture_input(
            target_path,
            model_id=target_model,
            author_chain_ids=target_author_chains,
            label_chain_ids=target_label_chains,
            altloc_policy=target_altloc_policy,
            assembly_scope=target_assembly_scope,
        )
        _require_pairwise_chains(reference.parsed, target.parsed)
        active = services or CliServices()
        request = AnalysisReportRequest(
            reference.snapshot,
            target.snapshot,
            reference.selection,
            target.selection,
            analysis_settings=AnalysisSettings(alignment_mode=alignment_mode),
        )
        report = active.report_service.analyze(request)
        if report.analysis is None:
            raise ValueError("pairwise analysis was unavailable; inspect the report diagnostics")
        reference_detection = active.pocket_service.analyze(reference.parsed)
        target_detection = active.pocket_service.analyze(target.parsed)
        reference_volumes = _volume_map(active.pocket_service, reference.parsed, reference_detection)
        target_volumes = _volume_map(active.pocket_service, target.parsed, target_detection)
        interaction_differences = _interaction_differences(report)
        correspondences = tuple(_correspondence_from_snapshot(item) for item in report.analysis.correspondences)
        local_displacements = tuple(report.displacement_vectors)
        comparison = active.comparison_service.compare(
            reference_detection.candidates,
            target_detection.candidates,
            correspondences,
            transform=report.analysis.transform,
            reference_volumes=reference_volumes,
            target_volumes=target_volumes,
            mutations=report.analysis.mutations,
            interaction_differences=interaction_differences,
            local_displacements=local_displacements,
            reference_qc=report.input_quality.reference,
            target_qc=report.input_quality.target,
        )
        workflow_status = _aggregate_availability(
            reference_detection.availability,
            target_detection.availability,
            *(item.availability for item in reference_volumes.values()),
            *(item.availability for item in target_volumes.values()),
            comparison.availability,
        )
        payload = _comparison_payload(
            reference,
            target,
            report,
            reference_detection,
            target_detection,
            reference_volumes,
            target_volumes,
            comparison,
            workflow_status,
        )
        if json_path is not None:
            _write_json(json_path, payload)
        return WorkflowResult(
            workflow_status.value,
            payload,
            render_pocket_compare(
                comparison,
                status=workflow_status,
                reference_detection=reference_detection,
                target_detection=target_detection,
                reference_volumes=tuple(reference_volumes.values()),
                target_volumes=tuple(target_volumes.values()),
                additional_diagnostics=report.diagnostics,
            ),
            exit_code=_exit_code(workflow_status),
            snapshots=(reference.snapshot, target.snapshot),
        )
    except (StructLensError, OSError, ValueError) as error:
        return _failure(error)


def _single_pocket_payload(
    captured: CapturedInput,
    detection: PocketDetectionReport,
    volumes: tuple[PocketVolumeResult, ...],
    workflow_status: Availability,
) -> dict[str, Any]:
    diagnostics = tuple(detection.diagnostics) + tuple(item for volume in volumes for item in volume.diagnostics)
    return {
        "selection": selection_payload(captured.selection),
        "source": source_payload(captured.snapshot),
        "pockets": {
            "availability": workflow_status.value,
            "detection": detection.to_json(),
            "detections": [
                {
                    "role": "reference",
                    "availability": detection.availability.value,
                    "candidates": [candidate.to_json() for candidate in detection.candidates],
                    "diagnostics": [item.to_json() for item in detection.diagnostics],
                    "counts": dict(detection.counts),
                    "provenance": detection.provenance.to_json() if detection.provenance else None,
                }
            ],
            "volumes": [item.to_json() for item in volumes],
            "units": {"volume": "angstrom^3", "length": "angstrom"},
            "diagnostics": [item.to_json() for item in diagnostics],
            "provenance": detection.provenance.to_json() if detection.provenance else None,
        },
    }


def _comparison_payload(
    reference: CapturedInput,
    target: CapturedInput,
    report: Any,
    reference_detection: PocketDetectionReport,
    target_detection: PocketDetectionReport,
    reference_volumes: dict[str, PocketVolumeResult],
    target_volumes: dict[str, PocketVolumeResult],
    comparison: PocketComparisonReport,
    workflow_status: Availability,
) -> dict[str, Any]:
    diagnostics = (
        tuple(report.diagnostics)
        + tuple(reference_detection.diagnostics)
        + tuple(target_detection.diagnostics)
        + tuple(item for volume in reference_volumes.values() for item in volume.diagnostics)
        + tuple(item for volume in target_volumes.values() for item in volume.diagnostics)
        + tuple(comparison.diagnostics)
    )
    return {
        "reference_selection": selection_payload(reference.selection),
        "target_selection": selection_payload(target.selection),
        "reference_source": source_payload(reference.snapshot),
        "target_source": source_payload(target.snapshot),
        "analysis_report": report.to_json(),
        "analysis": report.analysis.to_json() if report.analysis is not None else None,
        "pockets": {
            "availability": workflow_status.value,
            "units": {"volume": "angstrom^3", "length": "angstrom"},
            "reference": {
                "detection": reference_detection.to_json(),
                "volumes": [item.to_json() for item in reference_volumes.values()],
            },
            "target": {
                "detection": target_detection.to_json(),
                "volumes": [item.to_json() for item in target_volumes.values()],
            },
            "matching": comparison.matching.to_json(),
            "comparisons": [item.to_json() for item in comparison.comparisons],
            "concordance": comparison.concordance.to_json() if comparison.concordance else None,
            "comparison_report": comparison.to_json(),
            "diagnostics": [item.to_json() for item in diagnostics],
        },
        "provenance": report.provenance.to_json() if report.provenance else None,
    }


def _measure_volumes(
    service: StructurePocketService,
    parsed: ParsedStructure,
    detection: PocketDetectionReport,
) -> tuple[PocketVolumeResult, ...]:
    return tuple(service.measure_volume(parsed, candidate) for candidate in detection.candidates)


def _volume_map(
    service: StructurePocketService,
    parsed: ParsedStructure,
    detection: PocketDetectionReport,
) -> dict[str, PocketVolumeResult]:
    return {candidate.candidate_id: service.measure_volume(parsed, candidate) for candidate in detection.candidates}


def _interaction_differences(report: Any) -> tuple[Any, ...]:
    if isinstance(report.interactions, InteractionEvidence):
        return tuple(report.interactions.differences)
    return tuple(item for item in (report.interactions or ()) if hasattr(item, "change"))


def _correspondence_from_snapshot(item: Any) -> ResidueCorrespondence:
    return ResidueCorrespondence(
        item.alignment_index,
        item.reference,
        item.target,
        item.reference_one_letter,
        item.target_one_letter,
        item.status,
        item.sequence_score,
        item.ca_displacement_angstrom,
        item.backbone_rmsd_angstrom,
        item.sidechain_rmsd_angstrom,
        item.all_heavy_atom_rmsd_angstrom,
        item.is_outlier,
        item.is_key_residue,
        item.mapping_source,
        item.mapping_locked,
    )


def _alignment_mode(value: str | AlignmentMode) -> AlignmentMode:
    try:
        return value if isinstance(value, AlignmentMode) else AlignmentMode(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"unknown alignment mode: {value!r}") from error


def _require_pairwise_chains(reference: ParsedStructure, target: ParsedStructure) -> None:
    if len(reference.protein_structure.chains) != 1 or len(target.protein_structure.chains) != 1:
        raise ValueError(
            "pairwise workflows require exactly one selected protein chain per input; "
            "repeat the chain selector with an unambiguous choice"
        )


def _validate_outputs(*paths: str | Path | None) -> None:
    resolved = [Path(path).resolve() for path in paths if path is not None]
    if len(set(resolved)) != len(resolved):
        raise ValueError("output paths must be distinct")


def _write_json(path: str | Path, payload: dict[str, Any]) -> None:
    atomic_write_bytes(path, json_bytes(payload), max_bytes=100 * 1024 * 1024, label="CLI JSON output", idempotent=True)


def _failure(error: Exception) -> WorkflowResult:
    return WorkflowResult("invalid_input", exit_code=2, error=str(error))


def _exit_code(availability: Availability) -> int:
    return 0 if availability in _SUCCESS_STATES else 2


def _aggregate_availability(*statuses: Availability) -> Availability:
    """Preserve fatal evidence states while keeping explicit no-results successful."""

    values = tuple(statuses)
    for fatal in (Availability.INVALID_INPUT, Availability.DEPENDENCY_UNAVAILABLE, Availability.NUMERICAL_FAILURE):
        if fatal in values:
            return fatal
    if Availability.NOT_DETECTED in values:
        return Availability.NOT_DETECTED
    if Availability.NOT_APPLICABLE in values:
        return Availability.NOT_APPLICABLE
    return Availability.AVAILABLE


_SUCCESS_STATES = frozenset(
    {
        Availability.AVAILABLE,
        Availability.NOT_DETECTED,
        Availability.NOT_APPLICABLE,
    }
)


__all__ = [
    "CliServices",
    "WorkflowResult",
    "run_compare",
    "run_pocket_compare",
    "run_pockets",
    "run_qc",
]
