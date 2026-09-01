"""Pure deterministic English rendering for command results."""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from typing import Any

from structlens.application.pocket_comparison_service import PocketComparisonReport
from structlens.application.pocket_service import PocketDetectionReport
from structlens.core.evidence import Availability, Diagnostic
from structlens.core.pockets import PocketVolumeResult
from structlens.core.quality import StructureQualityReport
from structlens.core.reports import AnalysisReport
from structlens.core.reports.serialization import selection_json


def json_bytes(payload: Mapping[str, Any]) -> bytes:
    """Encode a JSON envelope with stable ordering and no non-finite values."""

    return json.dumps(payload, ensure_ascii=False, allow_nan=False, sort_keys=True, separators=(",", ":")).encode(
        "utf-8"
    )


def selection_payload(selection: Any) -> dict[str, Any]:
    """Serialize a typed selection without exposing mutable local paths."""

    return selection_json(selection)


def source_payload(snapshot: Any) -> dict[str, Any]:
    """Serialize the immutable source identity used by a workflow."""

    return {
        "display_name": snapshot.display_name,
        "logical_format": snapshot.logical_format,
        "is_gzip": snapshot.is_gzip,
        "raw_sha256": snapshot.raw_sha256,
        "content_id": snapshot.content_id,
        "decompressed_sha256": snapshot.decompressed_sha256,
    }


def render_quality(report: StructureQualityReport) -> str:
    lines = [f"Status: {report.availability.value}"]
    for key, value in sorted(report.counts.items()):
        lines.append(f"{key}: {value}")
    lines.extend(_diagnostic_lines(report.diagnostics))
    return "\n".join(lines)


def render_pockets(
    report: PocketDetectionReport,
    *,
    volumes: Sequence[PocketVolumeResult] = (),
    status: Availability | None = None,
) -> str:
    lines = [
        f"Status: {(status or report.availability).value}",
        f"Detection: {report.availability.value}",
        f"Candidates: {report.candidate_count}",
        f"Volumes: {len(volumes)}",
    ]
    for index, candidate in enumerate(report.candidates, start=1):
        lines.append(f"Candidate {index}: {candidate.candidate_id}")
    for index, volume in enumerate(volumes, start=1):
        candidate_id = volume.candidate_id or "-"
        lines.append(f"Volume {index} [{candidate_id}]: {volume.availability.value}")
    diagnostics = tuple(report.diagnostics) + tuple(item for volume in volumes for item in volume.diagnostics)
    lines.extend(_diagnostic_lines(diagnostics))
    return "\n".join(lines)


def render_compare(report: AnalysisReport) -> str:
    analysis = report.analysis
    lines = [f"Availability: {report.availability.analysis.value}"]
    if analysis is None:
        lines.append("Analysis: unavailable")
        lines.extend(_diagnostic_lines(report.diagnostics))
        return "\n".join(lines)
    lines.extend(
        (
            f"Reference: {analysis.reference_id}",
            f"Target: {analysis.target_id}",
            f"Sequence identity: {analysis.sequence_identity:.3f}",
            f"Sequence coverage: {analysis.sequence_coverage:.3f}",
            f"Strict Cα RMSD: {_metric(analysis.strict_rmsd_angstrom)} Å",
            f"Refined Cα RMSD: {_metric(analysis.refined_rmsd_angstrom)} Å",
            f"Mapped residues: {analysis.mapped_residue_count}",
            f"Mutations: {len(analysis.mutations)}",
            f"Alignment: {analysis.alignment_decision}",
        )
    )
    lines.extend(_diagnostic_lines(report.diagnostics))
    return "\n".join(lines)


def render_pocket_compare(
    report: PocketComparisonReport,
    *,
    status: Availability | None = None,
    reference_detection: PocketDetectionReport | None = None,
    target_detection: PocketDetectionReport | None = None,
    reference_volumes: Sequence[PocketVolumeResult] = (),
    target_volumes: Sequence[PocketVolumeResult] = (),
    additional_diagnostics: Sequence[Diagnostic] = (),
) -> str:
    lines = [f"Status: {(status or report.availability).value}"]
    if reference_detection is not None:
        lines.append(
            f"Reference detection: {reference_detection.availability.value} ({reference_detection.candidate_count})"
        )
    if target_detection is not None:
        lines.append(f"Target detection: {target_detection.availability.value} ({target_detection.candidate_count})")
    if reference_volumes:
        lines.append("Reference volumes: " + ", ".join(item.availability.value for item in reference_volumes))
    if target_volumes:
        lines.append("Target volumes: " + ", ".join(item.availability.value for item in target_volumes))
    for index, match in enumerate(report.matches, start=1):
        reference = match.reference_candidate.candidate_id if match.reference_candidate else "-"
        target = match.target_candidate.candidate_id if match.target_candidate else "-"
        lines.append(f"Match {index}: {match.state.value} ({reference} -> {target})")
        if match.score is not None:
            lines.append(f"  Score: {match.score:.3f}")
        if index <= len(report.comparisons):
            comparison = report.comparisons[index - 1]
            if comparison.volume_delta_angstrom3 is not None:
                lines.append(f"  Volume delta: {comparison.volume_delta_angstrom3:.3f} Å³")
            if comparison.surface_delta_angstrom2 is not None:
                lines.append(f"  Surface delta: {comparison.surface_delta_angstrom2:.3f} Å²")
            if comparison.ca_displacement_angstrom is not None:
                lines.append(f"  Cα displacement: {comparison.ca_displacement_angstrom:.3f} Å")
            if comparison.local_displacement_angstrom is not None:
                lines.append(f"  Local displacement: {comparison.local_displacement_angstrom:.3f} Å")
    diagnostics = list(report.diagnostics) + list(additional_diagnostics)
    if reference_detection is not None:
        diagnostics.extend(reference_detection.diagnostics)
    if target_detection is not None:
        diagnostics.extend(target_detection.diagnostics)
    diagnostics.extend(item for volume in reference_volumes for item in volume.diagnostics)
    diagnostics.extend(item for volume in target_volumes for item in volume.diagnostics)
    lines.extend(_diagnostic_lines(tuple(diagnostics)))
    return "\n".join(lines)


def _metric(value: float | None) -> str:
    return "unavailable" if value is None else f"{value:.3f}"


def _diagnostic_lines(diagnostics: Sequence[Diagnostic]) -> list[str]:
    return [f"Diagnostic [{item.severity.value}] {item.code}: {item.message}" for item in diagnostics]


__all__ = [
    "json_bytes",
    "render_compare",
    "render_pocket_compare",
    "render_pockets",
    "render_quality",
    "selection_payload",
    "source_payload",
]
