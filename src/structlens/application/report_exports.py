"""Canonical exports for typed ``AnalysisReport`` values."""

from __future__ import annotations

import csv
import io
import json
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

from openpyxl import Workbook
from openpyxl.styles import Font

from structlens.application.report_pocket_rows import pocket_export_rows
from structlens.application.report_serialization import serialize_report
from structlens.application.report_xlsx import append_safe, deterministic_xlsx, normalise_workbook
from structlens.core.evidence import InteractionEvidence, SiteEvidence
from structlens.core.interactions import InteractionDifference, InteractionRecord
from structlens.core.reports import AnalysisReport
from structlens.core.reports.safe_io import atomic_write_bytes
from structlens.core.reports.serialization import interaction_json, mutation_json, residue_json
from structlens.core.sites import SiteMetrics

_MAX_EXPORT_BYTES = 100 * 1024 * 1024
_MAX_EXPORT_ROWS = 100_000
_MAX_EXPORT_CELLS = 1_000_000
_MAX_CELL_STRING = 32_767
_FIELDS = ["section", "role", "metric", "value", "units", "status", "reason", "provenance"]


def _safe(value: Any) -> Any:
    if isinstance(value, str):
        if len(value) > _MAX_CELL_STRING:
            raise ValueError("export cell string exceeds the Excel cell length limit")
        safe_value = "'" + value if value.lstrip().startswith(("=", "+", "-", "@")) else value
        if len(safe_value) > _MAX_CELL_STRING:
            raise ValueError("formula-safe export cell string exceeds the Excel cell length limit")
        return safe_value
    return value


def _reason(status: str) -> str:
    return "" if status == "available" else status.replace("_", " ")


def _report_diagnostic_status(report: AnalysisReport) -> str:
    """Return the most informative native state for report-level diagnostics."""

    states = tuple(getattr(report.availability, name).value for name in report.availability.__dataclass_fields__)
    for status in ("invalid_input", "numerical_failure", "dependency_unavailable", "not_detected", "not_applicable"):
        if status in states:
            return status
    return "available"


def _json(value: Any) -> str:
    serialized = json.dumps(value, ensure_ascii=False, sort_keys=True, default=str)
    if len(serialized) > _MAX_CELL_STRING:
        raise ValueError("serialized export cell exceeds the Excel cell length limit")
    return serialized


def _provenance_hash(hashes: Mapping[str, str], role: str, kind: str) -> str | None:
    for key, value in hashes.items():
        normalized = key.casefold().replace("-", "_")
        if normalized in {f"{role}_{kind}", f"{role}_{kind}_sha256", f"{role}_source_{kind}"}:
            return value
    return None


def _append(
    rows: list[dict[str, Any]],
    section: str,
    role: str,
    metric: str,
    value: Any,
    units: str,
    status: str,
    *,
    provenance: Any = "",
) -> None:
    if isinstance(value, (Mapping, list, tuple)):
        value = _json(value)
    if isinstance(provenance, (Mapping, list, tuple)):
        provenance = _json(provenance)
    if len(rows) >= _MAX_EXPORT_ROWS:
        raise ValueError("export exceeds the maximum row count")
    if (len(rows) + 1) * len(_FIELDS) > _MAX_EXPORT_CELLS:
        raise ValueError("export exceeds the maximum cell count")
    rows.append(
        {
            "section": section,
            "role": role,
            "metric": metric,
            "value": value,
            "units": units,
            "status": status,
            "reason": _reason(status),
            "provenance": provenance,
        }
    )


def _rows(report: AnalysisReport) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    availability = report.availability
    analysis = report.analysis
    for summary_metric, value, units in (
        ("report_id", report.report_id, ""),
        ("sequence_identity", analysis.sequence_identity if analysis else None, "fraction"),
        ("sequence_coverage", analysis.sequence_coverage if analysis else None, "fraction"),
        ("strict_rmsd_angstrom", analysis.strict_rmsd_angstrom if analysis else None, "angstrom"),
        ("refined_rmsd_angstrom", analysis.refined_rmsd_angstrom if analysis else None, "angstrom"),
    ):
        status = "available" if summary_metric == "report_id" else availability.analysis.value
        _append(rows, "summary", "report", summary_metric, value, units, status)
    hashes = {} if report.provenance is None else report.provenance.input_hashes
    for role, selection in (("reference", report.reference_selection), ("target", report.target_selection)):
        _append(rows, "structures", role, "content_id", selection.content_id, "sha256", "available")
        _append(rows, "structures", role, "raw_sha256", _provenance_hash(hashes, role, "raw"), "sha256", "available")
        _append(rows, "structures", role, "format", selection.format.value, "", "available")
        _append(rows, "structures", role, "model_id", selection.model_id, "", "available")
        _append(rows, "structures", role, "selection_id", selection.selection_id, "sha256", "available")
        _append(rows, "structures", role, "display_name", selection.display_name, "", "available")
        _append(rows, "structures", role, "author_chain_ids", _json(selection.author_chain_ids), "", "available")
        _append(rows, "structures", role, "label_chain_ids", _json(selection.label_chain_ids), "", "available")
        _append(rows, "structures", role, "altloc_policy", selection.altloc_policy.value, "", "available")
        _append(rows, "structures", role, "assembly_scope", selection.assembly_scope.value, "", "available")
    for role, quality in (("reference", report.input_quality.reference), ("target", report.input_quality.target)):
        status = quality.availability.value
        _append(rows, "qc", role, "availability", status, "", status)
        _append(rows, "qc", role, "error_count", quality.error_count, "count", status)
        _append(rows, "qc", role, "warning_count", quality.warning_count, "count", status)
        for diagnostic in quality.diagnostics:
            _append(
                rows,
                "diagnostics",
                role,
                diagnostic.code,
                diagnostic.message,
                "",
                status,
                provenance=diagnostic.source_id or "",
            )
        if quality.provenance is not None:
            _append(rows, "provenance", role, "qc_artifact_id", quality.provenance.artifact_id, "sha256", status)
    site_values: tuple[SiteMetrics, ...]
    if isinstance(report.sites, SiteEvidence):
        site_values = report.sites.metrics
    elif report.sites is None:
        site_values = ()
    else:
        site_values = tuple(report.sites)
    for site_metric in site_values:
        status = availability.sites.value
        _append(rows, "site_metrics", site_metric.structure_id, "site_id", site_metric.site_id, "", status)
        _append(
            rows,
            "site_metrics",
            site_metric.structure_id,
            "mapped_residue_count",
            site_metric.mapped_residue_count,
            "count",
            status,
        )
        _append(
            rows,
            "site_metrics",
            site_metric.structure_id,
            "coverage_fraction",
            site_metric.coverage_fraction,
            "fraction",
            status,
        )
        _append(
            rows,
            "site_metrics",
            site_metric.structure_id,
            "atomic_envelope_volume",
            site_metric.atomic_envelope_volume_angstrom3,
            "angstrom^3",
            status,
        )
        _append(
            rows, "site_metrics", site_metric.structure_id, "sasa", site_metric.sasa_angstrom2, "angstrom^2", status
        )
    pockets = report.pockets
    if pockets is not None:
        for detection in pockets.detections:
            status = detection.availability.value
            _append(rows, "pockets", detection.role, "detection_evidence", detection.to_json(), "", status)
            _append(rows, "pockets", detection.role, "candidate_count", len(detection.candidates), "count", status)
            for candidate in detection.candidates:
                _append(rows, "pockets", detection.role, "candidate_evidence", candidate.to_json(), "", status)
                _append(
                    rows,
                    "pockets",
                    detection.role,
                    "candidate_id",
                    candidate.candidate_id,
                    "sha256",
                    status,
                    provenance=candidate.source_content_id or "",
                )
                _append(
                    rows,
                    "pockets",
                    detection.role,
                    "lining_residue_count",
                    len(candidate.lining_residues),
                    "count",
                    status,
                    provenance=candidate.selection_id or "",
                )
            for key, value in sorted(detection.counts.items()):
                _append(rows, "pockets", detection.role, key, value, "count", status)
            if detection.settings is not None:
                _append(
                    rows,
                    "methods",
                    detection.role,
                    "pocket_detection_settings",
                    _json(detection.settings.to_json()),
                    "",
                    status,
                )
            if detection.provenance is not None:
                _append(
                    rows,
                    "provenance",
                    detection.role,
                    "pocket_detection_artifact_id",
                    detection.provenance.artifact_id,
                    "sha256",
                    status,
                )
            for diagnostic in detection.diagnostics:
                _append(rows, "diagnostics", detection.role, diagnostic.code, diagnostic.message, "", status)
        for volume in pockets.volumes:
            status = volume.availability.value if volume.availability is not None else "not_applicable"
            result = volume.result
            _append(rows, "pockets", volume.role, "volume_evidence", volume.to_json(), "", status)
            _append(
                rows,
                "pockets",
                volume.role,
                "free_volume",
                result.fine_volume_angstrom3 if result else None,
                "angstrom^3",
                status,
            )
            _append(
                rows,
                "pockets",
                volume.role,
                "coarse_free_volume",
                result.coarse_volume_angstrom3 if result else None,
                "angstrom^3",
                status,
            )
            _append(
                rows,
                "pockets",
                volume.role,
                "volume_sensitivity",
                result.absolute_sensitivity_angstrom3 if result else None,
                "angstrom^3",
                status,
            )
            _append(
                rows,
                "pockets",
                volume.role,
                "volume_candidate_id",
                result.candidate_id if result else None,
                "sha256",
                status,
            )
            if result is not None:
                if result.settings is not None:
                    _append(
                        rows,
                        "methods",
                        volume.role,
                        "pocket_volume_settings",
                        _json(result.settings.to_json()),
                        "",
                        status,
                    )
                if result.provenance is not None:
                    _append(
                        rows,
                        "provenance",
                        volume.role,
                        "pocket_volume_artifact_id",
                        result.provenance.artifact_id,
                        "sha256",
                        status,
                    )
                for diagnostic in result.diagnostics:
                    _append(rows, "diagnostics", volume.role, diagnostic.code, diagnostic.message, "", status)
        if pockets.matching is not None:
            status = pockets.matching.availability.value
            matching_payload = pockets.matching.to_json()
            matching_rows: tuple[tuple[str, object, str], ...] = (
                ("matching_evidence", matching_payload, ""),
                ("matching_settings", matching_payload["settings"], ""),
                ("matching_transform", matching_payload["transform"], ""),
                ("matching_diagnostics", matching_payload["diagnostics"], ""),
                ("matching_optimal_score", pockets.matching.optimal_score, "fraction"),
                ("matching_second_best_score", pockets.matching.second_best_score, "fraction"),
            )
            for matching_metric, matching_value, matching_units in matching_rows:
                _append(rows, "pocket_matches", "matching", matching_metric, matching_value, matching_units, status)
            for match in pockets.matching.matches:
                assignment_score = (
                    match.assignment_score if match.assignment_score is not None else pockets.matching.optimal_score
                )
                alternative_assignment_score = (
                    match.alternative_assignment_score
                    if match.alternative_assignment_score is not None
                    else pockets.matching.second_best_score
                )
                for metric, value, units in (
                    (
                        "reference_candidate_id",
                        match.reference_candidate.candidate_id if match.reference_candidate else None,
                        "sha256",
                    ),
                    (
                        "target_candidate_id",
                        match.target_candidate.candidate_id if match.target_candidate else None,
                        "sha256",
                    ),
                    ("match_state", match.state.value, ""),
                    ("score", match.score, match.units.get("score", "fraction")),
                    ("alternative_score", match.alternative_score, match.units.get("score", "fraction")),
                    ("ambiguity_margin", match.ambiguity_margin, "fraction"),
                    ("assignment_score", assignment_score, "fraction"),
                    ("alternative_assignment_score", alternative_assignment_score, "fraction"),
                    ("lining_jaccard", match.lining_jaccard, match.units.get("lining_jaccard", "fraction")),
                    (
                        "centroid_distance",
                        match.centroid_distance_angstrom,
                        match.units.get("centroid_distance_angstrom", "angstrom"),
                    ),
                ):
                    _append(rows, "pocket_matches", "comparison", metric, value, units, status)
                _append(rows, "pocket_matches", "comparison", "units", _json(match.units), "", status)
                for diagnostic in match.diagnostics:
                    _append(rows, "diagnostics", "comparison", diagnostic.code, diagnostic.message, "", status)
            for diagnostic in pockets.matching.diagnostics:
                _append(rows, "diagnostics", "matching", diagnostic.code, diagnostic.message, "", status)
        for comparison in pockets.comparisons:
            status = comparison.availability.value
            match = comparison.match
            comparison_role = "{}->{}".format(
                match.reference_candidate.candidate_id if match.reference_candidate else "",
                match.target_candidate.candidate_id if match.target_candidate else "",
            )
            _append(
                rows,
                "pocket_comparisons",
                comparison_role,
                "comparison_evidence",
                comparison.to_json(),
                "",
                status,
            )
            for metric, value, units in (
                (
                    "reference_candidate_id",
                    match.reference_candidate.candidate_id if match.reference_candidate else None,
                    "sha256",
                ),
                (
                    "target_candidate_id",
                    match.target_candidate.candidate_id if match.target_candidate else None,
                    "sha256",
                ),
                ("reference_volume", comparison.reference_volume_angstrom3, "angstrom^3"),
                ("target_volume", comparison.target_volume_angstrom3, "angstrom^3"),
                ("volume_delta", comparison.volume_delta_angstrom3, "angstrom^3"),
                ("relative_volume_delta", comparison.relative_volume_delta_fraction, "fraction"),
                ("surface_delta", comparison.surface_delta_angstrom2, "angstrom^2"),
                ("relative_surface_delta", comparison.relative_surface_delta_fraction, "fraction"),
                ("reference_surface_area", comparison.reference_surface_area_angstrom2, "angstrom^2"),
                ("target_surface_area", comparison.target_surface_area_angstrom2, "angstrom^2"),
                ("reference_surface_method", comparison.reference_surface_method, ""),
                ("target_surface_method", comparison.target_surface_method, ""),
                ("local_displacement", comparison.local_displacement_angstrom, "angstrom"),
                ("ca_displacement", comparison.ca_displacement_angstrom, "angstrom"),
                ("sidechain_displacement", comparison.sidechain_displacement_angstrom, "angstrom"),
                ("qc_availability", comparison.qc_availability.value if comparison.qc_availability else None, ""),
                ("conserved_lining", _json([residue_json(item) for item in comparison.lining_residue_conserved]), ""),
                ("gained_lining", _json([residue_json(item) for item in comparison.lining_residue_gains]), ""),
                ("lost_lining", _json([residue_json(item) for item in comparison.lining_residue_losses]), ""),
                ("associated_mutations", _json([mutation_json(item) for item in comparison.associated_mutations]), ""),
                ("interaction_changes", _json([interaction_json(item) for item in comparison.interaction_changes]), ""),
                (
                    "local_displacements",
                    _json(
                        [
                            {"residue": residue_json(residue), "magnitude_angstrom": magnitude}
                            for residue, magnitude in comparison.local_displacements
                        ]
                    ),
                    "angstrom",
                ),
                ("units", _json(comparison.units), ""),
                (
                    "reference_surface_units",
                    _json(comparison.reference_surface_units or {}),
                    "",
                ),
                ("target_surface_units", _json(comparison.target_surface_units or {}), ""),
                (
                    "reference_surface_provenance",
                    _json(comparison.reference_surface_provenance or {}),
                    "",
                ),
                ("target_surface_provenance", _json(comparison.target_surface_provenance or {}), ""),
            ):
                _append(rows, "pocket_comparisons", comparison_role, metric, value, units, status)
            if comparison.volume is not None:
                comparison_volume = comparison.volume
                for metric, value, units in (
                    (
                        "reference_sensitivity",
                        (
                            _json(comparison_volume.reference_sensitivity.to_json())
                            if comparison_volume.reference_sensitivity
                            else None
                        ),
                        "",
                    ),
                    (
                        "target_sensitivity",
                        (
                            _json(comparison_volume.target_sensitivity.to_json())
                            if comparison_volume.target_sensitivity
                            else None
                        ),
                        "",
                    ),
                    ("volume_units", _json(comparison_volume.units), ""),
                    ("volume_provenance", _json([item.to_json() for item in comparison_volume.provenance]), ""),
                ):
                    _append(rows, "pocket_comparisons", comparison_role, metric, value, units, status)
            for diagnostic in comparison.qc_diagnostics:
                qc_status = comparison.qc_availability.value if comparison.qc_availability else status
                _append(rows, "diagnostics", comparison_role, diagnostic.code, diagnostic.message, "", qc_status)
            for diagnostic in comparison.diagnostics:
                _append(rows, "diagnostics", comparison_role, diagnostic.code, diagnostic.message, "", status)
        for lining in pockets.lining_residues:
            status = lining.availability.value
            _append(rows, "lining_residues", lining.role, "lining_evidence", lining.to_json(), "", status)
            _append(rows, "lining_residues", lining.role, lining.candidate_id, len(lining.residues), "count", status)
            for residue in lining.residues:
                _append(
                    rows,
                    "lining_residues",
                    lining.role,
                    "residue",
                    residue_json(residue),
                    "",
                    status,
                    provenance=lining.candidate_id,
                )
            for diagnostic in lining.diagnostics:
                _append(rows, "diagnostics", lining.role, diagnostic.code, diagnostic.message, "", status)
        for diagnostic in pockets.diagnostics:
            pocket_status = pockets.availability.value if pockets.availability is not None else "not_applicable"
            _append(rows, "diagnostics", "report", diagnostic.code, diagnostic.message, "", pocket_status)
        if pockets.provenance is not None:
            pocket_status = pockets.availability.value if pockets.availability is not None else "not_applicable"
            _append(rows, "provenance", "pockets", "pocket_method_id", pockets.provenance.method_id, "", pocket_status)
            _append(
                rows,
                "provenance",
                "pockets",
                "pocket_method_version",
                pockets.provenance.method_version,
                "",
                pocket_status,
            )
            _append(
                rows,
                "provenance",
                "pockets",
                "pocket_analyzed_representation",
                pockets.provenance.analyzed_representation,
                "",
                pocket_status,
            )
            _append(
                rows,
                "provenance",
                "pockets",
                "pocket_parameters",
                _json(pockets.provenance.parameters),
                "",
                pocket_status,
            )
            _append(
                rows,
                "provenance",
                "pockets",
                "pocket_units",
                _json(pockets.provenance.units),
                "",
                pocket_status,
            )
            _append(
                rows,
                "provenance",
                "pockets",
                "pocket_backend_versions",
                _json(pockets.provenance.backend_versions),
                "",
                pocket_status,
            )
            for key, value in sorted(pockets.provenance.input_hashes.items()):
                _append(rows, "provenance", "pockets", f"pocket_input_{key}", value, "sha256", pocket_status)
            _append(
                rows,
                "provenance",
                "pockets",
                "pocket_artifact_id",
                pockets.provenance.artifact_id,
                "sha256",
                pocket_status,
            )
        for pocket_row in pocket_export_rows(pockets, serialize=_json):
            _append(
                rows,
                pocket_row.section,
                pocket_row.role,
                pocket_row.metric,
                pocket_row.value,
                pocket_row.units,
                pocket_row.status,
                provenance=pocket_row.provenance,
            )
    interactions: Iterable[Any] = ()
    if isinstance(report.interactions, InteractionEvidence):
        interactions = report.interactions.differences
    elif report.interactions is not None:
        interactions = report.interactions
    for interaction in interactions:
        if isinstance(interaction, InteractionDifference):
            _append(
                rows,
                "matches",
                "comparison",
                interaction.key.reference_position_a,
                interaction.change.value,
                "",
                availability.interactions.value,
            )
        elif isinstance(interaction, InteractionRecord):
            _append(
                rows,
                "matches",
                interaction.structure_id,
                interaction.interaction_type.value,
                interaction.distance_angstrom,
                "angstrom",
                availability.interactions.value,
            )
    for diagnostic in report.diagnostics:
        _append(
            rows,
            "diagnostics",
            "report",
            diagnostic.code,
            diagnostic.message,
            "",
            _report_diagnostic_status(report),
            provenance=diagnostic.source_id or "",
        )
    if report.provenance is not None:
        _append(rows, "methods", "report", "method_id", report.provenance.method_id, "", "available")
        _append(rows, "methods", "report", "method_version", report.provenance.method_version, "", "available")
        _append(
            rows,
            "methods",
            "report",
            "analyzed_representation",
            report.provenance.analyzed_representation,
            "",
            "available",
        )
        _append(rows, "methods", "report", "parameters", _json(report.provenance.parameters), "", "available")
        for key, value in sorted(report.provenance.input_hashes.items()):
            _append(rows, "provenance", "report", key, value, "sha256", "available")
        _append(rows, "provenance", "report", "artifact_id", report.provenance.artifact_id, "sha256", "available")
    _validate_rows(rows)
    return rows


def _row_width(row: Mapping[str, Any]) -> int:
    if row["section"] == "summary":
        return 5
    if row["section"] == "pockets":
        return 9
    if row["section"] == "provenance":
        return 2
    return len(_FIELDS)


def _validate_rows(rows: list[dict[str, Any]]) -> None:
    """Bound all data before allocating an output writer/workbook."""

    if len(rows) > _MAX_EXPORT_ROWS:
        raise ValueError("export exceeds the maximum row count")
    header_cells = 9 * len(_FIELDS)
    cell_count = header_cells + sum(_row_width(row) for row in rows)
    if cell_count > _MAX_EXPORT_CELLS:
        raise ValueError("export exceeds the maximum cell count")
    estimated_bytes = header_cells * 16
    for row in rows:
        for value in row.values():
            safe_value = _safe(value)
            if isinstance(safe_value, str):
                estimated_bytes += len(safe_value.encode("utf-8")) + 16
            else:
                estimated_bytes += 32
    if estimated_bytes > _MAX_EXPORT_BYTES:
        raise ValueError("estimated export exceeds the maximum serialized size")


def _checked(
    report: Any, snapshots: tuple[Any, ...], source_paths: tuple[str | Path, ...], snapshot_dir: str | Path | None
) -> AnalysisReport:
    if not isinstance(report, AnalysisReport):
        raise TypeError("report must be an AnalysisReport")
    serialize_report(report, snapshots=snapshots, source_paths=source_paths, snapshot_dir=snapshot_dir)
    return report


def export_report_json(
    report: Any,
    path: str | Path,
    *,
    snapshots: tuple[Any, ...] = (),
    source_paths: tuple[str | Path, ...] = (),
    snapshot_dir: str | Path | None = None,
) -> None:
    from structlens.core.reports.safe_io import atomic_write_bytes

    checked = _checked(report, snapshots, source_paths, snapshot_dir)
    atomic_write_bytes(
        path,
        serialize_report(checked, snapshots=snapshots, source_paths=source_paths, snapshot_dir=snapshot_dir),
        max_bytes=_MAX_EXPORT_BYTES,
        label="report export",
        idempotent=True,
    )


def export_report_table(
    report: Any,
    path: str | Path,
    *,
    delimiter: str = ",",
    snapshots: tuple[Any, ...] = (),
    source_paths: tuple[str | Path, ...] = (),
    snapshot_dir: str | Path | None = None,
) -> None:
    checked = _checked(report, snapshots, source_paths, snapshot_dir)
    rows = _rows(checked)
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=_FIELDS, delimiter=delimiter, extrasaction="raise", lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({key: _safe(row[key]) for key in _FIELDS})
        if output.tell() > _MAX_EXPORT_BYTES // 4:
            raise ValueError("delimited export exceeds the maximum serialized size")
    atomic_write_bytes(
        path,
        output.getvalue().encode("utf-8"),
        max_bytes=_MAX_EXPORT_BYTES,
        label="delimited report export",
        idempotent=True,
    )


def export_report_csv(report: Any, path: str | Path, **kwargs: Any) -> None:
    export_report_table(report, path, delimiter=",", **kwargs)


def export_report_tsv(report: Any, path: str | Path, **kwargs: Any) -> None:
    export_report_table(report, path, delimiter="\t", **kwargs)


def export_report_xlsx(
    report: Any,
    path: str | Path,
    *,
    snapshots: tuple[Any, ...] = (),
    source_paths: tuple[str | Path, ...] = (),
    snapshot_dir: str | Path | None = None,
) -> None:
    checked = _checked(report, snapshots, source_paths, snapshot_dir)
    rows = _rows(checked)
    workbook = Workbook()
    default = workbook.active
    if default is None:
        raise RuntimeError("Workbook did not create a worksheet")
    workbook.remove(default)
    by_section: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        by_section.setdefault(row["section"], []).append(row)
    for name in (
        "Summary",
        "Structures",
        "QC",
        "Site Metrics",
        "Pockets",
        "Matches",
        "Lining Residues",
        "Diagnostics",
        "Methods",
        "Provenance",
    ):
        sheet = workbook.create_sheet(name)
        headers = _FIELDS
        if name == "Summary":
            headers = ["Metric", "Value", "Units", "Status", "Reason"]
        if name == "Pockets":
            headers = [
                "Section",
                "Role",
                "Metric",
                "Value",
                "Units",
                "Free volume (angstrom^3)",
                "Status",
                "Reason",
                "Provenance",
            ]
        if name == "Provenance":
            headers = ["Key", "Value"]
        append_safe(sheet, headers, _safe)
        section = {
            "Summary": "summary",
            "Structures": "structures",
            "QC": "qc",
            "Site Metrics": "site_metrics",
            "Pockets": "pockets",
            "Matches": "matches",
            "Lining Residues": "lining_residues",
            "Diagnostics": "diagnostics",
            "Methods": "methods",
            "Provenance": "provenance",
        }[name]
        sections: tuple[str, ...] = (section,)
        if name == "Matches":
            sections = ("matches", "pocket_matches", "pocket_comparisons")
        selected_rows = [item for selected in sections for item in by_section.get(selected, [])]
        if name == "Matches":
            selected_rows.extend(item for item in by_section.get("diagnostics", []) if item["role"] == "matching")
        for row in selected_rows:
            values = [row[field] for field in _FIELDS]
            if name == "Summary":
                values = [row["metric"], row["value"], row["units"], row["status"], row["reason"]]
            if name == "Pockets":
                values.insert(5, row["value"] if row["metric"] == "free_volume" else None)
            if name == "Provenance":
                values = [row["metric"], row["value"]]
            append_safe(sheet, values, _safe)
        sheet.freeze_panes = "A2"
        sheet.auto_filter.ref = sheet.dimensions
        for cell in sheet[1]:
            cell.font = Font(bold=True)
    normalise_workbook(workbook)
    stream = io.BytesIO()
    workbook.save(stream)
    atomic_write_bytes(
        path,
        deterministic_xlsx(stream.getvalue()),
        max_bytes=_MAX_EXPORT_BYTES,
        label="XLSX report export",
        idempotent=True,
    )


__all__ = ["export_report_csv", "export_report_json", "export_report_table", "export_report_tsv", "export_report_xlsx"]
