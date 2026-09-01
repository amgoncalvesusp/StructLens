"""Adversarial real-data fixtures for Task 13 final RED regressions."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from structlens.application.dto import AnalysisReportRequest
from structlens.application.report_service import ReportService
from structlens.core.evidence import Availability, Diagnostic, DiagnosticSeverity, InteractionEvidence
from structlens.core.interactions import InteractionThresholds
from structlens.core.models import AlignmentMode, AnalysisSettings
from structlens.core.msa import MSASettings
from structlens.core.parsing import capture_snapshot
from structlens.core.reports import AnalysisReport

from ._task13_round4_fixtures import multi_model_chain_pdb, rich_report, selection


def typed_report(tmp_path: Path) -> tuple[object, AnalysisReport]:
    """Return a real report carrying every typed value the GUI must expose."""

    request, report = rich_report(tmp_path)
    assert isinstance(report.interactions, InteractionEvidence)
    target = report.interactions.target_interactions[0]
    reference = report.interactions.reference_interactions[0] if report.interactions.reference_interactions else target
    target = replace(target, angle_degrees=121.5, distance_angstrom=2.71)
    reference = replace(reference, angle_degrees=119.0, distance_angstrom=2.82)
    differences = tuple(
        replace(item, reference_record=reference, target_record=target) for item in report.interactions.differences
    )
    interactions = replace(
        report.interactions,
        differences=differences,
        reference_interactions=(reference,),
        target_interactions=(target,),
    )
    assert report.sites
    site = replace(
        report.sites[0],
        global_frame_backbone_rmsd_angstrom=1.20,
        site_fitted_backbone_rmsd_angstrom=0.80,
        centroid_displacement_angstrom=0.50,
        radius_of_gyration_angstrom=2.30,
        atomic_envelope_volume_angstrom3=12.50,
        sasa_angstrom2=42.00,
        polar_residue_fraction=0.75,
        charged_residue_fraction=0.25,
    )
    return request, replace(report, interactions=interactions, sites=(site,))


def unavailable_typed_report(tmp_path: Path) -> tuple[object, AnalysisReport]:
    """Return native failure states with their public diagnostic reasons."""

    request, report = typed_report(tmp_path)
    diagnostics = report.diagnostics + (
        Diagnostic(
            "report.interactions.failed",
            DiagnosticSeverity.ERROR,
            "Interaction geometry could not be calculated.",
            remediation="Check coordinate completeness and geometry settings.",
        ),
        Diagnostic(
            "report.sites.invalid",
            DiagnosticSeverity.ERROR,
            "The selected site input is invalid.",
            remediation="Review the selected site residues.",
        ),
    )
    availability = replace(
        report.availability,
        interactions=Availability.NUMERICAL_FAILURE,
        sites=Availability.INVALID_INPUT,
    )
    return request, replace(report, interactions=None, sites=None, availability=availability, diagnostics=diagnostics)


def gapped_deviation_report(tmp_path: Path) -> tuple[object, AnalysisReport]:
    """Return non-consecutive and insertion-coded authoritative residues."""

    request, report = rich_report(tmp_path)
    assert report.analysis is not None
    updated = []
    labels = ("10", "21", "42")
    insertion_codes = (None, "A", None)
    for index, item in enumerate(report.analysis.correspondences):
        if index >= len(labels):
            updated.append(item)
            continue
        reference = (
            replace(
                item.reference,
                auth_seq_id=labels[index],
                insertion_code=insertion_codes[index],
            )
            if item.reference is not None
            else None
        )
        target = (
            replace(
                item.target,
                auth_seq_id=str(110 + index),
                insertion_code=None,
            )
            if item.target is not None
            else None
        )
        updated.append(replace(item, alignment_index=(2, 7, 11)[index], reference=reference, target=target))
    analysis = replace(report.analysis, correspondences=tuple(updated))
    return request, replace(report, analysis=analysis)


def mismatch_report(tmp_path: Path) -> tuple[object, AnalysisReport]:
    """Return a valid report whose request identity is intentionally different."""

    return rich_report(tmp_path)


def snapshot_identity_report(tmp_path: Path) -> tuple[object, AnalysisReport]:
    """Alias documenting tests that use the exact captured report snapshots."""

    return rich_report(tmp_path)


def multi_selected_report(tmp_path: Path) -> tuple[AnalysisReportRequest, AnalysisReport]:
    """Analyze only model 2 / chain B from two identical multi-model files."""

    reference_path = tmp_path / "multi-reference.pdb"
    target_path = tmp_path / "multi-target.pdb"
    payload = multi_model_chain_pdb()
    reference_path.write_bytes(payload)
    target_path.write_bytes(payload)
    reference_snapshot = capture_snapshot(reference_path)
    target_snapshot = capture_snapshot(target_path)
    reference_selection = selection(reference_snapshot, model="2", chain="B", path=reference_path)
    # Keep the request asymmetric so PyMOL state isolation is observable:
    # reference MODEL 2 is loaded as state 2, while target MODEL 1 is state 1.
    target_selection = selection(target_snapshot, model="1", chain="B", path=target_path)
    request = AnalysisReportRequest(
        reference_snapshot,
        target_snapshot,
        reference_selection,
        target_selection,
        analysis_settings=AnalysisSettings(alignment_mode=AlignmentMode.SEQUENCE),
        msa_settings=MSASettings(),
        interaction_thresholds=InteractionThresholds(),
    )
    return request, ReportService().analyze(request)


__all__ = [
    "gapped_deviation_report",
    "mismatch_report",
    "multi_selected_report",
    "snapshot_identity_report",
    "typed_report",
    "unavailable_typed_report",
]
