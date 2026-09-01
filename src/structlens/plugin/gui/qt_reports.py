"""Deterministic projections from canonical reports for the Qt view."""

from __future__ import annotations

from structlens.application.chart_data import ChartDataset, ChartSeries, MatrixDataset
from structlens.core.evidence import SiteEvidence
from structlens.core.models import AnalysisResult, ResidueCorrespondence
from structlens.core.reports import AnalysisReport
from structlens.core.sites import SiteMetrics


def report_site_metrics(report: AnalysisReport) -> tuple[SiteMetrics, ...]:
    value = report.sites
    if isinstance(value, SiteEvidence):
        return tuple(value.metrics)
    return tuple(value or ())


def report_chart_datasets(report: AnalysisReport) -> dict[str, ChartDataset | MatrixDataset]:
    datasets: dict[str, ChartDataset | MatrixDataset] = {}
    if report.msa is not None:
        datasets["MSA conservation profile"] = ChartDataset(
            "msa_conservation_profile",
            "MSA conservation profile",
            "Alignment column",
            "Conservation",
            "fraction",
            (
                ChartSeries(
                    "Alignment conservation",
                    tuple((float(item.index + 1), item.conservation_score) for item in report.msa.columns),
                ),
            ),
            "Authoritative alignment conservation; gaps and ambiguous residues remain explicit.",
        )
    if report.analysis is not None:
        datasets["Structural deviation profile"] = ChartDataset(
            "structural_deviation_profile",
            "Structural deviation profile",
            "Reference residue position",
            "Cα displacement (Å)",
            "Å",
            (
                ChartSeries(
                    "Target",
                    tuple(
                        (float(item.reference.auth_seq_id), item.ca_displacement_angstrom)
                        for item in report.analysis.correspondences
                        if item.reference is not None and item.reference.auth_seq_id.lstrip("-").isdigit()
                    ),
                    tuple(
                        f"{item.reference.chain_id}:{item.reference.auth_seq_id}{item.reference.insertion_code or ''}"
                        for item in report.analysis.correspondences
                        if item.reference is not None and item.reference.auth_seq_id.lstrip("-").isdigit()
                    ),
                    tuple(
                        {"alignment_index": item.alignment_index}
                        for item in report.analysis.correspondences
                        if item.reference is not None and item.reference.auth_seq_id.lstrip("-").isdigit()
                    ),
                ),
            ),
            "Each point is an explicit mapped residue metric; missing values remain empty.",
        )
    return datasets


def legacy_analysis(report: AnalysisReport | None) -> AnalysisResult | None:
    """Adapt one report snapshot for the compatibility PyMOL renderer."""

    if report is None or report.analysis is None:
        return None
    snapshot = report.analysis
    correspondences = tuple(
        ResidueCorrespondence(
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
        for item in snapshot.correspondences
    )
    return AnalysisResult(
        reference_id=snapshot.reference_id,
        target_id=snapshot.target_id,
        correspondences=correspondences,
        mutations=snapshot.mutations,
        sequence_identity=snapshot.sequence_identity,
        sequence_coverage=snapshot.sequence_coverage,
        alignment_decision=snapshot.alignment_decision,
        sequence_similarity=snapshot.sequence_similarity,
        strict_rmsd_angstrom=snapshot.strict_rmsd_angstrom,
        refined_rmsd_angstrom=snapshot.refined_rmsd_angstrom,
        mapped_residue_count=snapshot.mapped_residue_count,
        refined_residue_count=snapshot.refined_residue_count,
        excluded_alignment_indices=snapshot.excluded_alignment_indices,
        tm_score=snapshot.tm_score,
        provenance=dict(snapshot.legacy_provenance),
        transform=snapshot.transform,
        method_provenance=snapshot.method_provenance,
    )


__all__ = ["legacy_analysis", "report_chart_datasets", "report_site_metrics"]
