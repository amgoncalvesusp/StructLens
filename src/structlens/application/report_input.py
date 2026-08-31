"""Input normalization helpers for the canonical pairwise report workflow.

The helpers in this module only normalize immutable input values and describe
parser-boundary failures.  They deliberately do not perform scientific
analysis or access application services.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import replace

from structlens.application.dto import AnalysisReportRequest
from structlens.core.evidence import Availability, Diagnostic, DiagnosticSeverity
from structlens.core.models import ProteinChain, ResidueId
from structlens.core.parsing import InputSelection, ParsedStructure, SourceSnapshot
from structlens.core.quality import StructureQualityReport
from structlens.core.sites import SiteDefinition


def inputs_are_usable(
    reference: ParsedStructure | None,
    target: ParsedStructure | None,
    reference_qc: StructureQualityReport,
    target_qc: StructureQualityReport,
) -> bool:
    """Return whether both inputs passed QC and resolve to one protein chain."""

    from structlens.core.evidence import Availability

    return (
        reference is not None
        and target is not None
        and reference_qc.availability is Availability.AVAILABLE
        and target_qc.availability is Availability.AVAILABLE
        and len(reference.protein_structure.chains) == 1
        and len(target.protein_structure.chains) == 1
    )


def selection_diagnostics(
    request: AnalysisReportRequest,
    reference: ParsedStructure | None,
    target: ParsedStructure | None,
) -> tuple[Diagnostic, ...]:
    """Describe selections that resolve to an unsupported number of chains."""

    diagnostics: list[Diagnostic] = []
    for role, parsed, selection in (
        ("reference", reference, request.reference_selection),
        ("target", target, request.target_selection),
    ):
        if parsed is not None and len(parsed.protein_structure.chains) != 1:
            diagnostics.append(
                Diagnostic(
                    code=f"report.input.{role}.chain_count",
                    severity=DiagnosticSeverity.ERROR,
                    message=(
                        f"The {role} selection resolved to {len(parsed.protein_structure.chains)} protein chains; "
                        "pairwise reporting requires exactly one."
                    ),
                    source_id=f"{role}.{selection.format.value}",
                    remediation="Select exactly one author or label chain for this comparison.",
                )
            )
    return tuple(diagnostics)


def apply_selection_diagnostics(
    reference: StructureQualityReport,
    target: StructureQualityReport,
    diagnostics: Sequence[Diagnostic],
) -> tuple[StructureQualityReport, StructureQualityReport]:
    """Project invalid pairwise selections into the corresponding QC reports."""

    updated: list[StructureQualityReport] = []
    for role, report in (("reference", reference), ("target", target)):
        role_diagnostics = tuple(
            item for item in diagnostics if item.code.startswith(f"report.input.{role}.")
        )
        if role_diagnostics and report.availability is Availability.AVAILABLE:
            report = replace(
                report,
                availability=Availability.INVALID_INPUT,
                diagnostics=report.diagnostics + role_diagnostics,
            )
        updated.append(report)
    return updated[0], updated[1]


def single_chain(parsed: ParsedStructure) -> ProteinChain:
    """Return the only selected chain, raising for an invalid pairwise input."""

    if len(parsed.protein_structure.chains) != 1:
        raise ValueError("the report workflow requires exactly one selected chain per structure")
    return parsed.protein_structure.chains[0]


def canonical_snapshot(snapshot: SourceSnapshot, role: str) -> SourceSnapshot:
    """Remove filesystem naming from a snapshot before parser invocation."""

    suffix = "pdb" if snapshot.logical_format == "pdb" else "cif"
    return SourceSnapshot(
        decompressed_bytes=snapshot.decompressed_bytes,
        raw_sha256=snapshot.raw_sha256,
        content_id=snapshot.content_id,
        decompressed_sha256=snapshot.decompressed_sha256,
        display_name=f"{role}.{suffix}",
        logical_format=snapshot.logical_format,
        is_gzip=snapshot.is_gzip,
    )


def canonical_selection(selection: InputSelection, role: str) -> InputSelection:
    """Normalize a selection display name while preserving its content identity."""

    suffix = "pdb" if selection.format.value == "pdb" else "cif"
    return InputSelection(
        selection.content_id,
        f"{role}.{suffix}",
        selection.format,
        selection.model_id,
        author_chain_ids=selection.author_chain_ids,
        label_chain_ids=selection.label_chain_ids,
        altloc_policy=selection.altloc_policy,
        assembly_scope=selection.assembly_scope,
        chain_locators=selection.chain_locators,
    )


def canonical_residue_id(residue: ResidueId, role: str = "reference") -> ResidueId:
    """Set an explicit structure role on a residue identity."""

    return ResidueId(
        role,
        residue.model_id,
        residue.chain_id,
        residue.auth_seq_id,
        residue.insertion_code,
        residue.residue_name,
    )


def canonical_optional_residue_id(
    residue: ResidueId | None,
    role: str = "reference",
) -> ResidueId | None:
    return None if residue is None else canonical_residue_id(residue, role)


def canonical_site_definition(definition: SiteDefinition) -> SiteDefinition:
    """Normalize all residue roles in a site definition."""

    return SiteDefinition(
        definition.site_id,
        definition.name,
        definition.mode,
        tuple(canonical_residue_id(item) for item in definition.reference_residues),
        canonical_optional_residue_id(definition.center_residue),
        definition.ligand_id,
        definition.radius_angstrom,
    )


def canonical_manual_pairs(
    pairs: Sequence[tuple[ResidueId, ResidueId]],
) -> list[tuple[object, object]]:
    """Normalize manual pair identities for the legacy analysis service API."""

    return [
        (canonical_residue_id(reference, "reference"), canonical_residue_id(target, "target"))
        for reference, target in pairs
    ]


__all__ = [
    "apply_selection_diagnostics",
    "canonical_manual_pairs",
    "canonical_optional_residue_id",
    "canonical_residue_id",
    "canonical_selection",
    "canonical_site_definition",
    "canonical_snapshot",
    "inputs_are_usable",
    "selection_diagnostics",
    "single_chain",
]
