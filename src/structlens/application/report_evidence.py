"""Evidence-card and interaction/site evidence assembly for reports."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

from structlens.application.report_geometry import (
    analysis_sequence,
    homogeneous_transform,
    position_lookup,
    reference_position,
)
from structlens.application.report_input import single_chain
from structlens.application.site_service import define_site, measure_resolved_site
from structlens.core.evidence import (
    Availability,
    DiagnosticSeverity,
    EvidenceCard,
    InteractionEvidence,
    SequenceEvidence,
    SiteEvidence,
    StructureEvidence,
    build_evidence_card,
    quality_for_sections,
)
from structlens.core.models import AnalysisResult, AtomRecord, ProteinChain, ResidueId, ResidueRecord
from structlens.core.msa import MultipleSequenceAlignment, SequenceResidueRef
from structlens.core.parsing import ParsedStructure
from structlens.core.reports import InputQualityBundle
from structlens.core.sites import SiteDefinition, SiteMetrics


def reference_ligand_atoms(reference: ParsedStructure) -> dict[str, Sequence[AtomRecord]]:
    """Index ligands by exact ID and only unambiguous residue-name aliases."""

    ligands = tuple(component for component in reference.components if component.is_ligand)
    ligand_atoms: dict[str, Sequence[AtomRecord]] = {
        component.component_id: component.atoms for component in ligands
    }
    by_name: dict[str, list[Sequence[AtomRecord]]] = {}
    for component in ligands:
        if component.residue_name is not None:
            by_name.setdefault(component.residue_name, []).append(component.atoms)
    for residue_name, candidates in by_name.items():
        if len(candidates) == 1 and residue_name not in ligand_atoms:
            ligand_atoms[residue_name] = candidates[0]
    return ligand_atoms


def evidence_cards(
    result: AnalysisResult,
    reference: ProteinChain,
    msa: object | None,
    interactions: object | None,
    sites: object | None,
    definitions: Sequence[SiteDefinition],
    quality: InputQualityBundle,
    interaction_state: Availability,
    site_state: Availability,
    *,
    resolved: Mapping[SiteDefinition, Sequence[ResidueRecord]] | None = None,
    ligand_atoms: Mapping[str, Sequence[AtomRecord]] | None = None,
) -> tuple[EvidenceCard, ...]:
    """Assemble one card per mapped reference residue with explicit availability."""

    typed_msa = msa if isinstance(msa, MultipleSequenceAlignment) else None
    typed_interactions = interactions if isinstance(interactions, InteractionEvidence) else InteractionEvidence()
    typed_sites = tuple(sites) if isinstance(sites, tuple) and all(isinstance(item, SiteMetrics) for item in sites) else ()
    sequence_refs = {item.residue_id: item for item in analysis_sequence(reference).residues}
    columns = (
        {
            column.reference_residue: column
            for column in typed_msa.columns
            if column.reference_residue is not None
        }
        if typed_msa is not None
        else {}
    )
    position_for = position_lookup(result.correspondences)
    warnings = tuple(
        diagnostic.code
        for diagnostic in quality.reference.diagnostics + quality.target.diagnostics
        if diagnostic.severity is DiagnosticSeverity.WARNING
    )
    metrics_by_site = {metric.site_id: metric for metric in typed_sites}
    if resolved is None:
        resolved = {
            definition: define_site(definition, reference.residue_records, ligand_atoms=ligand_atoms)
            for definition in definitions
        }
    sites_by_member: dict[ResidueId, list[SiteMetrics]] = {}
    for definition, members in resolved.items():
        metric = metrics_by_site.get(definition.site_id)
        if metric is not None:
            for member in members:
                sites_by_member.setdefault(member.residue_id, []).append(metric)
    cards: list[EvidenceCard] = []
    for item in result.correspondences:
        if item.reference is None:
            continue
        residue_ref = sequence_refs.get(item.reference)
        if residue_ref is None:
            residue_ref = SequenceResidueRef(item.alignment_index, item.reference_one_letter or "X", item.reference)
        column = columns.get(item.reference)
        sequence = SequenceEvidence(
            reference_one_letter=item.reference_one_letter,
            target_one_letter=item.target_one_letter,
            alignment_index=item.alignment_index,
            sequence_identity=result.sequence_identity,
            conservation_fraction=column.conservation_score if column is not None else None,
            entropy_bits=column.entropy_bits if column is not None else None,
            gap_fraction=column.gap_fraction if column is not None else None,
            ambiguous_fraction=column.ambiguous_fraction if column is not None else None,
            source_refs=tuple(cell.residue for cell in column.cells if cell.residue is not None)
            if column is not None
            else (residue_ref,),
        )
        position = reference_position(item.reference)
        interaction_evidence = interactions_for_position(typed_interactions, position, position_for)
        site_values = tuple(sites_by_member.get(item.reference, ()))
        structure_available = item.target is not None and item.ca_displacement_angstrom is not None
        available = ["sequence"]
        unavailable: list[str] = []
        if structure_available:
            available.append("structure")
        else:
            unavailable.append("structure")
        if interaction_state is Availability.AVAILABLE:
            available.append("interactions")
        else:
            unavailable.append("interactions")
        if site_state is Availability.AVAILABLE:
            available.append("sites")
        else:
            unavailable.append("sites")
        cards.append(
            build_evidence_card(
                residue_ref,
                target_id=result.target_id,
                sequence=sequence,
                structure=StructureEvidence(
                    item.ca_displacement_angstrom,
                    item.backbone_rmsd_angstrom,
                    item.sidechain_rmsd_angstrom,
                    item.all_heavy_atom_rmsd_angstrom,
                    available=structure_available,
                ),
                interactions=interaction_evidence,
                site=SiteEvidence(site_values),
                quality=quality_for_sections(available, unavailable, warnings=warnings),
                provenance=("structlens.analysis_report", result.alignment_decision),
            )
        )
    return tuple(cards)


def interactions_for_position(
    evidence: InteractionEvidence,
    position: str,
    position_for: Mapping[ResidueId, str],
) -> InteractionEvidence:
    differences = tuple(
        item
        for item in evidence.differences
        if position in {item.key.reference_position_a, item.key.reference_position_b}
    )
    reference = tuple(
        item
        for item in evidence.reference_interactions
        if position
        in {
            position_for.get(item.residue_a),
            position_for.get(item.residue_b) if item.residue_b else None,
        }
    )
    target = tuple(
        item
        for item in evidence.target_interactions
        if position
        in {
            position_for.get(item.residue_a),
            position_for.get(item.residue_b) if item.residue_b else None,
        }
    )
    return InteractionEvidence(differences, reference, target)


def sites_for_residue(
    residue: ResidueId,
    metrics: Sequence[SiteMetrics],
    definitions: Sequence[SiteDefinition],
    reference_records: Sequence[ResidueRecord],
    *,
    ligand_atoms: Mapping[str, Sequence[AtomRecord]] | None = None,
) -> tuple[SiteMetrics, ...]:
    relevant = {
        definition.site_id
        for definition in definitions
        if any(
            record.residue_id == residue
            for record in define_site(definition, reference_records, ligand_atoms=ligand_atoms)
        )
    }
    return tuple(metric for metric in metrics if metric.site_id in relevant)


def site_metrics(
    definitions: Sequence[SiteDefinition],
    reference: ParsedStructure,
    target: ParsedStructure,
    result: AnalysisResult,
    *,
    resolved: Mapping[SiteDefinition, Sequence[ResidueRecord]] | None = None,
) -> tuple[SiteMetrics, ...]:
    """Calculate resolved site metrics using the fitted target-to-reference frame."""

    target_chain = single_chain(target)
    correspondence = {
        item.reference: item.target
        for item in result.correspondences
        if item.reference is not None and item.target is not None
    }
    matrix = homogeneous_transform(result.transform)
    if resolved is None:
        records = single_chain(reference).residue_records
        atoms = reference_ligand_atoms(reference)
        resolved = {definition: define_site(definition, records, ligand_atoms=atoms) for definition in definitions}
    return tuple(
        measure_resolved_site(
            definition,
            members,
            target_chain.residue_records,
            correspondence,
            target_structure_id=target_chain.structure_id,
            target_transform=matrix,
        )
        for definition, members in resolved.items()
    )


__all__ = [
    "evidence_cards",
    "interactions_for_position",
    "reference_ligand_atoms",
    "site_metrics",
    "sites_for_residue",
]
