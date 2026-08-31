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
from structlens.application.site_service import calculate_site_metrics, define_site
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
    """Index ligand components by component ID and residue name."""

    ligand_atoms: dict[str, Sequence[AtomRecord]] = {
        component.component_id: component.atoms
        for component in reference.components
        if component.is_ligand
    }
    ligand_atoms.update(
        {
            component.residue_name: component.atoms
            for component in reference.components
            if component.is_ligand and component.residue_name is not None
        }
    )
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
        site_values = sites_for_residue(item.reference, typed_sites, definitions, reference.residue_records)
        available = ["sequence", "structure"]
        unavailable: list[str] = []
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
                    available=item.target is not None and item.ca_displacement_angstrom is not None,
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
) -> tuple[SiteMetrics, ...]:
    relevant = {
        definition.site_id
        for definition in definitions
        if residue in {record.residue_id for record in define_site(definition, reference_records)}
    }
    return tuple(metric for metric in metrics if metric.site_id in relevant)


def site_metrics(
    definitions: Sequence[SiteDefinition],
    reference: ParsedStructure,
    target: ParsedStructure,
    result: AnalysisResult,
) -> tuple[SiteMetrics, ...]:
    """Calculate resolved site metrics using the fitted target-to-reference frame."""

    reference_chain = single_chain(reference)
    target_chain = single_chain(target)
    correspondence = {
        item.reference: item.target
        for item in result.correspondences
        if item.reference is not None and item.target is not None
    }
    matrix = homogeneous_transform(result.transform)
    ligand_atoms = reference_ligand_atoms(reference)
    return tuple(
        calculate_site_metrics(
            definition,
            reference_chain.residue_records,
            target_chain.residue_records,
            correspondence,
            target_structure_id=target_chain.structure_id,
            target_transform=matrix,
            ligand_atoms=ligand_atoms,
        )
        for definition in definitions
    )


__all__ = [
    "evidence_cards",
    "interactions_for_position",
    "reference_ligand_atoms",
    "site_metrics",
    "sites_for_residue",
]
