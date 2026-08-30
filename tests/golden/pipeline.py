"""Deterministic v0.3 scientific pipeline used by the golden regression test.

The module is import-only: it turns the two checked-in PDB fixtures into one
JSON-serializable snapshot covering every scientific lane named in the v0.3
release gates (MSA, conservation, interactions, sites, distance maps, vectors,
and Evidence Cards). Keeping it separate from the assertions lets the same code
regenerate the golden file and verify it.

The fixtures are synthetic: bond lengths and inter-atom distances are chosen to
exercise every branch deterministically, not to be chemically realistic. The
detected hydrogen bonds in particular sit well below a physiological 2.8-3.2 A,
so the golden pins the arithmetic rather than validating the chemistry. Replace
the fixtures with a curated real structure pair if the interaction thresholds
themselves ever need regression cover.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from structlens.application.difference_map_service import build_displacement_vectors, calculate_distance_difference
from structlens.application.interaction_service import InteractionAnalysisService
from structlens.application.msa_service import align_sequences
from structlens.application.site_service import calculate_site_metrics
from structlens.core.evidence import InteractionEvidence, SequenceEvidence, SiteEvidence, StructureEvidence
from structlens.core.evidence.builder import build_evidence_card
from structlens.core.evidence.completeness import quality_for_sections
from structlens.core.interactions.comparison import compare_interactions
from structlens.core.models import AtomRecord, ProteinChain, ResidueId, ResidueRecord
from structlens.core.msa import AnalysisSequence, MSASettings, MultipleSequenceAlignment, SequenceResidueRef
from structlens.core.parsing.pdb import load_pdb
from structlens.core.sites import SiteDefinition, SiteDefinitionMode

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "golden"
LIGAND_NAME = "LIG"


def _round(value: float | None, digits: int = 6) -> float | None:
    return None if value is None else round(float(value), digits)


def _protein_residues(chain: ProteinChain) -> tuple[ResidueRecord, ...]:
    return tuple(record for record in chain.residue_records if record.is_standard)


def _ligand_atoms(chain: ProteinChain) -> dict[str, tuple[AtomRecord, ...]]:
    return {record.residue_name: record.atoms for record in chain.residue_records if record.residue_name == LIGAND_NAME}


def _analysis_sequence(chain: ProteinChain) -> AnalysisSequence:
    residues = _protein_residues(chain)
    letters = "".join(record.one_letter or "X" for record in residues)
    refs = tuple(SequenceResidueRef(index, letters[index], record.residue_id) for index, record in enumerate(residues))
    return AnalysisSequence(chain.structure_id, chain.chain_id, letters, refs, "structure")


def _ca(record: ResidueRecord) -> tuple[float, float, float] | None:
    for atom in record.atoms:
        if atom.name == "CA":
            x, y, z = (float(value) for value in atom.coordinate)
            return (x, y, z)
    return None


def _correspondence(alignment: MultipleSequenceAlignment) -> dict[str, dict[str, ResidueId]]:
    """Map each fully aligned column label to the residues it links per structure."""

    mapped: dict[str, dict[str, ResidueId]] = {}
    for column in alignment.columns:
        entry: dict[str, ResidueId] = {}
        for cell in column.cells:
            if cell.residue is not None and cell.residue.residue_id is not None:
                entry[cell.structure_id] = cell.residue.residue_id
        if len(entry) == len(alignment.sequences):
            mapped[column.reference_label] = entry
    return mapped


def _msa_section(alignment: MultipleSequenceAlignment) -> dict[str, Any]:
    return {
        "algorithm": alignment.algorithm,
        "reference_structure_id": alignment.reference_structure_id,
        "aligned_rows": [[identifier, row] for identifier, row in alignment.aligned_rows],
        "columns": [
            {
                "index": column.index,
                "reference_label": column.reference_label,
                "reference_auth_seq_id": (
                    None if column.reference_residue is None else column.reference_residue.auth_seq_id
                ),
                "characters": [cell.character for cell in column.cells],
                "non_gap_count": column.non_gap_count,
            }
            for column in alignment.columns
        ],
    }


def _conservation_section(alignment: MultipleSequenceAlignment) -> list[dict[str, Any]]:
    return [
        {
            "index": column.index,
            "reference_label": column.reference_label,
            "conservation_score": _round(column.conservation_score),
            "entropy_bits": _round(column.entropy_bits),
            "gap_fraction": _round(column.gap_fraction),
            "ambiguous_fraction": _round(column.ambiguous_fraction),
        }
        for column in alignment.columns
    ]


def _interaction_row(record: Any) -> dict[str, Any]:
    return {
        "type": record.interaction_type.value,
        "evidence_mode": record.evidence_mode,
        "distance_angstrom": _round(record.distance_angstrom, 3),
    }


def _interaction_section(
    reference: tuple[ResidueRecord, ...],
    target: tuple[ResidueRecord, ...],
    correspondence: dict[str, dict[str, ResidueId]],
    reference_id: str,
    target_id: str,
) -> dict[str, Any]:
    service = InteractionAnalysisService()
    reference_records = service.detect(reference)
    target_records = service.detect(target)
    position_by_residue: dict[ResidueId, str] = {}
    for label, entry in correspondence.items():
        position_by_residue[entry[reference_id]] = label
        position_by_residue[entry[target_id]] = label
    differences = compare_interactions(
        reference_records,
        target_records,
        lambda residue: position_by_residue.get(residue),
    )
    return {
        "reference": [_interaction_row(record) for record in reference_records],
        "target": [_interaction_row(record) for record in target_records],
        "differences": [
            {
                "position_a": difference.key.reference_position_a,
                "position_b": difference.key.reference_position_b,
                "type": difference.key.interaction_type.value,
                "change": difference.change.value,
            }
            for difference in differences
        ],
    }


def _site_metrics_row(metrics: Any) -> dict[str, Any]:
    return {
        "site_id": metrics.site_id,
        "mapped_residue_count": metrics.mapped_residue_count,
        "coverage_fraction": _round(metrics.coverage_fraction),
        "global_frame_backbone_rmsd_angstrom": _round(metrics.global_frame_backbone_rmsd_angstrom),
        "site_fitted_backbone_rmsd_angstrom": _round(metrics.site_fitted_backbone_rmsd_angstrom),
        "atomic_envelope_volume_angstrom3": _round(metrics.atomic_envelope_volume_angstrom3),
        "sasa_angstrom2": _round(metrics.sasa_angstrom2),
        "polar_residue_fraction": _round(metrics.polar_residue_fraction),
        "charged_residue_fraction": _round(metrics.charged_residue_fraction),
    }


def _site_section(
    reference: tuple[ResidueRecord, ...],
    target: tuple[ResidueRecord, ...],
    residue_mapping: dict[ResidueId, ResidueId],
    ligand_atoms: dict[str, tuple[AtomRecord, ...]],
    target_id: str,
) -> list[dict[str, Any]]:
    definitions = (
        SiteDefinition(
            "key",
            "Key residues",
            SiteDefinitionMode.KEY_RESIDUES,
            tuple(record.residue_id for record in reference[:3]),
            None,
            None,
            None,
        ),
        SiteDefinition("ligand", "Ligand radius", SiteDefinitionMode.LIGAND_RADIUS, (), None, LIGAND_NAME, 6.0),
        SiteDefinition(
            "radius",
            "Residue radius",
            SiteDefinitionMode.RESIDUE_RADIUS,
            (),
            reference[1].residue_id,
            None,
            6.0,
        ),
    )
    return [
        _site_metrics_row(
            calculate_site_metrics(
                definition,
                reference,
                target,
                residue_mapping,
                target_structure_id=target_id,
                ligand_atoms=ligand_atoms,
            )
        )
        for definition in definitions
    ]


def build_snapshot() -> dict[str, Any]:
    """Run every v0.3 scientific lane over the checked-in fixtures."""

    reference_chain = load_pdb(FIXTURES / "reference.pdb").chains[0]
    target_chain = load_pdb(FIXTURES / "target.pdb").chains[0]
    reference_id = reference_chain.structure_id
    target_id = target_chain.structure_id

    reference_residues = _protein_residues(reference_chain)
    target_residues = _protein_residues(target_chain)

    alignment = align_sequences(
        (_analysis_sequence(reference_chain), _analysis_sequence(target_chain)),
        MSASettings(),
    )
    correspondence = _correspondence(alignment)
    residue_mapping = {entry[reference_id]: entry[target_id] for entry in correspondence.values()}

    reference_by_id = {record.residue_id: record for record in reference_residues}
    target_by_id = {record.residue_id: record for record in target_residues}
    reference_coordinates: dict[str, tuple[float, float, float]] = {}
    target_coordinates: dict[str, tuple[float, float, float]] = {}
    for label in sorted(correspondence):
        entry = correspondence[label]
        reference_ca = _ca(reference_by_id[entry[reference_id]])
        target_ca = _ca(target_by_id[entry[target_id]])
        if reference_ca is None or target_ca is None:
            continue
        reference_coordinates[label] = reference_ca
        target_coordinates[label] = target_ca
    labels = tuple(reference_coordinates)

    matrix = calculate_distance_difference(labels, reference_coordinates, target_coordinates)
    vectors = build_displacement_vectors(
        labels,
        {label: correspondence[label][reference_id] for label in labels},
        {label: correspondence[label][target_id] for label in labels},
        reference_coordinates,
        target_coordinates,
    )

    first_label = labels[0]
    first_column = next(column for column in alignment.columns if column.reference_label == first_label)
    card = build_evidence_card(
        correspondence[first_label][reference_id],
        target_id=target_id,
        sequence=SequenceEvidence(
            reference_one_letter=reference_by_id[correspondence[first_label][reference_id]].one_letter,
            target_one_letter=target_by_id[correspondence[first_label][target_id]].one_letter,
            alignment_index=first_column.index,
            conservation_fraction=first_column.conservation_score,
            entropy_bits=first_column.entropy_bits,
            gap_fraction=first_column.gap_fraction,
            ambiguous_fraction=first_column.ambiguous_fraction,
        ),
        structure=StructureEvidence(ca_displacement_angstrom=vectors[0].magnitude_angstrom, available=True),
        interactions=InteractionEvidence(),
        site=SiteEvidence(),
        quality=quality_for_sections(("sequence", "structure"), ("interactions", "site")),
        provenance=("golden-fixture",),
    )
    sequence_evidence = card.sequence
    structure_evidence = card.structure
    quality = card.quality
    if sequence_evidence is None or structure_evidence is None or quality is None:
        raise AssertionError("golden evidence card must carry sequence, structure, and quality sections")

    return {
        "msa": _msa_section(alignment),
        "conservation": _conservation_section(alignment),
        "interactions": _interaction_section(
            reference_residues, target_residues, correspondence, reference_id, target_id
        ),
        "sites": _site_section(
            reference_residues,
            target_residues,
            residue_mapping,
            _ligand_atoms(reference_chain),
            target_id,
        ),
        "distance_map": {
            "labels": list(matrix.reference_positions),
            "valid_mask": [[bool(value) for value in row] for row in matrix.valid_mask.tolist()],
            "delta_angstrom": [[_round(value, 6) for value in row] for row in matrix.delta_angstrom.tolist()],
        },
        "vectors": [
            {
                "reference_position": vector.reference_position,
                "reference_auth_seq_id": vector.reference_residue.auth_seq_id,
                "target_auth_seq_id": vector.target_residue.auth_seq_id,
                "vector_xyz": [_round(value) for value in vector.vector_xyz],
                "magnitude_angstrom": _round(vector.magnitude_angstrom),
            }
            for vector in vectors
        ],
        "evidence_card": {
            "reference_auth_seq_id": card.reference_residue.auth_seq_id,
            "target_id": card.target_id,
            "reference_one_letter": sequence_evidence.reference_one_letter,
            "target_one_letter": sequence_evidence.target_one_letter,
            "conservation_fraction": _round(sequence_evidence.conservation_fraction),
            "entropy_bits": _round(sequence_evidence.entropy_bits),
            "ca_displacement_angstrom": _round(structure_evidence.ca_displacement_angstrom),
            "quality": {
                "overall_status": quality.overall_status,
                "available_sections": list(quality.available_sections),
                "unavailable_sections": list(quality.unavailable_sections),
                "source_count": quality.source_count,
            },
            "provenance": list(card.provenance),
        },
    }


__all__ = ["FIXTURES", "build_snapshot"]
