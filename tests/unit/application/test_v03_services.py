from __future__ import annotations

from structlens.application.difference_map_service import build_displacement_vectors, calculate_distance_difference
from structlens.application.interaction_service import InteractionAnalysisService
from structlens.application.msa_service import align_sequences, parse_alignment
from structlens.application.site_service import calculate_site_metrics
from structlens.core.interactions import InteractionType
from structlens.core.models import AtomRecord, ResidueId, ResidueNumbering, ResidueRecord
from structlens.core.msa import AnalysisSequence, MSASettings, SequenceResidueRef
from structlens.core.sites import SiteDefinition, SiteDefinitionMode


def _rid(structure: str, number: int, name: str = "ALA") -> ResidueId:
    return ResidueId(structure, "1", "A", str(number), None, name)


def _seq(identifier: str, value: str) -> AnalysisSequence:
    return AnalysisSequence(identifier, "A", value, tuple(SequenceResidueRef(i, char, _rid(identifier, i + 1)) for i, char in enumerate(value)), "structure")


def _residue(structure: str, number: int, name: str, xyz: tuple[float, float, float]) -> ResidueRecord:
    rid = _rid(structure, number, name)
    return ResidueRecord(rid, ResidueNumbering(str(number), str(number), None), name, name[0], (AtomRecord("CA", "C", xyz),))


def test_msa_fallback_preserves_reference_gap_and_conservation() -> None:
    result = align_sequences((_seq("ref", "ABC"), _seq("target", "ABCGH")), MSASettings())
    assert result.aligned_rows[0][1] == "ABC--"
    assert result.columns[3].reference_residue is None
    assert result.columns[3].reference_label == "A:3+1"
    assert result.columns[0].conservation_score == 1.0


def test_muscle_parser_rejects_duplicate_missing_and_unexpected_ids() -> None:
    text = ">SLSEQ000001\nAB\n>SLSEQ000002\nAB\n"
    assert parse_alignment(text, (_seq("ref", "AB"), _seq("target", "AB")))[0][1] == "AB"
    try:
        parse_alignment(">SLSEQ000001\nAB\n>SLSEQ000001\nAB\n", (_seq("ref", "AB"), _seq("target", "AB")))
    except ValueError as error:
        assert "duplicate" in str(error)


def test_distance_difference_and_vectors_are_reference_labelled() -> None:
    matrix = calculate_distance_difference(("A:1", "A:2"), {"A:1": (0, 0, 0), "A:2": (2, 0, 0)}, {"A:1": (0, 0, 0), "A:2": (3, 0, 0)})
    assert matrix.delta_angstrom[0, 1] == 1.0
    vectors = build_displacement_vectors(("A:1",), {"A:1": _rid("ref", 1)}, {"A:1": _rid("tar", 1)}, {"A:1": (0, 0, 0)}, {"A:1": (1, 0, 0)})
    assert vectors[0].magnitude_angstrom == 1.0


def test_interaction_service_uses_explicit_heavy_atom_evidence() -> None:
    residues = (_residue("ref", 1, "LYS", (0, 0, 0)), _residue("ref", 2, "ASP", (0, 0, 3)))
    records = InteractionAnalysisService().detect(residues)
    assert records and records[0].interaction_type in {InteractionType.HBOND_GEOMETRIC, InteractionType.SALT_BRIDGE}
    assert records[0].evidence_mode == "heavy_atom_geometry"


def test_site_metrics_report_coverage_and_unavailable_sasa() -> None:
    reference = (_residue("ref", 1, "SER", (0, 0, 0)), _residue("ref", 2, "LYS", (2, 0, 0)))
    target = (_residue("tar", 1, "SER", (0, 0, 0)),)
    definition = SiteDefinition("active", "Active", SiteDefinitionMode.KEY_RESIDUES, tuple(item.residue_id for item in reference), None, None, None)
    metrics = calculate_site_metrics(definition, reference, target, {reference[0].residue_id: target[0].residue_id}, target_structure_id="tar")
    assert metrics.mapped_residue_count == 1
    assert metrics.coverage_fraction == 0.5
    assert metrics.sasa_angstrom2 is None


def test_ligand_radius_site_uses_explicit_ligand_atoms() -> None:
    reference = (_residue("ref", 1, "SER", (0, 0, 0)), _residue("ref", 2, "LYS", (10, 0, 0)))
    target = (_residue("tar", 1, "SER", (0, 0, 0)),)
    definition = SiteDefinition("ligand", "Ligand site", SiteDefinitionMode.LIGAND_RADIUS, (), None, "LIG1", 2.0)
    ligand = AtomRecord("C1", "C", (0.5, 0, 0))
    metrics = calculate_site_metrics(
        definition,
        reference,
        target,
        {reference[0].residue_id: target[0].residue_id},
        target_structure_id="tar",
        ligand_atoms={"LIG1": (ligand,)},
    )
    assert metrics.mapped_residue_count == 1
    assert metrics.coverage_fraction == 1.0


def test_site_composition_uses_target_residue_chemistry() -> None:
    reference = (_residue("ref", 1, "SER", (0, 0, 0)),)
    target = (_residue("tar", 1, "VAL", (0, 0, 0)),)
    definition = SiteDefinition("active", "Active", SiteDefinitionMode.KEY_RESIDUES, (reference[0].residue_id,), None, None, None)
    metrics = calculate_site_metrics(
        definition,
        reference,
        target,
        {reference[0].residue_id: target[0].residue_id},
        target_structure_id="tar",
    )
    assert metrics.polar_residue_fraction == 0.0
    assert metrics.charged_residue_fraction == 0.0
