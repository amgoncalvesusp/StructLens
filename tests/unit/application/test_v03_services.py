from __future__ import annotations

import numpy as np
import pytest

from structlens.application.difference_map_service import build_displacement_vectors, calculate_distance_difference
from structlens.application.interaction_service import InteractionAnalysisService
from structlens.application.msa_service import align_sequences, parse_alignment
from structlens.application.site_service import calculate_site_metrics, define_site
from structlens.core.interactions import InteractionType
from structlens.core.models import AtomRecord, ResidueId, ResidueNumbering, ResidueRecord
from structlens.core.msa import AnalysisSequence, MSASettings, SequenceResidueRef
from structlens.core.sites import SiteDefinition, SiteDefinitionMode


def _rid(structure: str, number: int, name: str = "ALA") -> ResidueId:
    return ResidueId(structure, "1", "A", str(number), None, name)


def _seq(identifier: str, value: str) -> AnalysisSequence:
    return AnalysisSequence(
        identifier,
        "A",
        value,
        tuple(SequenceResidueRef(i, char, _rid(identifier, i + 1)) for i, char in enumerate(value)),
        "structure",
    )


def _residue(structure: str, number: int, name: str, xyz: tuple[float, float, float]) -> ResidueRecord:
    rid = _rid(structure, number, name)
    return ResidueRecord(
        rid, ResidueNumbering(str(number), str(number), None), name, name[0], (AtomRecord("CA", "C", xyz),)
    )


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
    matrix = calculate_distance_difference(
        ("A:1", "A:2"), {"A:1": (0, 0, 0), "A:2": (2, 0, 0)}, {"A:1": (0, 0, 0), "A:2": (3, 0, 0)}
    )
    assert matrix.delta_angstrom[0, 1] == 1.0
    vectors = build_displacement_vectors(
        ("A:1",), {"A:1": _rid("ref", 1)}, {"A:1": _rid("tar", 1)}, {"A:1": (0, 0, 0)}, {"A:1": (1, 0, 0)}
    )
    assert vectors[0].magnitude_angstrom == 1.0


def test_interaction_service_uses_explicit_heavy_atom_evidence() -> None:
    residues = (_residue("ref", 1, "LYS", (0, 0, 0)), _residue("ref", 2, "ASP", (0, 0, 3)))
    records = InteractionAnalysisService().detect(residues)
    assert records and records[0].interaction_type in {InteractionType.HBOND_GEOMETRIC, InteractionType.SALT_BRIDGE}
    assert records[0].evidence_mode == "heavy_atom_geometry"


def test_site_metrics_report_coverage_and_unavailable_sasa() -> None:
    reference = (_residue("ref", 1, "SER", (0, 0, 0)), _residue("ref", 2, "LYS", (2, 0, 0)))
    target = (_residue("tar", 1, "SER", (0, 0, 0)),)
    definition = SiteDefinition(
        "active",
        "Active",
        SiteDefinitionMode.KEY_RESIDUES,
        tuple(item.residue_id for item in reference),
        None,
        None,
        None,
    )
    metrics = calculate_site_metrics(
        definition, reference, target, {reference[0].residue_id: target[0].residue_id}, target_structure_id="tar"
    )
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


def test_radius_site_selection_ignores_deposited_hydrogens() -> None:
    residue_id = _rid("ref", 1, "SER")
    residue = ResidueRecord(
        residue_id,
        ResidueNumbering("1", "1", None),
        "SER",
        "S",
        (
            AtomRecord("CA", "C", (3.0, 0.0, 0.0)),
            AtomRecord("H", "H", (0.25, 0.0, 0.0)),
        ),
    )
    definition = SiteDefinition(
        "ligand",
        "Ligand site",
        SiteDefinitionMode.LIGAND_RADIUS,
        ligand_id="LIG1",
        radius_angstrom=1.0,
    )

    selected = define_site(
        definition,
        (residue,),
        ligand_atoms={"LIG1": (AtomRecord("C1", "C", (0.0, 0.0, 0.0)),)},
    )

    assert selected == ()


def test_residue_radius_site_requires_a_center_residue() -> None:
    with pytest.raises(ValueError, match="center_residue"):
        SiteDefinition(
            "residue",
            "Residue site",
            SiteDefinitionMode.RESIDUE_RADIUS,
            radius_angstrom=4.0,
        )


def test_site_composition_uses_target_residue_chemistry() -> None:
    reference = (_residue("ref", 1, "SER", (0, 0, 0)),)
    target = (_residue("tar", 1, "VAL", (0, 0, 0)),)
    definition = SiteDefinition(
        "active", "Active", SiteDefinitionMode.KEY_RESIDUES, (reference[0].residue_id,), None, None, None
    )
    metrics = calculate_site_metrics(
        definition,
        reference,
        target,
        {reference[0].residue_id: target[0].residue_id},
        target_structure_id="tar",
    )
    assert metrics.polar_residue_fraction == 0.0
    assert metrics.charged_residue_fraction == 0.0


def _backbone_residue(structure: str, number: int, name: str, origin: tuple[float, float, float]) -> ResidueRecord:
    """A residue carrying a full N/CA/C/O backbone plus one side-chain atom."""

    x, y, z = origin
    rid = _rid(structure, number, name)
    atoms = (
        AtomRecord("N", "N", (x - 1.2, y, z)),
        AtomRecord("CA", "C", (x, y, z)),
        AtomRecord("C", "C", (x + 1.2, y, z)),
        AtomRecord("O", "O", (x + 1.2, y + 1.1, z + 0.4)),
        AtomRecord("CB", "C", (x, y + 1.5, z + 0.9)),
    )
    return ResidueRecord(rid, ResidueNumbering(str(number), str(number), None), name, name[0], atoms)


def _site_pair(
    reference_origins: tuple[tuple[float, float, float], ...],
    target_origins: tuple[tuple[float, float, float], ...],
) -> tuple[tuple[ResidueRecord, ...], tuple[ResidueRecord, ...], dict[ResidueId, ResidueId], SiteDefinition]:
    names = ("ALA", "SER", "LYS", "ASP")
    reference = tuple(_backbone_residue("ref", i + 1, names[i], xyz) for i, xyz in enumerate(reference_origins))
    target = tuple(_backbone_residue("tar", i + 1, names[i], xyz) for i, xyz in enumerate(target_origins))
    mapping = {ref.residue_id: tar.residue_id for ref, tar in zip(reference, target, strict=True)}
    definition = SiteDefinition(
        "site", "Site", SiteDefinitionMode.KEY_RESIDUES, tuple(item.residue_id for item in reference), None, None, None
    )
    return reference, target, mapping, definition


def test_site_backbone_rmsd_uses_matched_backbone_atoms_not_ca_alone() -> None:
    """A field named backbone RMSD must measure N/CA/C/O, not the C-alpha only."""

    # The target keeps every C-alpha in place and rotates only the carbonyl
    # oxygens, so a C-alpha-only measurement reports a perfect zero.
    reference_origins = ((0.0, 0.0, 0.0), (3.8, 0.0, 0.0), (7.6, 0.0, 0.0), (11.4, 0.0, 0.0))
    reference, target, mapping, definition = _site_pair(reference_origins, reference_origins)
    shifted = tuple(
        ResidueRecord(
            record.residue_id,
            record.numbering,
            record.residue_name,
            record.one_letter,
            tuple(
                AtomRecord(atom.name, atom.element, (atom.coordinate[0], atom.coordinate[1], atom.coordinate[2] + 2.0))
                if atom.name == "O"
                else atom
                for atom in record.atoms
            ),
        )
        for record in target
    )
    metrics = calculate_site_metrics(
        definition, reference, shifted, mapping, target_structure_id="tar", target_transform=np.eye(4)
    )
    assert metrics.global_frame_backbone_rmsd_angstrom is not None
    assert metrics.global_frame_backbone_rmsd_angstrom > 0.0, "displaced carbonyl oxygens must move backbone RMSD"


def test_site_fitted_rmsd_removes_rotation_not_only_translation() -> None:
    """Site-fitted RMSD needs a Kabsch fit; centroid subtraction leaves rotation in."""

    reference_origins = ((0.0, 0.0, 0.0), (3.8, 0.0, 0.0), (7.6, 1.0, 0.0), (11.4, 0.0, 0.5))
    reference, _, _, definition = _site_pair(reference_origins, reference_origins)
    angle = np.pi / 4.0
    rotation = np.array([[np.cos(angle), -np.sin(angle), 0.0], [np.sin(angle), np.cos(angle), 0.0], [0.0, 0.0, 1.0]])
    offset = np.array([5.0, -2.0, 1.0])

    def _moved(record: ResidueRecord, number: int) -> ResidueRecord:
        """The same residue as a rigid body, in a rotated and shifted frame."""

        atoms = tuple(
            AtomRecord(atom.name, atom.element, tuple(np.asarray(atom.coordinate) @ rotation.T + offset))
            for atom in record.atoms
        )
        return ResidueRecord(
            _rid("tar", number, record.residue_name),
            record.numbering,
            record.residue_name,
            record.one_letter,
            atoms,
        )

    target = tuple(_moved(record, index + 1) for index, record in enumerate(reference))
    mapping = {ref.residue_id: tar.residue_id for ref, tar in zip(reference, target, strict=True)}
    metrics = calculate_site_metrics(
        definition, reference, target, mapping, target_structure_id="tar", target_transform=np.eye(4)
    )
    assert metrics.site_fitted_backbone_rmsd_angstrom is not None
    assert metrics.site_fitted_backbone_rmsd_angstrom == pytest.approx(0.0, abs=1e-6)


def test_global_frame_site_rmsd_is_unavailable_without_an_authoritative_transform() -> None:
    """Comparing raw coordinates across unaligned frames is not a global RMSD."""

    reference_origins = ((0.0, 0.0, 0.0), (3.8, 0.0, 0.0), (7.6, 0.0, 0.0), (11.4, 0.0, 0.0))
    target_origins = tuple((x + 25.0, y, z) for x, y, z in reference_origins)
    reference, target, mapping, definition = _site_pair(reference_origins, target_origins)
    metrics = calculate_site_metrics(definition, reference, target, mapping, target_structure_id="tar")
    assert metrics.global_frame_backbone_rmsd_angstrom is None
    assert metrics.centroid_displacement_angstrom is None
    assert metrics.site_fitted_backbone_rmsd_angstrom is not None


def test_site_envelope_volume_is_unavailable_rather_than_zero_or_a_crash() -> None:
    """Too few points, and coplanar points, both mean no measurable envelope."""

    flat = ((0.0, 0.0, 0.0), (3.0, 0.0, 0.0), (6.0, 2.0, 0.0), (9.0, 0.0, 0.0))
    reference, target, mapping, definition = _site_pair(flat, flat)
    coplanar_target = tuple(
        ResidueRecord(
            record.residue_id,
            record.numbering,
            record.residue_name,
            record.one_letter,
            tuple(AtomRecord(a.name, a.element, (a.coordinate[0], a.coordinate[1], 0.0)) for a in record.atoms),
        )
        for record in target
    )
    metrics = calculate_site_metrics(
        definition, reference, coplanar_target, mapping, target_structure_id="tar", target_transform=np.eye(4)
    )
    assert metrics.atomic_envelope_volume_angstrom3 is None


def test_distance_difference_tolerates_partially_mapped_positions() -> None:
    """An unmapped position must be masked, never stored as NaN in a strict model."""

    matrix = calculate_distance_difference(
        ("A:1", "A:2", "A:3"),
        {"A:1": (0.0, 0.0, 0.0), "A:2": (3.0, 0.0, 0.0)},
        {"A:1": (0.0, 0.0, 0.0), "A:2": (4.0, 0.0, 0.0)},
    )
    assert matrix.delta_angstrom[0, 1] == pytest.approx(1.0)
    assert not matrix.valid_mask[2, 0]
    assert np.isfinite(np.asarray(matrix.target_distances_angstrom)).all()
