from __future__ import annotations

import pytest

from structlens.core.models import AtomRecord, ComponentKind, ResidueId, StructureComponent
from structlens.core.pockets import AlphaSphere, PocketCandidate
from structlens.core.pockets.ligands import (
    POCKET_LIGAND_RULES_VERSION,
    PocketLigandSupport,
    eligible_pocket_ligands,
    measure_ligand_support,
)


def _residue(index: int) -> ResidueId:
    return ResidueId("protein", "1", "A", str(index), None, "ALA")


def _sphere(index: int, center: tuple[float, float, float], radius: float) -> AlphaSphere:
    atom_ids = tuple(f"sphere-{index}-atom-{offset}" for offset in range(4))
    return AlphaSphere(
        center_xyz=center,
        radius_angstrom=radius,
        touching_atom_ids=atom_ids,
        lining_residues=(_residue(index), _residue(index + 10)),
        source_simplex_atom_ids=atom_ids,
    )


def _ligand(
    component_id: str,
    residue_name: str,
    *atoms: tuple[float, float, float],
    kind: ComponentKind = ComponentKind.LIGAND,
) -> StructureComponent:
    records = tuple(
        AtomRecord(
            name=f"C{index}",
            element="C",
            coordinate=coordinate,
            source_atom_id=f"{component_id}-atom-{index}",
        )
        for index, coordinate in enumerate(atoms, start=1)
    )
    return StructureComponent(
        component_id=component_id,
        kind=kind,
        atoms=records,
        metadata={"selected_for_analysis": True},
        model_id="1",
        author_chain_id="A",
        residue_name=residue_name,
        auth_seq_id=component_id,
    )


def test_eligible_pocket_ligands_exclude_water_ions_buffers_and_other_components() -> None:
    atp = _ligand("ATP-1", "ATP", (0.0, 0.0, 0.0))
    sulfate = _ligand("SO4-1", "SO4", (1.0, 0.0, 0.0))
    magnesium = _ligand("MG-1", "MG", (2.0, 0.0, 0.0), kind=ComponentKind.ION)
    water = _ligand("HOH-1", "HOH", (3.0, 0.0, 0.0), kind=ComponentKind.WATER)
    unknown = _ligand("UNK-1", "UNK", (4.0, 0.0, 0.0), kind=ComponentKind.OTHER)

    eligible = eligible_pocket_ligands((water, magnesium, sulfate, unknown, atp))

    assert POCKET_LIGAND_RULES_VERSION == "structlens-pocket-known-ligands-1"
    assert tuple(component.component_id for component in eligible) == ("ATP-1",)


def test_measure_ligand_support_reports_distance_coverage_and_lining_overlap() -> None:
    candidate = PocketCandidate((_sphere(1, (0.0, 0.0, 0.0), 2.0),))
    ligand = _ligand(
        "ATP-1",
        "ATP",
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (3.0, 0.0, 0.0),
    )

    support = measure_ligand_support(
        candidate,
        ligand,
        ligand_contact_residues=(_residue(1),),
    )

    assert isinstance(support, PocketLigandSupport)
    assert support.component_id == "ATP-1"
    assert support.residue_name == "ATP"
    assert support.atom_count == 3
    assert support.covered_atom_count == 2
    assert support.atom_coverage_fraction == pytest.approx(2.0 / 3.0)
    assert support.ligand_center_distance_angstrom == pytest.approx(4.0 / 3.0)
    assert support.lining_residue_overlap_fraction == pytest.approx(0.5)


def test_ligand_support_rejects_inconsistent_fraction_and_count_state() -> None:
    with pytest.raises(ValueError, match="atom_coverage_fraction"):
        PocketLigandSupport(
            component_id="ATP-1",
            residue_name="ATP",
            ligand_center_distance_angstrom=1.0,
            atom_coverage_fraction=0.1,
            lining_residue_overlap_fraction=0.0,
            covered_atom_count=2,
            atom_count=2,
        )


@pytest.mark.parametrize(
    ("field", "value", "message"),
    (
        ("component_id", "   ", "component_id"),
        ("residue_name", "   ", "residue_name"),
        ("ligand_center_distance_angstrom", float("inf"), "finite"),
        ("ligand_center_distance_angstrom", -1.0, "non-negative"),
        ("atom_coverage_fraction", 1.1, "must not exceed 1"),
        ("covered_atom_count", True, "non-negative integer"),
        ("atom_count", 0, "positive"),
    ),
)
def test_ligand_support_rejects_invalid_field_values(field: str, value: object, message: str) -> None:
    payload = {
        "component_id": "ATP-1",
        "residue_name": "ATP",
        "ligand_center_distance_angstrom": 1.0,
        "atom_coverage_fraction": 0.5,
        "lining_residue_overlap_fraction": 0.5,
        "covered_atom_count": 1,
        "atom_count": 2,
    }
    payload[field] = value

    with pytest.raises(ValueError, match=message):
        PocketLigandSupport(**payload)


def test_ligand_support_rejects_covered_atoms_above_atom_count() -> None:
    with pytest.raises(ValueError, match="covered_atom_count"):
        PocketLigandSupport(
            component_id="ATP-1",
            residue_name="ATP",
            ligand_center_distance_angstrom=1.0,
            atom_coverage_fraction=1.0,
            lining_residue_overlap_fraction=0.0,
            covered_atom_count=3,
            atom_count=2,
        )


def test_eligible_pocket_ligands_require_selected_non_excluded_ligands_and_sort_by_component_id() -> None:
    atp = _ligand("ATP-2", " atp ", (0.0, 0.0, 0.0))
    adp = _ligand("ADP-1", "ADP", (1.0, 0.0, 0.0))
    glycerol = _ligand("GOL-1", "GOL", (2.0, 0.0, 0.0))
    not_selected = StructureComponent(
        component_id="ATP-3",
        kind=ComponentKind.LIGAND,
        atoms=atp.atoms,
        metadata={"selected_for_analysis": False},
        model_id="1",
        author_chain_id="A",
        residue_name="ATP",
        auth_seq_id="ATP-3",
    )

    eligible = eligible_pocket_ligands((atp, glycerol, not_selected, adp))

    assert tuple(component.component_id for component in eligible) == ("ADP-1", "ATP-2")


def test_measure_ligand_support_rejects_untyped_inputs() -> None:
    candidate = PocketCandidate((_sphere(1, (0.0, 0.0, 0.0), 2.0),))
    ligand = _ligand("ATP-1", "ATP", (0.0, 0.0, 0.0))

    with pytest.raises(TypeError, match="PocketCandidate"):
        measure_ligand_support(object(), ligand)  # type: ignore[arg-type]

    with pytest.raises(TypeError, match="StructureComponent"):
        measure_ligand_support(candidate, object())  # type: ignore[arg-type]


def test_measure_ligand_support_rejects_empty_ligands() -> None:
    candidate = PocketCandidate((_sphere(1, (0.0, 0.0, 0.0), 2.0),))
    ligand = StructureComponent(
        component_id="ATP-1",
        kind=ComponentKind.LIGAND,
        atoms=(),
        metadata={"selected_for_analysis": True},
        model_id="1",
        author_chain_id="A",
        residue_name="ATP",
        auth_seq_id="ATP-1",
    )

    with pytest.raises(ValueError, match="at least one atom"):
        measure_ligand_support(candidate, ligand)


def test_public_ligand_api_uses_one_primary_heavy_altloc_representation() -> None:
    candidate = PocketCandidate((_sphere(1, (0.0, 0.0, 0.0), 2.0),))
    ligand = StructureComponent(
        component_id="ATP-altloc",
        kind=ComponentKind.LIGAND,
        atoms=(
            AtomRecord(
                "C1",
                "C",
                (0.0, 0.0, 0.0),
                altloc="A",
                occupancy=0.8,
                source_atom_id="atp-c1-a",
            ),
            AtomRecord(
                "C1",
                "C",
                (4.0, 0.0, 0.0),
                altloc="B",
                occupancy=0.2,
                source_atom_id="atp-c1-b",
            ),
            AtomRecord("H1", "H", (0.0, 0.0, 0.0), source_atom_id="atp-h1"),
        ),
        metadata={"selected_for_analysis": True},
        model_id="1",
        author_chain_id="A",
        residue_name="ATP",
        auth_seq_id="ATP-altloc",
    )

    eligible = eligible_pocket_ligands((ligand,))
    support = measure_ligand_support(candidate, ligand)

    assert tuple(atom.source_atom_id for atom in eligible[0].atoms) == ("atp-c1-a",)
    assert support.atom_count == 1
    assert support.covered_atom_count == 1
    assert support.atom_coverage_fraction == 1.0
    assert support.to_json()["atom_scope"] == "primary_deposited_heavy_atoms"


def test_ligand_support_validates_contact_residue_types() -> None:
    candidate = PocketCandidate((_sphere(1, (0.0, 0.0, 0.0), 2.0),))
    ligand = _ligand("ATP-1", "ATP", (0.0, 0.0, 0.0))

    with pytest.raises(TypeError, match="ligand_contact_residues"):
        measure_ligand_support(
            candidate,
            ligand,
            ligand_contact_residues=(object(),),  # type: ignore[arg-type]
        )
