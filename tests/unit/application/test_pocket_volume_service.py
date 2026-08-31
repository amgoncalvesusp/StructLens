from __future__ import annotations

from dataclasses import FrozenInstanceError

import pytest

from structlens.application.pocket_service import StructurePocketService
from structlens.core.evidence import Availability
from structlens.core.models import (
    AtomRecord,
    ComponentKind,
    ProteinChain,
    ProteinStructure,
    ResidueId,
    ResidueNumbering,
    ResidueRecord,
    StructureComponent,
)
from structlens.core.parsing import (
    AssemblyScope,
    InputSelection,
    ParsedStructure,
    StructureFormat,
    StructureMetadata,
)
from structlens.core.pockets import POCKET_RADII_VERSION, AlphaSphere, PocketCandidate
from structlens.core.pockets.volume import PocketVolumeSettings


def _candidate(*, radius: float = 2.0) -> PocketCandidate:
    residue = ResidueId("synthetic", "1", "A", "10", None, "ALA")
    atom_ids = tuple(f"sphere-atom-{index}" for index in range(4))
    return PocketCandidate(
        (
            AlphaSphere(
                center_xyz=(0.0, 0.0, 0.0),
                radius_angstrom=radius,
                touching_atom_ids=atom_ids,
                lining_residues=(residue,),
                source_simplex_atom_ids=atom_ids,
            ),
        )
    )


def _atom(
    coordinate: tuple[float, float, float],
    *,
    element: str = "H",
    source_atom_id: str,
) -> AtomRecord:
    return AtomRecord(
        name=element,
        element=element,
        coordinate=coordinate,
        source_atom_id=source_atom_id,
    )


def _component(
    component_id: str,
    kind: ComponentKind,
    atom: AtomRecord,
    *,
    selected: bool = True,
) -> StructureComponent:
    return StructureComponent(
        component_id,
        kind,
        (atom,),
        metadata={"selected_for_analysis": selected},
        model_id="1",
        author_chain_id="A",
        residue_name="LIG" if kind is ComponentKind.LIGAND else "HOH",
        auth_seq_id=component_id,
    )


def _parsed(*, include_unselected_ligand: bool = False) -> ParsedStructure:
    # Keep the polymer atom just outside the synthetic cavity's sampled cells.
    # The retained ligand is at the center and is used to distinguish the two
    # explicit component policies in a holo structure.
    protein_atom = _atom((0.0, 0.0, 3.0), source_atom_id="protein-h")
    residue_id = ResidueId("protein", "1", "A", "1", None, "ALA")
    residue = ResidueRecord(
        residue_id,
        ResidueNumbering("1", "1", None),
        "ALA",
        "A",
        (protein_atom,),
    )
    chain = ProteinChain(
        "protein",
        "1",
        "A",
        residues=(residue_id,),
        sequence="A",
        residue_records=(residue,),
        author_chain_id="A",
    )
    polymer = StructureComponent(
        "polymer-1",
        ComponentKind.POLYMER_RESIDUE,
        (protein_atom,),
        residue_id=residue_id,
        model_id="1",
        author_chain_id="A",
        residue_name="ALA",
        auth_seq_id="1",
    )
    ligand = _component(
        "ligand-1",
        ComponentKind.LIGAND,
        _atom((0.0, 0.0, 0.0), source_atom_id="ligand-h"),
    )
    water = _component(
        "water-1",
        ComponentKind.WATER,
        _atom((0.0, 0.0, 0.0), source_atom_id="water-h"),
    )
    components: tuple[StructureComponent, ...] = (polymer, ligand, water)
    if include_unselected_ligand:
        unselected = _component(
            "ligand-unselected",
            ComponentKind.LIGAND,
            _atom((1.5, 1.5, 1.5), source_atom_id="unselected-h"),
            selected=False,
        )
        components += (unselected,)
    selection = InputSelection(
        "a" * 64,
        "holo-pocket.pdb",
        StructureFormat.PDB,
        "1",
        author_chain_ids=("A",),
    )
    metadata = StructureMetadata(
        StructureFormat.PDB,
        "1",
        ("1",),
        ("A",),
        assembly_scope=AssemblyScope.ASYMMETRIC_UNIT,
        author_chain_ids=("A",),
    )
    return ParsedStructure(
        ProteinStructure("protein", (chain,)),
        components,
        selection,
        metadata,
        "b" * 64,
    )


def _settings(
    *,
    policy: str = "protein_only",
    max_voxels: int = 100_000,
) -> PocketVolumeSettings:
    return PocketVolumeSettings(
        coarse_grid_spacing_angstrom=1.0,
        fine_grid_spacing_angstrom=0.5,
        boundary_margin_angstrom=0.0,
        component_exclusion_policy=policy,
        max_voxel_count=max_voxels,
        voxel_chunk_size=8,
    )


def test_volume_service_returns_typed_result_with_complete_provenance() -> None:
    parsed = _parsed()
    candidate = _candidate()
    service = StructurePocketService()

    first = service.measure_volume(parsed, candidate, settings=_settings())
    second = service.measure_volume(parsed, candidate, settings=_settings())

    assert first.availability is Availability.AVAILABLE
    assert first == second
    assert first.units["volume"] == "angstrom^3"
    assert first.provenance is not None
    assert first.provenance.method_id == "structlens.pocket.volume"
    assert first.provenance.input_hashes == {
        "raw_source": "b" * 64,
        "logical_content": "a" * 64,
    }
    assert first.provenance.parameters["radii_version"] == POCKET_RADII_VERSION
    assert first.provenance.parameters["component_exclusion_policy"] == "protein_only"
    assert first.provenance.parameters["grid_phase"] == "cell_center"
    assert first.provenance.parameters["selection_id"] == parsed.selection.selection_id
    assert first.provenance.artifact_id == second.provenance.artifact_id


def test_service_applies_protein_only_and_unoccupied_policies_to_holo_components() -> None:
    parsed = _parsed()
    candidate = _candidate()
    service = StructurePocketService()

    protein_only = service.measure_volume(parsed, candidate, settings=_settings(policy="protein_only"))
    unoccupied = service.measure_volume(parsed, candidate, settings=_settings(policy="unoccupied"))

    # The polymer atom is outside the sampled sphere cells. The ligand at the
    # center removes the eight innermost fine-grid cells only in unoccupied.
    assert protein_only.coarse_voxel_count == 32
    assert unoccupied.coarse_voxel_count == 24
    assert unoccupied.fine_volume_angstrom3 < protein_only.fine_volume_angstrom3
    assert protein_only.provenance.parameters["component_exclusion_policy"] == "protein_only"
    assert unoccupied.provenance.parameters["component_exclusion_policy"] == "unoccupied"


def test_service_ignores_unselected_ligand_and_water_for_unoccupied_policy() -> None:
    candidate = _candidate()
    service = StructurePocketService()
    selected_only = _parsed()
    with_unselected = _parsed(include_unselected_ligand=True)

    baseline = service.measure_volume(selected_only, candidate, settings=_settings(policy="unoccupied"))
    report = service.measure_volume(with_unselected, candidate, settings=_settings(policy="unoccupied"))

    assert report.coarse_voxel_count == baseline.coarse_voxel_count
    assert report.fine_voxel_count == baseline.fine_voxel_count
    assert report.fine_volume_angstrom3 == baseline.fine_volume_angstrom3


def test_service_applies_the_selected_altloc_policy_to_retained_ligands() -> None:
    parsed = _parsed()
    primary = AtomRecord(
        "H",
        "H",
        (0.0, 0.0, 0.0),
        altloc="A",
        occupancy=0.8,
        source_atom_id="ligand-h-a",
    )
    alternate = AtomRecord(
        "H",
        "H",
        (1.5, 0.0, 0.0),
        altloc="B",
        occupancy=0.2,
        source_atom_id="ligand-h-b",
    )
    ligand_with_altlocs = StructureComponent(
        "ligand-1",
        ComponentKind.LIGAND,
        (primary, alternate),
        metadata={"selected_for_analysis": True},
        model_id="1",
        author_chain_id="A",
        residue_name="LIG",
        auth_seq_id="ligand-1",
    )
    with_altlocs = ParsedStructure(
        parsed.protein_structure,
        tuple(
            ligand_with_altlocs if component.component_id == "ligand-1" else component
            for component in parsed.components
        ),
        parsed.selection,
        parsed.metadata,
        parsed.raw_source_hash,
    )
    service = StructurePocketService()

    baseline = service.measure_volume(parsed, _candidate(), settings=_settings(policy="unoccupied"))
    report = service.measure_volume(
        with_altlocs,
        _candidate(),
        settings=_settings(policy="unoccupied"),
    )

    assert report.coarse_voxel_count == baseline.coarse_voxel_count
    assert report.fine_voxel_count == baseline.fine_voxel_count


def test_service_preserves_unavailable_state_for_missing_candidate() -> None:
    report = StructurePocketService().measure_volume(_parsed(), None, settings=_settings())

    assert report.availability is Availability.NOT_APPLICABLE
    assert report.coarse_volume_angstrom3 is None
    assert report.fine_volume_angstrom3 is None
    assert report.provenance is not None
    assert any(item.code == "pocket.volume.no_candidate" for item in report.diagnostics)


def test_service_propagates_voxel_limit_as_typed_failure_with_remediation() -> None:
    report = StructurePocketService().measure_volume(
        _parsed(),
        _candidate(),
        settings=_settings(max_voxels=4),
    )

    assert report.availability is Availability.NUMERICAL_FAILURE
    assert report.coarse_volume_angstrom3 is None
    diagnostic = next(item for item in report.diagnostics if item.code == "pocket.volume.resource_limit")
    assert diagnostic.remediation
    assert "voxel" in diagnostic.message.lower()


@pytest.mark.parametrize("policy", ("protein_only", "unoccupied"))
def test_service_result_is_immutable_for_each_component_policy(policy: str) -> None:
    report = StructurePocketService().measure_volume(
        _parsed(),
        _candidate(),
        settings=_settings(policy=policy),
    )

    with pytest.raises(FrozenInstanceError):
        report.fine_volume_angstrom3 = 1.0  # type: ignore[misc]


@pytest.mark.parametrize(
    "bad_policy",
    ("all_atoms", "ligand_only", "protein-plus-ligand"),
)
def test_service_rejects_unknown_component_policy_at_the_boundary(bad_policy: str) -> None:
    with pytest.raises(ValueError, match="component_exclusion_policy"):
        _settings(policy=bad_policy)
