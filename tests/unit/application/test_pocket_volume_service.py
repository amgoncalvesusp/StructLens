from __future__ import annotations

from dataclasses import FrozenInstanceError, replace

import pytest

from structlens.application.pocket_comparison_service import PocketComparisonService
from structlens.application.pocket_service import StructurePocketService
from structlens.core.evidence import Availability
from structlens.core.models import (
    AtomRecord,
    ComponentKind,
    ProteinChain,
    ProteinStructure,
    ResidueCorrespondence,
    ResidueId,
    ResidueNumbering,
    ResidueRecord,
    StructureComponent,
)
from structlens.core.parsing import (
    AltlocPolicy,
    AssemblyScope,
    InputSelection,
    ParsedStructure,
    StructureFormat,
    StructureMetadata,
)
from structlens.core.pockets import (
    POCKET_RADII_VERSION,
    AlphaSphere,
    PocketCandidate,
    compare_pocket_volumes,
)
from structlens.core.pockets.volume import PocketVolumeSettings


def _candidate(*, radius: float = 2.0) -> PocketCandidate:
    residue = ResidueId("synthetic", "1", "A", "10", None, "ALA")
    atom_ids = tuple(f"sphere-atom-{index}" for index in range(4))
    selection = _selection()
    return PocketCandidate(
        (
            AlphaSphere(
                center_xyz=(0.0, 0.0, 0.0),
                radius_angstrom=radius,
                touching_atom_ids=atom_ids,
                lining_residues=(residue,),
                source_simplex_atom_ids=atom_ids,
            ),
        ),
        source_content_id=selection.content_id,
        selection_id=selection.selection_id,
    )


def _selection() -> InputSelection:
    return InputSelection(
        "a" * 64,
        "holo-pocket.pdb",
        StructureFormat.PDB,
        "1",
        author_chain_ids=("A",),
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
        _atom((0.0, 0.0, 0.0), element="O", source_atom_id="ligand-o"),
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
    selection = _selection()
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
        "candidate": candidate.candidate_id,
        "source_content": candidate.source_content_id,
    }
    assert first.provenance.parameters["radii_version"] == POCKET_RADII_VERSION
    assert first.provenance.parameters["component_exclusion_policy"] == "protein_only"
    assert first.provenance.parameters["grid_phase"] == "cell_center"
    assert first.provenance.parameters["selection_id"] == parsed.selection.selection_id
    assert first.provenance_parameters["raw_source_hash"] == parsed.raw_source_hash
    assert first.provenance_parameters["logical_content"] == parsed.selection.content_id
    assert first.provenance.artifact_id == second.provenance.artifact_id


@pytest.mark.parametrize(
    "tamper",
    ("method_raw_source", "result_raw_source", "method_logical_content", "result_logical_content"),
)
def test_comparison_rejects_isolated_source_snapshot_tampering_from_real_service(tamper: str) -> None:
    parsed = _parsed()
    candidate = _candidate()
    valid = StructurePocketService().measure_volume(parsed, candidate, settings=_settings())
    assert valid.provenance is not None
    if tamper.startswith("method_"):
        hash_name = tamper.removeprefix("method_")
        hashes = dict(valid.provenance.input_hashes)
        hashes[hash_name] = "f" * 64
        tampered = replace(valid, provenance=replace(valid.provenance, input_hashes=hashes))
    else:
        result_name = {
            "result_raw_source": "raw_source_hash",
            "result_logical_content": "logical_content",
        }[tamper]
        tampered = replace(
            valid,
            provenance_parameters={**valid.provenance_parameters, result_name: "f" * 64},
        )

    with pytest.raises(ValueError, match="source|logical|hash|provenance"):
        PocketComparisonService().compare(
            candidate,
            candidate,
            (),
            reference_volume=tampered,
            target_volume=valid,
        )


def test_structure_volume_service_results_are_accepted_by_comparison_service() -> None:
    parsed = _parsed()
    target_selection = replace(parsed.selection, content_id="c" * 64, display_name="target-pocket.pdb")
    target_parsed = replace(parsed, selection=target_selection, raw_source_hash="d" * 64)
    reference = _candidate(radius=2.0)
    target = replace(
        _candidate(radius=2.25),
        source_content_id=target_selection.content_id,
        selection_id=target_selection.selection_id,
    )
    volume_service = StructurePocketService()
    reference_volume = volume_service.measure_volume(parsed, reference, settings=_settings())
    target_volume = volume_service.measure_volume(target_parsed, target, settings=_settings())
    residue = reference.alpha_spheres[0].lining_residues[0]
    correspondence = ResidueCorrespondence(0, residue, residue, "A", "A", "conserved")

    report = PocketComparisonService().compare(
        reference,
        target,
        (correspondence,),
        reference_volume=reference_volume,
        target_volume=target_volume,
    )

    assert report.availability is Availability.AVAILABLE
    assert len(report.comparisons) == 1
    assert report.comparisons[0].volume_delta_angstrom3 is not None


@pytest.mark.parametrize("container", ("result", "method"))
@pytest.mark.parametrize("field", ("sphere_ids", "sphere_radii_angstrom", "sphere_count"))
def test_application_volume_boundary_rejects_each_isolated_sphere_metadata_tampering(
    container: str,
    field: str,
) -> None:
    parsed = _parsed()
    candidate = _candidate()
    valid = StructurePocketService().measure_volume(parsed, candidate, settings=_settings())
    assert valid.provenance is not None
    if field == "sphere_ids":
        tampered_value: object = ("f" * 64,)
    elif field == "sphere_radii_angstrom":
        tampered_value = (9.0,)
    else:
        tampered_value = 2
    if container == "result":
        tampered = replace(
            valid,
            provenance_parameters={**valid.provenance_parameters, field: tampered_value},
        )
    else:
        tampered = replace(
            valid,
            provenance=replace(
                valid.provenance,
                parameters={**valid.provenance.parameters, field: tampered_value},
            ),
        )

    with pytest.raises(ValueError, match="sphere|candidate|provenance|coherent"):
        PocketComparisonService().compare(
            candidate,
            candidate,
            (),
            reference_volume=tampered,
            target_volume=valid,
        )


@pytest.mark.parametrize("field", ("sphere_ids", "sphere_radii_angstrom", "sphere_count"))
def test_application_volume_boundary_checks_coordinated_sphere_metadata_against_candidate(field: str) -> None:
    parsed = _parsed()
    candidate = _candidate()
    valid = StructurePocketService().measure_volume(parsed, candidate, settings=_settings())
    assert valid.provenance is not None
    if field == "sphere_ids":
        tampered_value: object = ("f" * 64,)
    elif field == "sphere_radii_angstrom":
        tampered_value = (9.0,)
    else:
        tampered_value = 2
    tampered = replace(
        valid,
        provenance_parameters={**valid.provenance_parameters, field: tampered_value},
        provenance=replace(
            valid.provenance,
            parameters={**valid.provenance.parameters, field: tampered_value},
        ),
    )

    with pytest.raises(ValueError, match="sphere|candidate"):
        PocketComparisonService().compare(
            candidate,
            candidate,
            (),
            reference_volume=tampered,
            target_volume=valid,
        )


@pytest.mark.parametrize("container", ("result", "method"))
@pytest.mark.parametrize("field", ("candidate_id", "source_content_id", "selection_id"))
def test_application_volume_boundary_checks_each_identity_field_against_candidate(
    container: str,
    field: str,
) -> None:
    parsed = _parsed()
    candidate = _candidate()
    valid = StructurePocketService().measure_volume(parsed, candidate, settings=_settings())
    assert valid.provenance is not None
    if container == "result":
        tampered = replace(
            valid,
            provenance_parameters={**valid.provenance_parameters, field: "f" * 64},
        )
    else:
        tampered = replace(
            valid,
            provenance=replace(
                valid.provenance,
                parameters={**valid.provenance.parameters, field: "f" * 64},
            ),
        )

    with pytest.raises(ValueError, match="candidate|lineage|provenance"):
        PocketComparisonService().compare(
            candidate,
            candidate,
            (),
            reference_volume=tampered,
            target_volume=valid,
        )


@pytest.mark.parametrize(
    "required_parameter",
    (
        "selection",
        "model_id",
        "author_chain_ids",
        "label_chain_ids",
        "altloc_policy",
        "assembly_scope",
        "atom_scope",
        "component_scope",
        "component_rules_version",
        "hydrogen_policy",
        "polymer_atom_count",
        "component_ids",
    ),
)
def test_application_volume_boundary_requires_complete_application_schema(required_parameter: str) -> None:
    parsed = _parsed()
    candidate = _candidate()
    valid = StructurePocketService().measure_volume(parsed, candidate, settings=_settings())
    assert valid.provenance is not None
    parameters = dict(valid.provenance.parameters)
    parameters.pop(required_parameter)
    tampered = replace(valid, provenance=replace(valid.provenance, parameters=parameters))

    with pytest.raises(ValueError, match="application|schema|selection|policy|provenance"):
        PocketComparisonService().compare(
            candidate,
            candidate,
            (),
            reference_volume=tampered,
            target_volume=valid,
        )


@pytest.mark.parametrize("required_hash", ("raw_source", "logical_content", "candidate", "source_content"))
def test_application_volume_boundary_requires_complete_input_hash_schema(required_hash: str) -> None:
    parsed = _parsed()
    candidate = _candidate()
    valid = StructurePocketService().measure_volume(parsed, candidate, settings=_settings())
    assert valid.provenance is not None
    hashes = dict(valid.provenance.input_hashes)
    hashes.pop(required_hash)
    tampered = replace(valid, provenance=replace(valid.provenance, input_hashes=hashes))

    with pytest.raises(ValueError, match="hash|candidate|source|schema|provenance"):
        PocketComparisonService().compare(
            candidate,
            candidate,
            (),
            reference_volume=tampered,
            target_volume=valid,
        )


def test_application_volume_rejects_core_alias_with_application_only_schema() -> None:
    parsed = _parsed()
    candidate = _candidate()
    valid = StructurePocketService().measure_volume(parsed, candidate, settings=_settings())
    assert valid.provenance is not None
    aliased = replace(
        valid,
        provenance=replace(
            valid.provenance,
            method_id="structlens.pocket_free_volume",
            method_version="0.4",
        ),
    )

    with pytest.raises(ValueError, match="core|schema|producer|hash|provenance"):
        PocketComparisonService().compare(
            candidate,
            candidate,
            (),
            reference_volume=aliased,
            target_volume=valid,
        )


@pytest.mark.parametrize(
    "tamper",
    (
        "selection_not_mapping",
        "selection_missing_field",
        "selection_mismatch",
        "selection_format",
        "altloc_policy",
        "assembly_scope",
        "analyzed_representation",
        "author_chain_ids",
        "chain_locators",
        "resource_setting",
        "atom_scope",
        "polymer_atom_count",
        "component_ids",
        "logical_content_hash",
        "backend_schema",
        "units_schema",
    ),
)
def test_application_volume_strict_schema_rejects_incoherent_producer_fields(tamper: str) -> None:
    parsed = _parsed()
    candidate = _candidate()
    valid = StructurePocketService().measure_volume(parsed, candidate, settings=_settings())
    assert valid.provenance is not None
    method_parameters = dict(valid.provenance.parameters)
    result_parameters = dict(valid.provenance_parameters)
    provenance_kwargs: dict[str, object] = {}
    signature = valid.compatibility_signature
    if tamper == "selection_not_mapping":
        method_parameters["selection"] = "invalid"
    else:
        selection = dict(method_parameters["selection"])  # type: ignore[arg-type]
        if tamper == "selection_missing_field":
            selection.pop("format")
        elif tamper == "selection_mismatch":
            selection["model_id"] = "2"
        elif tamper == "selection_format":
            selection["format"] = "mol2"
        elif tamper == "altloc_policy":
            selection["altloc_policy"] = method_parameters["altloc_policy"] = "invalid"
            result_parameters["altloc_policy"] = "invalid"
            signature = signature[:-2] + (("altloc_policy", "invalid"), signature[-1])
        elif tamper == "assembly_scope":
            selection["assembly_scope"] = method_parameters["assembly_scope"] = "invalid"
            provenance_kwargs["analyzed_representation"] = "invalid"
        elif tamper == "author_chain_ids":
            selection["author_chain_ids"] = method_parameters["author_chain_ids"] = (1,)
        elif tamper == "chain_locators":
            selection["chain_locators"] = ({"author_chain_id": "A"},)
        method_parameters["selection"] = selection
    if tamper == "analyzed_representation":
        provenance_kwargs["analyzed_representation"] = "biological_assembly"
    elif tamper == "resource_setting":
        method_parameters["max_voxel_count"] = valid.settings.max_voxel_count + 1  # type: ignore[union-attr]
    elif tamper == "atom_scope":
        method_parameters["atom_scope"] = "all_atoms"
    elif tamper == "polymer_atom_count":
        method_parameters["polymer_atom_count"] = -1
    elif tamper == "component_ids":
        method_parameters["component_ids"] = (1,)
    elif tamper == "logical_content_hash":
        hashes = dict(valid.provenance.input_hashes)
        hashes["logical_content"] = "f" * 64
        provenance_kwargs["input_hashes"] = hashes
    elif tamper == "backend_schema":
        provenance_kwargs["backend_versions"] = {"numpy": "1", "scipy": "1"}
    elif tamper == "units_schema":
        provenance_kwargs["units"] = {**valid.provenance.units, "extra": "count"}
    tampered = replace(
        valid,
        provenance_parameters=result_parameters,
        compatibility_signature=signature,
        provenance=replace(valid.provenance, parameters=method_parameters, **provenance_kwargs),
    )

    with pytest.raises(ValueError, match="application|selection|policy|schema|setting|count|hash|unit|_ids"):
        PocketComparisonService().compare(
            candidate,
            candidate,
            (),
            reference_volume=tampered,
            target_volume=valid,
        )


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
        "O",
        "O",
        (0.0, 0.0, 0.0),
        altloc="A",
        occupancy=0.8,
        source_atom_id="ligand-o-a",
    )
    alternate = AtomRecord(
        "O",
        "O",
        (1.5, 0.0, 0.0),
        altloc="B",
        occupancy=0.2,
        source_atom_id="ligand-o-b",
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


def test_service_disables_cross_altloc_policy_volume_deltas() -> None:
    parsed = _parsed()
    first_choice = AtomRecord(
        "O",
        "O",
        (1.5, 0.0, 0.0),
        altloc="B",
        occupancy=0.2,
        source_atom_id="ligand-o-b",
    )
    highest_choice = AtomRecord(
        "O",
        "O",
        (0.0, 0.0, 0.0),
        altloc="A",
        occupancy=0.8,
        source_atom_id="ligand-o-a",
    )
    alternate_ligand = replace(
        next(component for component in parsed.components if component.component_id == "ligand-1"),
        atoms=(first_choice, highest_choice),
    )
    alternate_components = tuple(
        alternate_ligand if component.component_id == "ligand-1" else component for component in parsed.components
    )
    highest_parsed = ParsedStructure(
        parsed.protein_structure,
        alternate_components,
        parsed.selection,
        parsed.metadata,
        parsed.raw_source_hash,
    )
    first_selection = replace(highest_parsed.selection, altloc_policy=AltlocPolicy.FIRST)
    first_parsed = ParsedStructure(
        highest_parsed.protein_structure,
        alternate_components,
        first_selection,
        highest_parsed.metadata,
        highest_parsed.raw_source_hash,
    )
    settings = _settings(policy="unoccupied")
    service = StructurePocketService()

    highest = service.measure_volume(highest_parsed, _candidate(), settings=settings)
    first_candidate = replace(
        _candidate(),
        source_content_id=first_parsed.selection.content_id,
        selection_id=first_parsed.selection.selection_id,
    )
    first = service.measure_volume(first_parsed, first_candidate, settings=settings)
    comparison = compare_pocket_volumes(highest, first)

    assert highest.fine_volume_angstrom3 != first.fine_volume_angstrom3
    assert highest.compatibility_signature != first.compatibility_signature
    assert comparison.availability is Availability.NOT_APPLICABLE
    assert comparison.delta_angstrom3 is None


def test_service_uses_heavy_atoms_for_both_polymer_and_retained_components() -> None:
    parsed = _parsed()
    hydrogen_ligand = _component(
        "ligand-1",
        ComponentKind.LIGAND,
        _atom((0.0, 0.0, 0.0), element="H", source_atom_id="ligand-h"),
    )
    hydrogen_only = ParsedStructure(
        parsed.protein_structure,
        tuple(
            hydrogen_ligand if component.component_id == "ligand-1" else component for component in parsed.components
        ),
        parsed.selection,
        parsed.metadata,
        parsed.raw_source_hash,
    )
    service = StructurePocketService()

    protein_only = service.measure_volume(
        hydrogen_only,
        _candidate(),
        settings=_settings(policy="protein_only"),
    )
    unoccupied = service.measure_volume(
        hydrogen_only,
        _candidate(),
        settings=_settings(policy="unoccupied"),
    )

    assert unoccupied.coarse_voxel_count == protein_only.coarse_voxel_count
    assert unoccupied.fine_voxel_count == protein_only.fine_voxel_count
    assert unoccupied.provenance.parameters["hydrogen_policy"] == "deposited_heavy_atoms_only"


def test_service_preserves_unavailable_state_for_missing_candidate() -> None:
    report = StructurePocketService().measure_volume(_parsed(), None, settings=_settings())

    assert report.availability is Availability.NOT_APPLICABLE
    assert report.coarse_volume_angstrom3 is None
    assert report.fine_volume_angstrom3 is None
    assert report.provenance is not None
    assert any(item.code == "pocket.volume.no_candidate" for item in report.diagnostics)


def test_volume_service_rejects_candidate_from_another_selection() -> None:
    parsed = _parsed()
    foreign = replace(_candidate(), source_content_id="c" * 64)

    with pytest.raises(ValueError, match="candidate lineage"):
        StructurePocketService().measure_volume(parsed, foreign, settings=_settings())


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
