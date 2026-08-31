from __future__ import annotations

from dataclasses import replace
from threading import Event

import pytest

from structlens.application.pocket_service import (
    PocketDetectionReport,
    StructurePocketService,
    detect_blind_pockets,
)
from structlens.core.errors import AnalysisCancelledError
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
from structlens.core.pockets import (
    POCKET_LIGAND_RULES_VERSION,
    AlphaSphere,
    PocketCandidate,
    PocketDetectionSettings,
)
from structlens.core.sites import SiteDefinition, SiteDefinitionMode


def _parsed_shell(*, open_face: bool = False) -> ParsedStructure:
    residues: list[ResidueRecord] = []
    components: list[StructureComponent] = []
    auth_seq_id = 1
    radius = 5.2
    coordinates = (-radius, -radius / 2.0, 0.0, radius / 2.0, radius)
    for x in coordinates:
        for y in coordinates:
            for z in coordinates:
                if max(abs(x), abs(y), abs(z)) != radius:
                    continue
                if open_face and x == radius:
                    continue
                residue_id = ResidueId("protein", "1", "A", str(auth_seq_id), None, "ALA")
                atom = AtomRecord("S", "S", (x, y, z), source_atom_id=f"atom-{auth_seq_id}")
                record = ResidueRecord(
                    residue_id,
                    ResidueNumbering(str(auth_seq_id), str(auth_seq_id), None),
                    "ALA",
                    "A",
                    (atom,),
                )
                residues.append(record)
                components.append(
                    StructureComponent(
                        f"component-{auth_seq_id}",
                        ComponentKind.POLYMER_RESIDUE,
                        (atom,),
                        residue_id=residue_id,
                        model_id="1",
                        author_chain_id="A",
                        residue_name="ALA",
                        auth_seq_id=str(auth_seq_id),
                    )
                )
                auth_seq_id += 1
    chain = ProteinChain(
        "protein",
        "1",
        "A",
        residues=tuple(item.residue_id for item in residues),
        sequence="A" * len(residues),
        residue_records=tuple(residues),
        author_chain_id="A",
    )
    selection = InputSelection(
        "a" * 64,
        "synthetic-shell.pdb",
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
        tuple(components),
        selection,
        metadata,
        "b" * 64,
    )


def _bind_candidate(candidate: PocketCandidate, parsed: ParsedStructure) -> PocketCandidate:
    return replace(
        candidate,
        source_content_id=parsed.selection.content_id,
        selection_id=parsed.selection.selection_id,
    )


def test_pocket_service_detects_buried_pockets_and_records_provenance() -> None:
    progress: list[str] = []
    parsed = _parsed_shell()

    report = detect_blind_pockets(
        parsed,
        progress_callback=progress.append,
    )

    assert isinstance(report, PocketDetectionReport)
    assert report.availability is Availability.AVAILABLE
    assert report.candidates
    assert report.counts["input_atoms"] > 0
    assert report.counts["alpha_spheres"] >= len(report.candidates[0].alpha_spheres)
    assert report.provenance is not None
    assert report.provenance.method_id == "structlens.pocket.detect"
    assert report.provenance.parameters["selection_id"] == parsed.selection.selection_id
    assert report.candidates[0].source_content_id == parsed.selection.content_id
    assert report.candidates[0].selection_id == parsed.selection.selection_id
    assert report.provenance.parameters["selection"] == {
        "model_id": "1",
        "author_chain_ids": ("A",),
        "label_chain_ids": (),
        "chain_locators": ({"author_chain_id": "A", "label_chain_id": None, "entity_id": None},),
        "altloc_policy": "highest_occupancy",
        "assembly_scope": "asymmetric_unit",
    }
    assert progress == ["prepare", "tessellate", "cluster", "rank"]


def test_pocket_service_reports_not_detected_for_open_cavity() -> None:
    report = StructurePocketService(PocketDetectionSettings(solvent_grid_spacing_angstrom=0.75)).analyze(
        _parsed_shell(open_face=True)
    )

    assert report.availability is Availability.NOT_DETECTED
    assert report.candidates == ()
    assert any(item.code == "pocket.detect.solvent_exposed_cluster" for item in report.diagnostics)


def test_pocket_service_honors_cancellation_before_work_starts() -> None:
    cancel_event = Event()
    cancel_event.set()

    with pytest.raises(AnalysisCancelledError):
        detect_blind_pockets(_parsed_shell(), cancel_event=cancel_event)


def test_pocket_service_surfaces_resource_limit_as_a_typed_failure() -> None:
    report = StructurePocketService(PocketDetectionSettings(max_estimated_simplices=10)).analyze(_parsed_shell())

    assert report.availability is Availability.NUMERICAL_FAILURE
    assert report.candidates == ()
    assert [item.code for item in report.diagnostics if item.code == "pocket.detect.resource_limit"] == [
        "pocket.detect.resource_limit"
    ]


def test_pocket_service_uses_only_normalized_selected_residue_atoms() -> None:
    parsed = _parsed_shell()
    unselected_residue = ResidueId("protein", "1", "B", "999", None, "ALA")
    alternate = AtomRecord("S", "S", (100.0, 100.0, 100.0), source_atom_id="unselected-atom")
    unselected_component = StructureComponent(
        "unselected-component",
        ComponentKind.POLYMER_RESIDUE,
        (alternate,),
        residue_id=unselected_residue,
        model_id="1",
        author_chain_id="B",
        residue_name="ALA",
        auth_seq_id="999",
        metadata={"selected_for_analysis": False},
    )
    with_extra_component = ParsedStructure(
        parsed.protein_structure,
        parsed.components + (unselected_component,),
        parsed.selection,
        parsed.metadata,
        parsed.raw_source_hash,
    )

    baseline = detect_blind_pockets(parsed)
    report = detect_blind_pockets(with_extra_component)

    assert report.counts["input_atoms"] == baseline.counts["input_atoms"]
    assert tuple(item.candidate_id for item in report.candidates) == tuple(
        item.candidate_id for item in baseline.candidates
    )


def test_pocket_service_does_not_widen_an_explicit_chain_selection() -> None:
    parsed = _parsed_shell()
    atom = AtomRecord("S", "S", (100.0, 100.0, 100.0), source_atom_id="chain-b-atom")
    residue_id = ResidueId("protein", "1", "B", "1", None, "ALA")
    residue = ResidueRecord(
        residue_id,
        ResidueNumbering("1", "1", None),
        "ALA",
        "A",
        (atom,),
    )
    chain_b = ProteinChain(
        "protein",
        "1",
        "B",
        residues=(residue_id,),
        sequence="A",
        residue_records=(residue,),
        author_chain_id="B",
    )
    broad_metadata = StructureMetadata(
        StructureFormat.PDB,
        "1",
        ("1",),
        ("A", "B"),
        assembly_scope=AssemblyScope.ASYMMETRIC_UNIT,
        author_chain_ids=("A", "B"),
    )
    broad_structure = ParsedStructure(
        ProteinStructure("protein", parsed.protein_structure.chains + (chain_b,)),
        parsed.components,
        parsed.selection,
        broad_metadata,
        parsed.raw_source_hash,
    )

    baseline = detect_blind_pockets(parsed)
    report = detect_blind_pockets(broad_structure)

    assert report.counts["input_atoms"] == baseline.counts["input_atoms"]
    assert tuple(item.candidate_id for item in report.candidates) == tuple(
        item.candidate_id for item in baseline.candidates
    )
    selection_b = InputSelection(
        parsed.selection.content_id,
        parsed.selection.display_name,
        parsed.selection.format,
        parsed.selection.model_id,
        author_chain_ids=("B",),
    )
    chain_b_report = detect_blind_pockets(
        ParsedStructure(
            broad_structure.protein_structure,
            broad_structure.components,
            selection_b,
            broad_metadata,
            broad_structure.raw_source_hash,
        )
    )
    assert report.provenance is not None
    assert chain_b_report.provenance is not None
    assert report.provenance.artifact_id != chain_b_report.provenance.artifact_id


def test_pocket_service_reports_fewer_than_four_atoms_as_invalid_input() -> None:
    parsed = _parsed_shell()
    residues = parsed.protein_structure.chains[0].residue_records[:3]
    chain = ProteinChain(
        "protein",
        "1",
        "A",
        residues=tuple(item.residue_id for item in residues),
        sequence="AAA",
        residue_records=residues,
        author_chain_id="A",
    )
    underspecified = ParsedStructure(
        ProteinStructure("protein", (chain,)),
        parsed.components[:3],
        parsed.selection,
        parsed.metadata,
        parsed.raw_source_hash,
    )

    report = detect_blind_pockets(underspecified)

    assert report.availability is Availability.INVALID_INPUT
    assert [item.code for item in report.diagnostics] == ["pocket.detect.insufficient_atoms"]


def test_pocket_service_excludes_unknown_radii_with_a_typed_diagnostic() -> None:
    parsed = _parsed_shell()
    records = parsed.protein_structure.chains[0].residue_records
    first = records[0]
    unknown_atom = AtomRecord(
        "XX",
        "XX",
        first.atoms[0].coordinate,
        source_atom_id=first.atoms[0].source_atom_id,
    )
    first_with_unknown = ResidueRecord(
        first.residue_id,
        first.numbering,
        first.residue_name,
        first.one_letter,
        (unknown_atom,),
    )
    selected_records = (first_with_unknown,) + records[1:]
    chain = ProteinChain(
        "protein",
        "1",
        "A",
        residues=tuple(item.residue_id for item in selected_records),
        sequence="A" * len(selected_records),
        residue_records=selected_records,
        author_chain_id="A",
    )
    with_unknown = ParsedStructure(
        ProteinStructure("protein", (chain,)),
        parsed.components,
        parsed.selection,
        parsed.metadata,
        parsed.raw_source_hash,
    )

    report = detect_blind_pockets(with_unknown)

    assert report.counts["input_atoms"] == len(selected_records) - 1
    assert report.counts["excluded_unknown_radii"] == 1
    assert [item.code for item in report.diagnostics].count("pocket.detect.unknown_radius") == 1


def test_pocket_service_honors_cancellation_between_stages() -> None:
    cancel_event = Event()

    def cancel_after_prepare(stage: str) -> None:
        if stage == "prepare":
            cancel_event.set()

    with pytest.raises(AnalysisCancelledError):
        detect_blind_pockets(
            _parsed_shell(),
            progress_callback=cancel_after_prepare,
            cancel_event=cancel_event,
        )


def test_pocket_service_selects_the_best_candidate_for_a_key_residue_site() -> None:
    parsed = _parsed_shell()
    first = PocketCandidate(
        (
            AlphaSphere(
                center_xyz=(0.0, 0.0, 0.0),
                radius_angstrom=2.0,
                touching_atom_ids=("a1", "a2", "a3", "a4"),
                lining_residues=(parsed.protein_structure.chains[0].residue_records[0].residue_id,),
                source_simplex_atom_ids=("a1", "a2", "a3", "a4"),
            ),
        )
    )
    second = PocketCandidate(
        (
            AlphaSphere(
                center_xyz=(4.0, 0.0, 0.0),
                radius_angstrom=2.0,
                touching_atom_ids=("b1", "b2", "b3", "b4"),
                lining_residues=(parsed.protein_structure.chains[0].residue_records[1].residue_id,),
                source_simplex_atom_ids=("b1", "b2", "b3", "b4"),
            ),
        )
    )
    definition = SiteDefinition(
        "focus",
        "Focus",
        SiteDefinitionMode.KEY_RESIDUES,
        (first.lining_residues[0],),
    )

    selection = StructurePocketService().select_focused_candidate(
        parsed,
        (_bind_candidate(second, parsed), _bind_candidate(first, parsed)),
        definition,
    )

    assert selection.availability is Availability.AVAILABLE
    assert selection.candidate == _bind_candidate(first, parsed)


def test_pocket_service_rejects_foreign_lineage_in_focused_selection() -> None:
    parsed = _parsed_shell()
    residue_id = parsed.protein_structure.chains[0].residue_records[0].residue_id
    candidate = _bind_candidate(
        PocketCandidate(
            (
                AlphaSphere(
                    center_xyz=(0.0, 0.0, 0.0),
                    radius_angstrom=2.0,
                    touching_atom_ids=("a1", "a2", "a3", "a4"),
                    lining_residues=(residue_id,),
                    source_simplex_atom_ids=("a1", "a2", "a3", "a4"),
                ),
            )
        ),
        parsed,
    )
    definition = SiteDefinition(
        "focus",
        "Focus",
        SiteDefinitionMode.KEY_RESIDUES,
        (residue_id,),
    )

    with pytest.raises(ValueError, match="candidate lineage"):
        StructurePocketService().select_focused_candidate(
            parsed,
            (replace(candidate, source_content_id="c" * 64),),
            definition,
        )


def test_pocket_service_reports_no_eligible_ligand_for_buffer_only_sites() -> None:
    parsed = _parsed_shell()
    sulfate = StructureComponent(
        component_id="SO4-1",
        kind=ComponentKind.LIGAND,
        atoms=(AtomRecord("S", "S", (0.0, 0.0, 0.0), source_atom_id="so4-s"),),
        metadata={"selected_for_analysis": True},
        model_id="1",
        author_chain_id="A",
        residue_name="SO4",
        auth_seq_id="SO4-1",
    )
    with_buffer = ParsedStructure(
        parsed.protein_structure,
        parsed.components + (sulfate,),
        parsed.selection,
        parsed.metadata,
        parsed.raw_source_hash,
    )
    candidate = PocketCandidate(
        (
            AlphaSphere(
                center_xyz=(0.0, 0.0, 0.0),
                radius_angstrom=2.0,
                touching_atom_ids=("c1", "c2", "c3", "c4"),
                lining_residues=(parsed.protein_structure.chains[0].residue_records[0].residue_id,),
                source_simplex_atom_ids=("c1", "c2", "c3", "c4"),
            ),
        )
    )
    definition = SiteDefinition(
        "ligand-focus",
        "Ligand Focus",
        SiteDefinitionMode.LIGAND_RADIUS,
        (),
        None,
        "SO4-1",
        4.0,
    )

    selection = StructurePocketService().select_focused_candidate(
        with_buffer,
        (_bind_candidate(candidate, with_buffer),),
        definition,
    )

    assert selection.availability is Availability.NOT_APPLICABLE
    assert selection.candidate is None
    assert any(item.code == "pocket.focus.no_eligible_ligand" for item in selection.diagnostics)


def test_pocket_service_retains_ligand_support_and_focus_provenance() -> None:
    parsed = _parsed_shell()
    atp = StructureComponent(
        component_id="ATP-1",
        kind=ComponentKind.LIGAND,
        atoms=(AtomRecord("P", "P", (0.0, 0.0, 0.0), source_atom_id="atp-p"),),
        metadata={"selected_for_analysis": True},
        model_id="1",
        author_chain_id="A",
        residue_name="ATP",
        auth_seq_id="ATP-1",
    )
    with_atp = ParsedStructure(
        parsed.protein_structure,
        parsed.components + (atp,),
        parsed.selection,
        parsed.metadata,
        parsed.raw_source_hash,
    )
    candidate = PocketCandidate(
        (
            AlphaSphere(
                center_xyz=(0.0, 0.0, 0.0),
                radius_angstrom=2.0,
                touching_atom_ids=("d1", "d2", "d3", "d4"),
                lining_residues=(parsed.protein_structure.chains[0].residue_records[0].residue_id,),
                source_simplex_atom_ids=("d1", "d2", "d3", "d4"),
            ),
        )
    )
    definition = SiteDefinition(
        "atp-focus",
        "ATP focus",
        SiteDefinitionMode.LIGAND_RADIUS,
        ligand_id="ATP-1",
        radius_angstrom=4.0,
    )

    selection = StructurePocketService().select_focused_candidate(
        with_atp,
        (_bind_candidate(candidate, with_atp),),
        definition,
    )

    assert selection.availability is Availability.AVAILABLE
    assert selection.ligand_support is not None
    assert selection.ligand_support.component_id == "ATP-1"
    assert selection.provenance is not None
    assert selection.provenance.method_id == "structlens.pocket.focus"
    assert selection.provenance.parameters["ligand_rules_version"] == POCKET_LIGAND_RULES_VERSION
    assert selection.to_json()["provenance"]["input_hashes"]["raw_source"] == "b" * 64
