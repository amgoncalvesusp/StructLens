from __future__ import annotations

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
from structlens.core.pockets import PocketDetectionSettings


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


def test_pocket_service_detects_buried_pockets_and_records_provenance() -> None:
    progress: list[str] = []

    report = detect_blind_pockets(
        _parsed_shell(),
        progress_callback=progress.append,
    )

    assert isinstance(report, PocketDetectionReport)
    assert report.availability is Availability.AVAILABLE
    assert report.candidates
    assert report.counts["input_atoms"] > 0
    assert report.counts["alpha_spheres"] >= len(report.candidates[0].alpha_spheres)
    assert report.provenance is not None
    assert report.provenance.method_id == "structlens.pocket.detect"
    assert report.provenance.parameters["selection_id"] == _parsed_shell().selection.selection_id
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
    report = StructurePocketService(
        PocketDetectionSettings(solvent_grid_spacing_angstrom=0.75)
    ).analyze(_parsed_shell(open_face=True))

    assert report.availability is Availability.NOT_DETECTED
    assert report.candidates == ()
    assert any(item.code == "pocket.detect.solvent_exposed_cluster" for item in report.diagnostics)


def test_pocket_service_honors_cancellation_before_work_starts() -> None:
    cancel_event = Event()
    cancel_event.set()

    with pytest.raises(AnalysisCancelledError):
        detect_blind_pockets(_parsed_shell(), cancel_event=cancel_event)


def test_pocket_service_surfaces_resource_limit_as_a_typed_failure() -> None:
    report = StructurePocketService(
        PocketDetectionSettings(max_estimated_simplices=10)
    ).analyze(_parsed_shell())

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
