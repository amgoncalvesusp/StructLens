from __future__ import annotations

from structlens.application.quality_service import StructureQualityService, assess_structure_quality
from structlens.core.evidence import Availability, Diagnostic, DiagnosticSeverity
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


def _residue(number: int, origin: float, *, overlap_sidechain: bool = False) -> ResidueRecord:
    residue_id = ResidueId("protein", "1", "A", str(number), None, "ALA")
    atoms = (
        AtomRecord("N", "N", (origin - 1.3, 0.0, 0.0), source_atom_id=f"{number}-n"),
        AtomRecord("CA", "C", (origin, 0.0, 0.0), source_atom_id=f"{number}-ca"),
        AtomRecord("C", "C", (origin + 1.3, 0.0, 0.0), source_atom_id=f"{number}-c"),
        AtomRecord("O", "O", (origin + 1.3, 1.2, 0.0), source_atom_id=f"{number}-o"),
        AtomRecord(
            "CB",
            "C",
            (0.0, 3.0, 0.0) if overlap_sidechain else (origin, 1.5, 0.0),
            source_atom_id=f"{number}-cb",
        ),
    )
    return ResidueRecord(
        residue_id,
        ResidueNumbering(str(number), str(number), None),
        "ALA",
        "A",
        atoms,
    )


def _parsed_with_overlap() -> ParsedStructure:
    residues = (
        _residue(1, 0.0, overlap_sidechain=True),
        _residue(2, 3.8),
        _residue(3, 7.6, overlap_sidechain=True),
    )
    chain = ProteinChain(
        "protein",
        "1",
        "A",
        residues=tuple(item.residue_id for item in residues),
        sequence="AAA",
        residue_records=residues,
        author_chain_id="A",
    )
    components = tuple(
        StructureComponent(
            f"component-{index}",
            ComponentKind.POLYMER_RESIDUE,
            item.atoms,
            residue_id=item.residue_id,
            model_id="1",
            author_chain_id="A",
            residue_name="ALA",
            auth_seq_id=str(index),
        )
        for index, item in enumerate(residues, start=1)
    ) + (
        StructureComponent(
            "ligand",
            ComponentKind.LIGAND,
            (AtomRecord("C1", "C", (0.0, 3.0, 0.0), source_atom_id="ligand-1"),),
            model_id="1",
            author_chain_id="B",
            residue_name="ATP",
            auth_seq_id="10",
        ),
    )
    selection = InputSelection(
        "a" * 64,
        "protein.pdb",
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
        diagnostics=(
            Diagnostic(
                "parser.warning",
                DiagnosticSeverity.WARNING,
                "Parser retained a recoverable source warning.",
                source_id="protein.pdb",
            ),
        ),
    )


def _empty_parsed() -> ParsedStructure:
    selection = InputSelection("c" * 64, "empty.pdb", StructureFormat.PDB, "1")
    metadata = StructureMetadata(StructureFormat.PDB, "1", ("1",), ())
    return ParsedStructure(ProteinStructure("empty"), (), selection, metadata, "d" * 64)


def test_quality_service_combines_checks_counts_and_typed_provenance() -> None:
    parsed = _parsed_with_overlap()

    report = assess_structure_quality(parsed)

    assert report.availability is Availability.AVAILABLE
    assert report.counts["chains"] == 1
    assert report.counts["polymer_residues"] == 3
    assert report.counts["components_ligand"] == 1
    assert report.counts["heavy_atom_overlaps"] >= 1
    assert "parser.warning" in {item.code for item in report.diagnostics}
    assert "quality.clashes.heavy_atom_overlap" in {item.code for item in report.diagnostics}
    assert all("ligand-1" not in item.message for item in report.diagnostics if item.code.startswith("quality.clashes"))
    assert report.provenance is not None
    assert report.provenance.input_hashes["raw_source"] == "b" * 64
    assert report.provenance.parameters["atom_scope"] == "selected_primary_polymer_heavy_atoms"
    assert report.provenance.parameters["selection_id"] == parsed.selection.selection_id
    assert report.provenance.parameters["source_connection_records"] == "not_retained_in_v0.4"
    assert report.provenance.backend_versions["overlap_screen"] == "structlens-heavy-atom-overlap-1"


def test_zero_overlaps_is_an_available_measurement_and_service_is_repeatable() -> None:
    parsed = _parsed_with_overlap()
    separated = []
    for component in parsed.components:
        if component.kind is ComponentKind.POLYMER_RESIDUE and component.auth_seq_id == "3":
            continue
        separated.append(component)
    chain = parsed.protein_structure.chains[0]
    two_residue_chain = ProteinChain(
        chain.structure_id,
        chain.model_id,
        chain.chain_id,
        residues=chain.residues[:2],
        sequence="AA",
        residue_records=chain.residue_records[:2],
        author_chain_id="A",
    )
    value = ParsedStructure(
        ProteinStructure("protein", (two_residue_chain,)),
        tuple(separated),
        parsed.selection,
        parsed.metadata,
        parsed.raw_source_hash,
    )

    first = StructureQualityService().analyze(value)
    second = StructureQualityService().analyze(value)

    assert first.availability is Availability.AVAILABLE
    assert first.counts["heavy_atom_overlaps"] == 0
    assert first == second
    assert first.provenance is not None and second.provenance is not None
    assert first.provenance.artifact_id == second.provenance.artifact_id


def test_empty_structure_quality_is_not_applicable_not_a_zero_score() -> None:
    report = assess_structure_quality(_empty_parsed())

    assert report.availability is Availability.NOT_APPLICABLE
    assert report.counts["analysis_atoms"] == 0
    assert report.counts["heavy_atom_overlaps"] == 0
    assert all("score" not in item.message.lower() for item in report.diagnostics)
