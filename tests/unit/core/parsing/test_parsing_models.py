"""Immutable contracts for parser selections and enriched structure results."""

from dataclasses import FrozenInstanceError

import pytest

from structlens.core.models import ProteinStructure
from structlens.core.parsing.models import (
    AltlocPolicy,
    AssemblyScope,
    ChainLocator,
    InputSelection,
    ParsedStructure,
    StructureFormat,
    StructureMetadata,
)

_CONTENT_ID = "a" * 64


def test_input_selection_is_deterministic_and_path_is_display_only() -> None:
    first = InputSelection(
        content_id=_CONTENT_ID,
        display_name="protein.pdb",
        format=StructureFormat.PDB,
        model_id="1",
        author_chain_ids=("A",),
        label_chain_ids=("A",),
        altloc_policy=AltlocPolicy.HIGHEST_OCCUPANCY,
        assembly_scope=AssemblyScope.ASYMMETRIC_UNIT,
        path="C:/incoming/protein.pdb",
    )
    second = InputSelection(
        content_id=_CONTENT_ID,
        display_name="protein.pdb",
        format=StructureFormat.PDB,
        model_id="1",
        author_chain_ids=("A",),
        label_chain_ids=("A",),
        altloc_policy=AltlocPolicy.HIGHEST_OCCUPANCY,
        assembly_scope=AssemblyScope.ASYMMETRIC_UNIT,
        path="D:/replacement/protein.pdb",
    )

    assert first.selection_id == second.selection_id
    assert first.identifier == first.selection_id
    with pytest.raises(FrozenInstanceError):
        first.path = "other.pdb"  # type: ignore[misc]


def test_input_selection_normalizes_collections_and_rejects_invalid_values() -> None:
    selection = InputSelection(
        content_id=_CONTENT_ID.upper(),
        display_name="protein.cif",
        format="mmcif",
        model_id=1,
        author_chain_ids=["A"],
        label_chain_ids=("A",),
        altloc_policy="all",
        assembly_scope="asymmetric_unit",
    )

    assert selection.content_id == _CONTENT_ID
    assert selection.format is StructureFormat.MMCIF
    assert selection.model_id == "1"
    assert selection.author_chain_ids == ("A",)
    assert selection.altloc_policy is AltlocPolicy.ALL
    with pytest.raises(ValueError):
        InputSelection(
            content_id="not-a-sha256",
            display_name="protein.pdb",
            format=StructureFormat.PDB,
            model_id="1",
        )


def test_metadata_and_parsed_structure_retain_typed_evidence_without_qt() -> None:
    selection = InputSelection(
        content_id=_CONTENT_ID,
        display_name="protein.pdb",
        format=StructureFormat.PDB,
        model_id="1",
        altloc_policy=AltlocPolicy.FIRST,
        assembly_scope=AssemblyScope.ASYMMETRIC_UNIT,
    )
    metadata = StructureMetadata(
        format=StructureFormat.PDB,
        selected_model="1",
        available_models=("1", "2"),
        analyzed_chain_ids=(),
        experimental_method="X-RAY DIFFRACTION",
        resolution_angstrom=1.8,
        assembly_scope=AssemblyScope.ASYMMETRIC_UNIT,
    )
    parsed = ParsedStructure(
        structure=ProteinStructure("structure"),
        components=(),
        selection=selection,
        metadata=metadata,
        raw_source_hash="b" * 64,
        diagnostics=(),
    )

    assert parsed.protein_structure.structure_id == "structure"
    assert parsed.structure is parsed.protein_structure
    assert parsed.raw_source_hash == "b" * 64
    assert parsed.source_hash == "b" * 64
    assert parsed.content_id == _CONTENT_ID
    with pytest.raises(FrozenInstanceError):
        parsed.metadata = metadata  # type: ignore[misc]


def test_parsed_structure_requires_raw_source_hash_or_explicit_raw_alias() -> None:
    selection = InputSelection(
        content_id=_CONTENT_ID,
        display_name="protein.pdb",
        format=StructureFormat.PDB,
        model_id="1",
    )
    metadata = StructureMetadata(format=StructureFormat.PDB, selected_model="1")

    with pytest.raises(TypeError, match="raw_source_hash"):
        ParsedStructure(
            structure=ProteinStructure("structure"),
            selection=selection,
            metadata=metadata,
        )
    parsed = ParsedStructure(
        structure=ProteinStructure("structure"),
        selection=selection,
        metadata=metadata,
        source_hash="b" * 64,
    )
    assert parsed.raw_source_hash == "b" * 64


def test_input_selection_requires_a_real_model_and_explicit_multi_chain_pairs() -> None:
    with pytest.raises(ValueError, match="model_id"):
        InputSelection(
            content_id=_CONTENT_ID,
            display_name="protein.pdb",
            format=StructureFormat.PDB,
            model_id=None,
        )
    with pytest.raises(ValueError, match="ambiguous"):
        InputSelection(
            content_id=_CONTENT_ID,
            display_name="protein.pdb",
            format=StructureFormat.PDB,
            model_id="1",
            author_chain_ids=("A", "B"),
            label_chain_ids=("X", "Y"),
        )

    selection = InputSelection(
        content_id=_CONTENT_ID,
        display_name="protein.pdb",
        format=StructureFormat.PDB,
        model_id="1",
        chain_locators=(
            ChainLocator("A", "Y", "2"),
            ChainLocator("B", "X", "1"),
        ),
    )
    assert selection.author_chain_ids == ("A", "B")
    assert selection.label_chain_ids == ("X", "Y")
    assert selection.chain_locators == (
        ChainLocator("A", "Y", "2"),
        ChainLocator("B", "X", "1"),
    )


def test_parsed_structure_rejects_selection_and_metadata_not_present_in_structure() -> None:
    from structlens.core.models import ProteinChain

    chain = ProteinChain(
        "structure",
        "1",
        "canonical",
        author_chain_id="AUTH",
        label_chain_id="LABEL",
        entity_id="entity-1",
    )
    structure = ProteinStructure("structure", chains=(chain,))
    selection = InputSelection(
        content_id=_CONTENT_ID,
        display_name="protein.pdb",
        format=StructureFormat.PDB,
        model_id="1",
        author_chain_ids=("AUTH",),
        label_chain_ids=("LABEL",),
    )
    metadata = StructureMetadata(
        format=StructureFormat.PDB,
        selected_model="1",
        available_models=("1",),
            analyzed_chain_ids=("canonical",),
            author_chain_ids=("AUTH",),
            label_chain_ids=("LABEL",),
            entity_ids=("entity-1",),
    )
    ParsedStructure(
        structure=structure,
        selection=selection,
        metadata=metadata,
        raw_source_hash="b" * 64,
    )

    with pytest.raises(ValueError, match="model"):
        ParsedStructure(
            structure=structure,
            selection=InputSelection(
                content_id=_CONTENT_ID,
                display_name="protein.pdb",
                format=StructureFormat.PDB,
                model_id="2",
                author_chain_ids=("AUTH",),
                label_chain_ids=("LABEL",),
            ),
            metadata=StructureMetadata(format=StructureFormat.PDB, selected_model="2"),
            raw_source_hash="b" * 64,
        )
    with pytest.raises(ValueError, match="chain"):
        ParsedStructure(
            structure=structure,
            selection=InputSelection(
                content_id=_CONTENT_ID,
                display_name="protein.pdb",
                format=StructureFormat.PDB,
                model_id="1",
                author_chain_ids=("MISSING",),
            ),
            metadata=StructureMetadata(format=StructureFormat.PDB, selected_model="1"),
            raw_source_hash="b" * 64,
        )
    with pytest.raises(ValueError, match="metadata"):
        ParsedStructure(
            structure=structure,
            selection=selection,
            metadata=StructureMetadata(
                format=StructureFormat.PDB,
                selected_model="1",
                analyzed_chain_ids=("MISSING",),
            ),
            raw_source_hash="b" * 64,
        )


def test_blank_chain_ids_and_missing_label_namespace_remain_distinct() -> None:
    from structlens.core.models import ProteinChain

    blank_chain = ProteinChain("structure", "1", "", author_chain_id="", label_chain_id="")
    structure = ProteinStructure("structure", chains=(blank_chain,))
    selection = InputSelection(
        content_id=_CONTENT_ID,
        display_name="blank-chain.pdb",
        format=StructureFormat.PDB,
        model_id="1",
        author_chain_ids=("",),
        label_chain_ids=("",),
    )
    metadata = StructureMetadata(
        format=StructureFormat.PDB,
        selected_model="1",
        analyzed_chain_ids=("",),
        author_chain_ids=("",),
        label_chain_ids=("",),
    )
    parsed = ParsedStructure(
        structure=structure,
        selection=selection,
        metadata=metadata,
        raw_source_hash="b" * 64,
    )
    assert parsed.selection.author_chain_id == ""
    assert parsed.selection.label_chain_id == ""

    absent_label_chain = ProteinChain("structure", "1", "A", author_chain_id="A")
    absent_structure = ProteinStructure("structure", chains=(absent_label_chain,))
    contradictory_selection = InputSelection(
        content_id=_CONTENT_ID,
        display_name="absent-label.pdb",
        format=StructureFormat.PDB,
        model_id="1",
        label_chain_ids=("A",),
    )
    with pytest.raises(ValueError, match="label chain"):
        ParsedStructure(
            structure=absent_structure,
            selection=contradictory_selection,
            metadata=StructureMetadata(
                format=StructureFormat.PDB,
                selected_model="1",
                analyzed_chain_ids=("A",),
                author_chain_ids=("A",),
            ),
            raw_source_hash="b" * 64,
        )

    blank_author_chain = ProteinChain("structure", "1", "A", author_chain_id="")
    blank_author_structure = ProteinStructure("structure", chains=(blank_author_chain,))
    blank_author_selection = InputSelection(
        content_id=_CONTENT_ID,
        display_name="blank-author.pdb",
        format=StructureFormat.PDB,
        model_id="1",
        author_chain_ids=("",),
    )
    ParsedStructure(
        structure=blank_author_structure,
        selection=blank_author_selection,
        metadata=StructureMetadata(
            format=StructureFormat.PDB,
            selected_model="1",
            analyzed_chain_ids=("A",),
            author_chain_ids=("",),
        ),
        raw_source_hash="b" * 64,
    )


def test_metadata_chain_sets_must_exactly_match_selected_structure_chains() -> None:
    from structlens.core.models import ProteinChain

    chain = ProteinChain(
        "structure",
        "1",
        "canonical",
        author_chain_id="AUTH",
        label_chain_id="LABEL",
        entity_id="entity-1",
    )
    structure = ProteinStructure("structure", chains=(chain,))
    selection = InputSelection(
        content_id=_CONTENT_ID,
        display_name="protein.pdb",
        format=StructureFormat.PDB,
        model_id="1",
    )
    with pytest.raises(ValueError, match="metadata"):
        ParsedStructure(
            structure=structure,
            selection=selection,
            metadata=StructureMetadata(
                format=StructureFormat.PDB,
                selected_model="1",
                analyzed_chain_ids=("AUTH",),
                author_chain_ids=("AUTH",),
                label_chain_ids=("LABEL",),
                entity_ids=("entity-1",),
            ),
            raw_source_hash="b" * 64,
        )


def test_explicit_raw_hash_and_alias_conflict_is_rejected() -> None:
    selection = InputSelection(
        content_id=_CONTENT_ID,
        display_name="protein.pdb",
        format=StructureFormat.PDB,
        model_id="1",
    )
    metadata = StructureMetadata(format=StructureFormat.PDB, selected_model="1")
    with pytest.raises(ValueError, match="aliases"):
        ParsedStructure(
            structure=ProteinStructure("structure"),
            selection=selection,
            metadata=metadata,
            raw_source_hash="b" * 64,
            source_hash="c" * 64,
        )
    with pytest.raises(ValueError, match="aliases"):
        ParsedStructure(
            structure=ProteinStructure("structure"),
            selection=selection,
            metadata=metadata,
            raw_source_hash="b" * 64,
            source_hash="b" * 64,
            raw_hash="c" * 64,
        )
