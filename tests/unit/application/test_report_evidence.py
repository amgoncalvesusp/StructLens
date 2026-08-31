from __future__ import annotations

from types import SimpleNamespace

from structlens.application.report_evidence import evidence_cards, reference_ligand_atoms
from structlens.core.evidence import Availability
from structlens.core.models import (
    AnalysisResult,
    AtomRecord,
    ComponentKind,
    CorrespondenceStatus,
    ProteinChain,
    ResidueCorrespondence,
    ResidueId,
    ResidueNumbering,
    ResidueRecord,
    StructureComponent,
)
from structlens.core.quality import StructureQualityReport
from structlens.core.reports import InputQualityBundle


def _ligand(component_id: str, residue_name: str, x: float) -> StructureComponent:
    return StructureComponent(
        component_id,
        ComponentKind.LIGAND,
        (AtomRecord("C1", "C", (x, 0.0, 0.0), source_atom_id=f"{component_id}-C1"),),
        model_id="1",
        author_chain_id="L",
        residue_name=residue_name,
        auth_seq_id=component_id,
    )


def test_ligand_name_alias_is_available_only_for_one_unambiguous_component() -> None:
    first = _ligand("ligand-1", "ATP", 0.0)
    second = _ligand("ligand-2", "ATP", 20.0)

    ambiguous = reference_ligand_atoms(SimpleNamespace(components=(first, second)))  # type: ignore[arg-type]
    unique = reference_ligand_atoms(SimpleNamespace(components=(first,)))  # type: ignore[arg-type]

    assert ambiguous["ligand-1"] == first.atoms
    assert ambiguous["ligand-2"] == second.atoms
    assert "ATP" not in ambiguous
    assert unique["ATP"] == first.atoms


def test_ligand_name_alias_never_overwrites_an_exact_component_id() -> None:
    exact = _ligand("ATP", "ADP", 0.0)
    alias_candidate = _ligand("ligand-2", "ATP", 20.0)

    indexed = reference_ligand_atoms(  # type: ignore[arg-type]
        SimpleNamespace(components=(exact, alias_candidate))
    )

    assert indexed["ATP"] == exact.atoms


def test_evidence_card_marks_unmeasured_structure_section_unavailable() -> None:
    residue_id = ResidueId("reference", "1", "A", "1", None, "ALA")
    record = ResidueRecord(
        residue_id,
        ResidueNumbering("1", "1", None),
        "ALA",
        "A",
        (AtomRecord("CA", "C", (0.0, 0.0, 0.0), source_atom_id="1-CA"),),
    )
    reference = ProteinChain(
        "reference",
        "1",
        "A",
        residues=(residue_id,),
        sequence="A",
        residue_records=(record,),
    )
    deletion = ResidueCorrespondence(
        0,
        residue_id,
        None,
        "A",
        None,
        CorrespondenceStatus.DELETION,
    )
    result = AnalysisResult(
        "reference",
        "target",
        (deletion,),
        (),
        0.0,
        0.0,
        "sequence",
    )
    quality = InputQualityBundle(
        StructureQualityReport(Availability.AVAILABLE),
        StructureQualityReport(Availability.AVAILABLE),
    )

    card = evidence_cards(
        result,
        reference,
        None,
        None,
        None,
        (),
        quality,
        Availability.NOT_APPLICABLE,
        Availability.NOT_APPLICABLE,
    )[0]

    assert card.structure.available is False
    assert "structure" not in card.quality.available_sections
    assert "structure" in card.quality.unavailable_sections
