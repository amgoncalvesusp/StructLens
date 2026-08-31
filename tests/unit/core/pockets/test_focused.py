from __future__ import annotations

import pytest

from structlens.core.evidence import Availability, Diagnostic, DiagnosticSeverity
from structlens.core.models import ResidueId
from structlens.core.pockets import AlphaSphere, PocketCandidate
from structlens.core.pockets.focused import (
    FocusedPocketSelection,
    select_candidate_by_ligand_support,
    select_candidate_by_seed_residues,
)
from structlens.core.pockets.ligands import PocketLigandSupport
from structlens.core.sites import SiteDefinition, SiteDefinitionMode


def _residue(index: int) -> ResidueId:
    return ResidueId("protein", "1", "A", str(index), None, "ALA")


def _candidate(index: int, *lining: int) -> PocketCandidate:
    atom_ids = tuple(f"sphere-{index}-atom-{offset}" for offset in range(4))
    return PocketCandidate(
        (
            AlphaSphere(
                center_xyz=(float(index), 0.0, 0.0),
                radius_angstrom=2.0,
                touching_atom_ids=atom_ids,
                lining_residues=tuple(_residue(item) for item in lining),
                source_simplex_atom_ids=atom_ids,
            ),
        )
    )


def test_key_residue_focus_prefers_the_candidate_with_the_best_seed_overlap() -> None:
    first = _candidate(1, 10, 11)
    second = _candidate(2, 11, 12)

    selection = select_candidate_by_seed_residues((second, first), (_residue(10), _residue(11)))

    assert isinstance(selection, FocusedPocketSelection)
    assert selection.availability is Availability.AVAILABLE
    assert selection.candidate == first
    assert selection.selection_mode == "key_residues"


def test_ligand_focus_prefers_coverage_before_distance() -> None:
    first = _candidate(1, 10, 11)
    second = _candidate(2, 10, 12)

    selection = select_candidate_by_ligand_support(
        (
            (
                first,
                PocketLigandSupport(
                    component_id="ATP-1",
                    residue_name="ATP",
                    ligand_center_distance_angstrom=1.0,
                    atom_coverage_fraction=0.5,
                    lining_residue_overlap_fraction=0.5,
                    covered_atom_count=1,
                    atom_count=2,
                ),
            ),
            (
                second,
                PocketLigandSupport(
                    component_id="ATP-1",
                    residue_name="ATP",
                    ligand_center_distance_angstrom=5.0,
                    atom_coverage_fraction=1.0,
                    lining_residue_overlap_fraction=0.5,
                    covered_atom_count=2,
                    atom_count=2,
                ),
            ),
        ),
        ligand_id="ATP-1",
    )

    assert selection.availability is Availability.AVAILABLE
    assert selection.candidate == second
    assert selection.selection_mode == "ligand_radius"


def test_ligand_focus_reports_no_eligible_ligand_state() -> None:
    selection = select_candidate_by_ligand_support((), ligand_id="SO4-1")

    assert selection.availability is Availability.NOT_APPLICABLE
    assert selection.candidate is None
    assert any(item.code == "pocket.focus.no_eligible_ligand" for item in selection.diagnostics)


def test_focused_selection_normalizes_availability_diagnostics_and_mode() -> None:
    candidate = _candidate(1, 10, 11)

    selection = FocusedPocketSelection(
        availability="available",
        candidate=candidate,
        diagnostics=[
            Diagnostic(
                code="focused.selection",
                severity=DiagnosticSeverity.INFO,
                message="normalized",
            )
        ],
        selection_mode=" ligand_radius ",
    )

    assert selection.availability is Availability.AVAILABLE
    assert selection.diagnostics == (
        Diagnostic(
            code="focused.selection",
            severity=DiagnosticSeverity.INFO,
            message="normalized",
        ),
    )
    assert selection.selection_mode == "ligand_radius"


def test_focused_selection_requires_a_candidate_for_available_results() -> None:
    with pytest.raises(ValueError, match="candidate"):
        FocusedPocketSelection(availability=Availability.AVAILABLE)


def test_focused_selection_rejects_invalid_candidate_diagnostics_and_mode() -> None:
    candidate = _candidate(1, 10)

    with pytest.raises(TypeError, match="PocketCandidate"):
        FocusedPocketSelection(availability=Availability.AVAILABLE, candidate=object())  # type: ignore[arg-type]

    with pytest.raises(TypeError, match="Diagnostic"):
        FocusedPocketSelection(
            availability=Availability.NOT_DETECTED,
            candidate=candidate,
            diagnostics=(object(),),  # type: ignore[arg-type]
        )

    with pytest.raises(ValueError, match="selection_mode"):
        FocusedPocketSelection(
            availability=Availability.NOT_DETECTED,
            candidate=candidate,
            selection_mode="   ",
        )


def test_ligand_focus_rejects_untyped_candidate_support_pairs() -> None:
    with pytest.raises(TypeError, match="candidate_support"):
        select_candidate_by_ligand_support(((object(), object()),), ligand_id="ATP-1")  # type: ignore[arg-type]


def test_key_residue_focus_reports_not_detected_when_no_overlap_or_inputs_are_empty() -> None:
    first = _candidate(1, 10, 11)

    selection_without_overlap = select_candidate_by_seed_residues((first,), (_residue(99),))
    selection_without_candidates = select_candidate_by_seed_residues((), (_residue(10),), selection_mode="manual")
    selection_without_seeds = select_candidate_by_seed_residues((first,), (), selection_mode="manual")

    assert selection_without_overlap.availability is Availability.NOT_DETECTED
    assert selection_without_overlap.candidate is None
    assert selection_without_overlap.diagnostics[0].code == "pocket.focus.no_seed_overlap"
    assert selection_without_candidates.selection_mode == "manual"
    assert selection_without_candidates.diagnostics[0].code == "pocket.focus.no_seed_overlap"
    assert selection_without_seeds.selection_mode == "manual"


def test_ligand_focus_rejects_empty_identifier_and_reports_missing_supported_candidate() -> None:
    candidate = _candidate(1, 10, 11)
    support = PocketLigandSupport(
        component_id="ATP-1",
        residue_name="ATP",
        ligand_center_distance_angstrom=4.0,
        atom_coverage_fraction=0.0,
        lining_residue_overlap_fraction=0.0,
        covered_atom_count=0,
        atom_count=2,
    )

    with pytest.raises(ValueError, match="ligand_id"):
        select_candidate_by_ligand_support(((candidate, support),), ligand_id="  ")

    selection = select_candidate_by_ligand_support(((candidate, support),), ligand_id="ATP-1")

    assert selection.availability is Availability.NOT_DETECTED
    assert selection.candidate is None
    assert selection.diagnostics[0].code == "pocket.focus.no_ligand_supported_candidate"


def test_ligand_focus_retains_the_support_that_determined_selection() -> None:
    candidate = _candidate(1, 10, 11)
    support = PocketLigandSupport(
        component_id="ATP-1",
        residue_name="ATP",
        ligand_center_distance_angstrom=1.0,
        atom_coverage_fraction=1.0,
        lining_residue_overlap_fraction=0.5,
        covered_atom_count=2,
        atom_count=2,
    )

    selection = select_candidate_by_ligand_support(((candidate, support),), ligand_id="ATP-1")

    assert selection.ligand_support == support
    assert selection.to_json()["ligand_support"]["component_id"] == "ATP-1"
    assert selection.to_json()["candidate"]["candidate_id"] == candidate.candidate_id


def test_ligand_focus_rejects_support_for_a_different_ligand() -> None:
    candidate = _candidate(1, 10)
    support = PocketLigandSupport(
        component_id="ATP-1",
        residue_name="ATP",
        ligand_center_distance_angstrom=1.0,
        atom_coverage_fraction=1.0,
        lining_residue_overlap_fraction=0.0,
        covered_atom_count=1,
        atom_count=1,
    )

    with pytest.raises(ValueError, match="ligand_id"):
        select_candidate_by_ligand_support(((candidate, support),), ligand_id="ADP-1")


def test_focused_selection_rejects_candidate_on_unavailable_state() -> None:
    with pytest.raises(ValueError, match="unavailable"):
        FocusedPocketSelection(
            availability=Availability.NOT_DETECTED,
            candidate=_candidate(1, 10),
        )


def test_key_residue_focus_validates_residue_types() -> None:
    with pytest.raises(TypeError, match="seed_residues"):
        select_candidate_by_seed_residues(
            (_candidate(1, 10),),
            (object(),),  # type: ignore[arg-type]
        )


def test_ligand_site_definition_rejects_blank_identifiers() -> None:
    with pytest.raises(ValueError, match="ligand_id"):
        SiteDefinition(
            "focus",
            "Focus",
            SiteDefinitionMode.LIGAND_RADIUS,
            ligand_id="   ",
            radius_angstrom=4.0,
        )


def test_key_residue_focus_rejects_untyped_candidates() -> None:
    with pytest.raises(TypeError, match="PocketCandidate"):
        select_candidate_by_seed_residues((object(),), (_residue(10),))  # type: ignore[arg-type]
