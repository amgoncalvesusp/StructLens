from __future__ import annotations

from dataclasses import FrozenInstanceError, replace

import pytest

import structlens.application.pocket_comparison_service as service_module
from structlens.application.pocket_comparison_service import PocketComparisonService
from structlens.core.evidence import (
    Availability,
    Diagnostic,
    DiagnosticSeverity,
    EvidenceCardBuilder,
    PocketEvidenceChannel,
    PocketEvidenceSections,
    build_evidence_card,
)
from structlens.core.models import ResidueCorrespondence, ResidueId
from structlens.core.pockets import AlphaSphere, FocusedPocketSelection, PocketCandidate, PocketMatchingSettings
from structlens.core.pockets.comparison import PocketSurfaceMeasurement
from structlens.core.pockets.volume_models import (
    PocketVolumeComparison,
    PocketVolumeResult,
    PocketVolumeSensitivity,
    PocketVolumeSettings,
)
from structlens.core.provenance import MethodProvenance


def _candidate(structure: str) -> PocketCandidate:
    return _lineaged_candidate(structure)


def _bare_candidate(structure: str) -> PocketCandidate:
    atom_ids = tuple(f"{structure}-atom-{index}" for index in range(4))
    residue = ResidueId(structure, "1", "A", "10", None, "ALA")
    return PocketCandidate((AlphaSphere((0.0, 0.0, 0.0), 3.0, atom_ids, (residue,), atom_ids),))


def _lineaged_candidate(structure: str, selection: str = "a" * 64, source: str = "b" * 64) -> PocketCandidate:
    atom_ids = tuple(f"{structure}-atom-{index}" for index in range(4))
    residue = ResidueId(structure, "1", "A", "10", None, "ALA")
    return PocketCandidate(
        (AlphaSphere((0.0, 0.0, 0.0), 3.0, atom_ids, (residue,), atom_ids),),
        source_content_id=source,
        selection_id=selection,
    )


def _volume(candidate: PocketCandidate, fine: float, settings: PocketVolumeSettings) -> PocketVolumeResult:
    sphere_ids = tuple(sphere.sphere_id for sphere in candidate.alpha_spheres)
    sphere_radii = tuple(sphere.radius_angstrom for sphere in candidate.alpha_spheres)
    method_parameters = {
        **dict(settings.compatibility_parameters),
        "candidate_id": candidate.candidate_id,
        "source_content_id": candidate.source_content_id,
        "selection_id": candidate.selection_id,
        "candidate_lineage": {
            "candidate_id": candidate.candidate_id,
            "source_content_id": candidate.source_content_id,
            "selection_id": candidate.selection_id,
        },
        "sphere_count": len(candidate.alpha_spheres),
        "sphere_ids": sphere_ids,
        "sphere_radii_angstrom": sphere_radii,
    }
    return PocketVolumeResult(
        Availability.AVAILABLE,
        8,
        8,
        8.0,
        fine,
        (2, 2, 2),
        (0.0, 0.0, 0.0),
        sensitivity=PocketVolumeSensitivity(0.0, 0.0),
        candidate_id=candidate.candidate_id,
        compatibility_signature=(("method", "pocket_free_volume"),) + settings.compatibility_signature,
        settings=settings,
        provenance_parameters={
            **dict(settings.compatibility_parameters),
            "candidate_id": candidate.candidate_id,
            "source_content_id": candidate.source_content_id,
            "selection_id": candidate.selection_id,
            "candidate_lineage": {
                "candidate_id": candidate.candidate_id,
                "source_content_id": candidate.source_content_id,
                "selection_id": candidate.selection_id,
            },
            "sphere_count": len(candidate.alpha_spheres),
            "sphere_ids": sphere_ids,
            "sphere_radii_angstrom": sphere_radii,
        },
        provenance=MethodProvenance(
            "structlens.pocket_free_volume",
            "0.4",
            parameters=method_parameters,
            units={
                "volume": "angstrom^3",
                "coarse_grid_spacing_angstrom": "angstrom",
                "fine_grid_spacing_angstrom": "angstrom",
                "boundary_margin_angstrom": "angstrom",
                "sphere_radii_angstrom": "angstrom",
            },
            input_hashes={
                "candidate": candidate.candidate_id,
                "source_content": candidate.source_content_id or "0" * 64,
            },
        ),
    )


def _surface(candidate: PocketCandidate, area: float) -> PocketSurfaceMeasurement:
    return PocketSurfaceMeasurement(
        candidate.candidate_id,
        candidate.source_content_id or "",
        candidate.selection_id or "",
        area,
        "tests.surface",
        {"source": "test"},
    )


def test_service_fails_closed_when_volume_settings_are_incompatible() -> None:
    reference = _candidate("ref")
    target = _candidate("target")
    first_settings = PocketVolumeSettings(coarse_grid_spacing_angstrom=1.0, fine_grid_spacing_angstrom=0.5)
    second_settings = PocketVolumeSettings(coarse_grid_spacing_angstrom=0.5, fine_grid_spacing_angstrom=0.25)

    report = PocketComparisonService().compare(
        (reference,),
        (target,),
        (
            ResidueCorrespondence(
                0,
                ResidueId("ref", "1", "A", "10", None, "ALA"),
                ResidueId("target", "1", "A", "10", None, "ALA"),
                "A",
                "A",
                "conserved",
            ),
        ),
        reference_volumes={reference.candidate_id: _volume(reference, 8.0, first_settings)},
        target_volumes={target.candidate_id: _volume(target, 12.0, second_settings)},
    )

    assert report.availability is Availability.NOT_APPLICABLE
    assert report.comparisons[0].volume_delta_angstrom3 is None
    assert any(item.code == "pocket.volume.incompatible_settings" for item in report.comparisons[0].diagnostics)


def test_builder_exposes_independent_pocket_concordance_channels() -> None:
    card = EvidenceCardBuilder(ResidueId("ref", "1", "A", "10", None, "ALA")).build(
        pocket_sections={
            "geometry": Availability.AVAILABLE,
            "ligand_support": Availability.NOT_APPLICABLE,
            "parameter_persistence": Availability.NOT_DETECTED,
            "volume_sensitivity": Availability.AVAILABLE,
            "match_ambiguity": Availability.AVAILABLE,
            "qc": Availability.NOT_APPLICABLE,
            "interactions": Availability.AVAILABLE,
        }
    )

    assert card.quality.available_sections == ("geometry", "volume_sensitivity", "match_ambiguity", "interactions")
    assert card.quality.unavailable_sections == ("ligand_support", "parameter_persistence", "qc")


def test_builder_preserves_native_channel_measure_units_diagnostics_and_provenance() -> None:
    warning = Diagnostic("pocket.warning", DiagnosticSeverity.WARNING, "not measured")
    card = EvidenceCardBuilder(ResidueId("ref", "1", "A", "10", None, "ALA")).build(
        pocket_sections={
            "geometry": PocketEvidenceChannel(
                Availability.NOT_DETECTED,
                measure={"detected": False},
                units={"count": "count"},
                diagnostics=(warning,),
                provenance={"method": "detector-v1"},
            ),
        }
    )
    assert card.pocket_sections is not None
    channel = card.pocket_sections.geometry
    assert channel.availability is Availability.NOT_DETECTED
    assert channel.measure == {"detected": False}
    assert channel.units["count"] == "count"
    assert channel.diagnostics == (warning,)
    assert channel.provenance == {"method": "detector-v1"}
    assert set(card.pocket_sections.to_mapping()) == {
        "geometry",
        "ligand_support",
        "parameter_persistence",
        "volume_sensitivity",
        "match_ambiguity",
        "qc",
        "interactions",
    }


def test_builder_serializes_native_measure_collections_without_dataclass_objects() -> None:
    measurement = PocketVolumeComparison(
        Availability.AVAILABLE,
        delta_angstrom3=1.0,
        reference_volume_angstrom3=2.0,
        target_volume_angstrom3=3.0,
    )
    card = EvidenceCardBuilder(ResidueId("ref", "1", "A", "10", None, "ALA")).build(
        pocket_sections={"volume_sensitivity": PocketEvidenceChannel(Availability.AVAILABLE, measure=(measurement,))}
    )

    serialized = card.pocket_sections.to_json()["volume_sensitivity"]  # type: ignore[union-attr]

    assert serialized["measure"][0]["delta_angstrom3"] == 1.0  # type: ignore[index]


def test_builder_validates_typed_channels_and_replaces_stale_quality_sections() -> None:
    residue = ResidueId("ref", "1", "A", "10", None, "ALA")
    with pytest.raises(ValueError):
        PocketEvidenceChannel("bad")  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        PocketEvidenceChannel(Availability.AVAILABLE, units={"": "angstrom"})
    with pytest.raises(TypeError):
        PocketEvidenceChannel(Availability.AVAILABLE, diagnostics=(object(),))  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        PocketEvidenceSections(*([Availability.AVAILABLE] * 6 + [object()]))  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        EvidenceCardBuilder(residue).build(pocket_sections=object())  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        EvidenceCardBuilder(residue).build(pocket_sections={"unknown": Availability.AVAILABLE})
    sections = PocketEvidenceSections(*([Availability.NOT_APPLICABLE] * 7))
    card = EvidenceCardBuilder(residue).build(
        quality=__import__("structlens.core.evidence", fromlist=["EvidenceQuality"]).EvidenceQuality(
            "available",
            ("geometry",),
            (),
            (),
            None,
            1,
        ),
        pocket_sections=sections,
    )
    assert card.quality.available_sections == ()
    assert card.quality.unavailable_sections == (
        "geometry",
        "ligand_support",
        "parameter_persistence",
        "volume_sensitivity",
        "match_ambiguity",
        "qc",
        "interactions",
    )
    assert build_evidence_card(residue, pocket_sections=sections).pocket_concordance is not None


def test_service_accepts_focused_selection_and_surface_mapping() -> None:
    reference = _candidate("ref")
    target = _candidate("target")
    rows = (
        ResidueCorrespondence(
            0,
            ResidueId("ref", "1", "A", "10", None, "ALA"),
            ResidueId("target", "1", "A", "10", None, "ALA"),
            "A",
            "A",
            "conserved",
        ),
    )
    service = PocketComparisonService(PocketMatchingSettings(minimum_lining_jaccard=0.0))
    report = service.analyze(
        FocusedPocketSelection(Availability.AVAILABLE, reference),
        FocusedPocketSelection(Availability.AVAILABLE, target),
        rows,
        reference_surface_areas_angstrom2={reference.candidate_id: _surface(reference, 20.0)},
        target_surface_areas_angstrom2={target.candidate_id: _surface(target, 30.0)},
    )
    assert report.status is Availability.AVAILABLE
    assert report.pocket_comparisons[0].surface_delta_angstrom2 == 10.0
    assert report.concordance is not None
    assert tuple(report.concordance.to_mapping()) == (
        "geometry",
        "ligand_support",
        "parameter_persistence",
        "volume_sensitivity",
        "match_ambiguity",
        "qc",
        "interactions",
    )
    assert report.matches == report.matching.matches
    unavailable = service.compare(
        FocusedPocketSelection(Availability.NOT_DETECTED),
        FocusedPocketSelection(Availability.NOT_DETECTED),
    )
    assert unavailable.availability is Availability.NOT_APPLICABLE
    assert unavailable.comparisons == ()


def test_service_preserves_unmatched_and_resource_failure_states() -> None:
    reference = _candidate("ref")
    target = _candidate("target")
    unmatched = PocketComparisonService().compare((reference,), (target,))
    assert unmatched.availability is Availability.NOT_APPLICABLE
    assert unmatched.comparisons[0].match.state.value == "unmatched_reference"
    with pytest.raises(ValueError, match="collision"):
        PocketComparisonService(PocketMatchingSettings(maximum_candidate_count=1)).compare(
            (reference, reference),
            (),
        )
    with pytest.raises(TypeError, match="PocketVolumeResult"):
        PocketComparisonService(PocketMatchingSettings(minimum_lining_jaccard=0.0)).compare(
            reference,
            target,
            (),
            reference_volume=PocketVolumeComparison(Availability.NUMERICAL_FAILURE),
            target_volume=PocketVolumeComparison(Availability.NUMERICAL_FAILURE),
        )


def test_service_validates_inputs_and_report_contract() -> None:
    candidate = _candidate("ref")
    with pytest.raises(TypeError):
        PocketComparisonService(object())  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        PocketComparisonService().compare((object(),), ())  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        PocketComparisonService().compare((candidate,), (), reference_volumes={candidate.candidate_id: object()})  # type: ignore[dict-item]
    with pytest.raises((TypeError, ValueError)):
        PocketComparisonService().compare((candidate,), (), reference_surface_areas_angstrom2="bad")  # type: ignore[arg-type]
    matching = PocketComparisonService().compare((), ()).matching
    diagnostic = Diagnostic("service.info", DiagnosticSeverity.INFO, "empty")
    from structlens.application.pocket_comparison_service import PocketComparisonReport

    service = PocketComparisonService()
    assert service.settings is service.settings
    with pytest.raises(TypeError, match="typed"):
        service.compare(
            candidate,
            candidate,
            (),
            settings=PocketMatchingSettings(minimum_lining_jaccard=0.0),
            reference_volume=1.0,
            target_volume=2.0,
        )
    with pytest.raises(TypeError):
        service.compare(candidate, candidate, (), settings=object())  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        service.compare(candidate, candidate, (), reference_volumes=object())  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        service.compare(candidate, candidate, (), qc_diagnostics=(object(),))  # type: ignore[arg-type]

    checked = PocketComparisonReport(Availability.NOT_APPLICABLE, matching, diagnostics=(diagnostic,))
    assert checked.to_json()["diagnostics"][0]["code"] == "service.info"
    with pytest.raises(FrozenInstanceError):
        checked.availability = Availability.AVAILABLE  # type: ignore[misc]
    for kwargs in (
        {"availability": "available", "matching": object()},
        {"availability": Availability.AVAILABLE, "matching": matching, "comparisons": (object(),)},
        {"availability": Availability.AVAILABLE, "matching": matching, "diagnostics": (object(),)},
    ):
        with pytest.raises(TypeError):
            PocketComparisonReport(**kwargs)  # type: ignore[arg-type]


def test_service_requires_complete_homogeneous_lineage_and_rejects_raw_numeric_evidence() -> None:
    reference = _lineaged_candidate("ref")
    target = _lineaged_candidate("target")
    row = ResidueCorrespondence(
        0,
        ResidueId("ref", "1", "A", "10", None, "ALA"),
        ResidueId("target", "1", "A", "10", None, "ALA"),
        "A",
        "A",
        "conserved",
    )
    with pytest.raises(ValueError, match="lineage"):
        PocketComparisonService().compare(_bare_candidate("ref"), target, (row,))
    with pytest.raises(TypeError, match="typed"):
        PocketComparisonService().compare(
            reference,
            target,
            (row,),
            reference_volume=1.0,
            target_volume=2.0,
        )
    with pytest.raises(TypeError, match="surface"):
        PocketComparisonService().compare(
            reference,
            target,
            (row,),
            reference_surface_area_angstrom2=1.0,
            target_surface_area_angstrom2=2.0,
        )
    second_selection = _lineaged_candidate("other", selection="c" * 64)
    with pytest.raises(ValueError, match="selection"):
        PocketComparisonService().compare((reference, second_selection), target, (row,))


def test_service_propagates_focused_selection_input_and_numerical_failures() -> None:
    invalid = FocusedPocketSelection(
        Availability.INVALID_INPUT,
        diagnostics=(Diagnostic("focus.invalid", DiagnosticSeverity.ERROR, "bad focus"),),
    )
    numerical = FocusedPocketSelection(
        Availability.NUMERICAL_FAILURE,
        diagnostics=(Diagnostic("focus.numeric", DiagnosticSeverity.ERROR, "numeric failure"),),
    )
    reference = _lineaged_candidate("ref")
    target = _lineaged_candidate("target")
    assert (
        PocketComparisonService().analyze(invalid, FocusedPocketSelection(Availability.AVAILABLE, target)).availability
        is Availability.INVALID_INPUT
    )
    assert (
        PocketComparisonService()
        .analyze(FocusedPocketSelection(Availability.AVAILABLE, reference), numerical)
        .availability
        is Availability.NUMERICAL_FAILURE
    )


def test_service_rejects_volume_mapping_identity_mismatch_and_duplicate_candidates() -> None:
    reference = _lineaged_candidate("ref")
    target = _lineaged_candidate("target")
    row = ResidueCorrespondence(
        0,
        ResidueId("ref", "1", "A", "10", None, "ALA"),
        ResidueId("target", "1", "A", "10", None, "ALA"),
        "A",
        "A",
        "conserved",
    )
    wrong = _volume(target, 4.0, PocketVolumeSettings())
    with pytest.raises(ValueError, match="candidate"):
        PocketComparisonService().compare(
            reference,
            target,
            (row,),
            reference_volumes={reference.candidate_id: wrong},
        )
    with pytest.raises(ValueError, match="collision"):
        PocketComparisonService().compare((reference, reference), target, (row,))


def test_service_boundary_helpers_and_native_report_states() -> None:
    reference = _lineaged_candidate("ref")
    target = _lineaged_candidate("target")
    row = ResidueCorrespondence(
        0,
        ResidueId("ref", "1", "A", "10", None, "ALA"),
        ResidueId("target", "1", "A", "10", None, "ALA"),
        "A",
        "A",
        "conserved",
    )
    diagnostic = Diagnostic("service.qc", DiagnosticSeverity.WARNING, "incomplete local QC")
    with pytest.raises(TypeError):
        service_module._candidates(iter((reference,)))  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="source"):
        PocketComparisonService().compare(
            (reference, _lineaged_candidate("other", selection="a" * 64, source="c" * 64)),
            target,
            (row,),
        )
    valid = _volume(reference, 8.0, PocketVolumeSettings())
    with pytest.raises(ValueError, match="provenance"):
        service_module._validate_volume_result(
            PocketVolumeResult(
                Availability.AVAILABLE,
                1,
                1,
                1.0,
                1.0,
                (1, 1, 1),
                (0.0, 0.0, 0.0),
                sensitivity=PocketVolumeSensitivity(0.0, 0.0),
                candidate_id=reference.candidate_id,
            ),
            reference,
        )
    with pytest.raises(ValueError, match="candidate"):
        service_module._validate_volume_result(valid, target)
    bad_provenance = MethodProvenance(
        "tests.volume",
        "1",
        parameters={
            "source_content_id": reference.source_content_id,
            "selection_id": "d" * 64,
        },
        input_hashes={"candidate": reference.candidate_id},
    )
    with pytest.raises(ValueError, match="lineage"):
        service_module._validate_volume_result(
            replace(valid, provenance=bad_provenance),
            reference,
        )
    bad_candidate_provenance = MethodProvenance(
        "tests.volume",
        "1",
        parameters={
            "source_content_id": reference.source_content_id,
            "selection_id": reference.selection_id,
        },
        input_hashes={"candidate": target.candidate_id},
    )
    with pytest.raises(ValueError, match="identify candidate"):
        service_module._validate_volume_result(
            replace(valid, provenance=bad_candidate_provenance),
            reference,
        )
    assert service_module._volume_for(valid, reference) is valid
    with pytest.raises(TypeError):
        service_module._volume_for(PocketVolumeComparison(Availability.AVAILABLE), reference)
    with pytest.raises(TypeError):
        service_module._volume_for(object(), reference)  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        service_module._volume_for({reference.candidate_id: object()}, reference)  # type: ignore[dict-item]
    with pytest.raises(TypeError):
        service_module._surface_for(1.0, reference, "surface")  # type: ignore[arg-type]
    assert service_module._surface_for({}, reference, "surface") is None
    with pytest.raises(TypeError):
        service_module._surface_for({reference.candidate_id: object()}, reference, "surface")  # type: ignore[dict-item]
    wrong_surface = PocketSurfaceMeasurement(
        reference.candidate_id, "c" * 64, reference.selection_id or "", 1.0, "m", {"x": "y"}
    )
    with pytest.raises(ValueError, match="lineage"):
        service_module._surface_for(wrong_surface, reference, "surface")
    with pytest.raises(ValueError, match="unknown"):
        service_module._validate_mapping_identity({"not-a-candidate": valid}, (reference,), "volumes")
    with pytest.raises(ValueError, match="unambiguous"):
        service_module._validate_bound_volume(valid, (reference, target), "volumes")
    with pytest.raises(TypeError):
        service_module._validate_bound_volume(
            PocketVolumeComparison(Availability.AVAILABLE), (reference, target), "volumes"
        )
    with pytest.raises(ValueError, match="unambiguous"):
        service_module._validate_bound_surface(_surface(reference, 1.0), (reference, target), "surfaces")
    with pytest.raises(ValueError, match="lineage"):
        service_module._validate_bound_surface(wrong_surface, (reference,), "surfaces")
    with pytest.raises(TypeError):
        service_module._validate_volume_input(1.0, "volumes")
    with pytest.raises(TypeError):
        service_module._validate_volume_input({reference.candidate_id: object()}, "volumes")
    with pytest.raises(TypeError):
        service_module._validate_surface_input(1.0, "surfaces")
    with pytest.raises(TypeError):
        service_module._validate_surface_input({reference.candidate_id: object()}, "surfaces")
    with pytest.raises(TypeError):
        service_module._validate_displacement_input(object())
    with pytest.raises(TypeError):
        service_module._validate_displacement_input((object(),))
    service_module._validate_qc_input((diagnostic,), "qc")

    class TypedQC:
        availability = Availability.AVAILABLE
        diagnostics = (diagnostic,)

    service_module._validate_qc_input(TypedQC(), "qc")

    class TypedQCWithoutDiagnostics:
        availability = Availability.NOT_DETECTED

    service_module._validate_qc_input(TypedQCWithoutDiagnostics(), "qc")
    with pytest.raises(TypeError):
        service_module._validate_qc_input((object(),), "qc")
    with pytest.raises(TypeError):
        service_module._validate_qc_input(object(), "qc")

    class BadAvailabilityQC:
        availability = "invalid-state"
        diagnostics = ()

    with pytest.raises(TypeError):
        service_module._validate_qc_input(BadAvailabilityQC(), "qc")

    class UnboundedDiagnosticsQC:
        availability = Availability.AVAILABLE
        diagnostics = iter((diagnostic,))

    with pytest.raises(TypeError):
        service_module._validate_qc_input(UnboundedDiagnosticsQC(), "qc")
    with pytest.raises(TypeError):
        service_module._validate_diagnostics_input(object(), "qc_diagnostics")
    with pytest.raises(TypeError):
        service_module._validate_diagnostics_input((object(),), "qc_diagnostics")
    with pytest.raises(TypeError, match="qc"):
        PocketComparisonService().compare((), (), qc=object())
    with pytest.raises(TypeError):
        service_module.PocketComparisonReport(
            Availability.AVAILABLE,
            service_module.PocketMatchingResult(Availability.NOT_APPLICABLE),
            concordance=object(),
        )  # type: ignore[arg-type]
    assert (
        PocketComparisonService()
        .compare(reference, target, (row,), reference_qc=Availability.AVAILABLE, target_qc=Availability.AVAILABLE)
        .concordance.qc.availability
        is Availability.AVAILABLE
    )  # type: ignore[union-attr]
