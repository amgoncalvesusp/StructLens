from __future__ import annotations

import json
from collections.abc import Sequence
from dataclasses import replace

import pytest

import structlens.application.pocket_comparison_service as service_module
import structlens.core.evidence.builder as evidence_builder_module
from structlens.application.pocket_comparison_service import PocketComparisonService
from structlens.core.difference_maps import ResidueDisplacementVector
from structlens.core.evidence import Availability, Diagnostic, DiagnosticSeverity, PocketEvidenceChannel
from structlens.core.models import ResidueCorrespondence, ResidueId
from structlens.core.pockets import AlphaSphere, FocusedPocketSelection, PocketCandidate
from structlens.core.pockets.comparison import PocketSurfaceMeasurement
from structlens.core.pockets.volume_models import (
    PocketVolumeComparison,
    PocketVolumeResult,
    PocketVolumeSensitivity,
    PocketVolumeSettings,
)
from structlens.core.provenance import MethodProvenance


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


def test_service_rejects_paired_volume_comparison_as_per_candidate_evidence() -> None:
    reference = _lineaged_candidate("ref")
    target = _lineaged_candidate("target")
    with pytest.raises(TypeError, match="PocketVolumeResult"):
        PocketComparisonService().compare(
            reference,
            target,
            (),
            reference_volume=PocketVolumeComparison(Availability.AVAILABLE),
        )


@pytest.mark.parametrize(
    "tamper", ("candidate_id", "parameters_candidate_id", "candidate_lineage", "source", "source_hash")
)
def test_service_rejects_every_volume_lineage_tampering(tamper: str) -> None:
    reference = _lineaged_candidate("ref")
    settings = PocketVolumeSettings()
    valid = _volume(reference, 8.0, settings)
    assert valid.provenance is not None
    if tamper == "candidate_id":
        tampered = replace(valid, candidate_id="f" * 64)
    elif tamper == "parameters_candidate_id":
        parameters = dict(valid.provenance.parameters)
        parameters["candidate_id"] = "f" * 64
        provenance = replace(valid.provenance, parameters=parameters)
        tampered = replace(valid, provenance=provenance)
    elif tamper == "candidate_lineage":
        parameters = dict(valid.provenance.parameters)
        parameters["candidate_lineage"] = {
            "candidate_id": reference.candidate_id,
            "source_content_id": "f" * 64,
            "selection_id": reference.selection_id,
        }
        provenance = replace(valid.provenance, parameters=parameters)
        tampered = replace(valid, provenance=provenance)
    elif tamper == "source":
        parameters = dict(valid.provenance.parameters)
        parameters["source_content_id"] = "f" * 64
        provenance = replace(valid.provenance, parameters=parameters)
        tampered = replace(valid, provenance=provenance)
    else:
        provenance = replace(
            valid.provenance, input_hashes={"candidate": reference.candidate_id, "source_content": "f" * 64}
        )
        tampered = replace(valid, provenance=provenance)
    with pytest.raises(ValueError, match="provenance|lineage|candidate"):
        service_module._validate_volume_result(tampered, reference)


@pytest.mark.parametrize(
    "tamper",
    (
        "settings",
        "signature",
        "result_parameters",
        "method_parameters",
        "method_id",
        "method_version",
        "result_units",
        "method_units",
    ),
)
def test_service_rejects_every_volume_method_coherence_tampering(tamper: str) -> None:
    reference = _lineaged_candidate("ref")
    settings = PocketVolumeSettings()
    valid = _volume(reference, 8.0, settings)
    assert valid.provenance is not None
    if tamper == "settings":
        tampered = replace(valid, settings=replace(settings, fine_grid_spacing_angstrom=0.25))
    elif tamper == "signature":
        tampered = replace(
            valid,
            compatibility_signature=(("method", "pocket_free_volume"), ("fine_grid_spacing_angstrom", 0.25)),
        )
    elif tamper == "result_parameters":
        tampered = replace(
            valid,
            provenance_parameters={
                **valid.provenance_parameters,
                "fine_grid_spacing_angstrom": 0.25,
            },
        )
    elif tamper == "method_parameters":
        tampered = replace(
            valid,
            provenance=replace(
                valid.provenance,
                parameters={**valid.provenance.parameters, "fine_grid_spacing_angstrom": 0.25},
            ),
        )
    elif tamper == "method_id":
        tampered = replace(valid, provenance=replace(valid.provenance, method_id="unknown.volume"))
    elif tamper == "method_version":
        tampered = replace(valid, provenance=replace(valid.provenance, method_version="9.9"))
    elif tamper == "result_units":
        tampered = replace(valid, units={"volume": "angstrom"})
    else:
        tampered = replace(
            valid,
            provenance=replace(valid.provenance, units={**valid.provenance.units, "volume": "angstrom"}),
        )

    with pytest.raises(ValueError, match="coherent|method|settings|provenance|unit"):
        service_module._validate_volume_result(tampered, reference)


def test_service_preserves_typed_surface_measurement_and_rejects_invalid_unilateral_surface() -> None:
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
    surface = _surface(reference, 10.0)
    report = PocketComparisonService().compare(
        reference,
        target,
        (row,),
        reference_surface_area_angstrom2=surface,
    )
    assert report.comparisons[0].reference_surface_method == "tests.surface"
    assert report.comparisons[0].to_json()["reference_surface_provenance"] == {"source": "test"}
    with pytest.raises(ValueError):
        PocketComparisonService().compare(
            reference, target, (row,), reference_surface_area_angstrom2=replace(surface, surface_area_angstrom2=-1.0)
        )


def test_service_validates_every_input_before_focused_failure_and_bounds_materialization() -> None:
    invalid = FocusedPocketSelection(
        Availability.INVALID_INPUT, diagnostics=(Diagnostic("focus.invalid", DiagnosticSeverity.ERROR, "bad"),)
    )
    with pytest.raises(TypeError, match="correspondence"):
        PocketComparisonService().compare(invalid, (), (object(),))  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="displacement"):
        PocketComparisonService().compare(invalid, (), (), local_displacements=object())  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="qc"):
        PocketComparisonService().compare(invalid, (), (), qc=object())

    class Oversized(Sequence[PocketCandidate]):
        def __len__(self) -> int:
            return service_module.MAX_INPUT_COLLECTION_SIZE + 1

        def __getitem__(self, index: int) -> PocketCandidate:
            raise AssertionError("oversized input must not be materialized")

    report = PocketComparisonService().compare(Oversized(), (), ())
    assert report.availability is Availability.NUMERICAL_FAILURE


def test_service_rejects_duplicate_displacements_and_json_payload_is_plain() -> None:
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
    vector = ResidueDisplacementVector(
        "A:10", row.reference, row.target, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 0.0, 1.0), 1.0
    )  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="duplicate"):
        PocketComparisonService().compare(reference, target, (row,), local_displacements=(vector, vector))
    channel = PocketEvidenceChannel(
        Availability.AVAILABLE, measure={"nested": [{"value": 1.0}]}, provenance={"run": {"id": "x"}}
    )
    payload = channel.to_json()
    assert json.loads(json.dumps(payload)) == payload
    source = {"nested": [{"value": 1.0}]}
    frozen = PocketEvidenceChannel(Availability.AVAILABLE, measure=source)
    source["nested"][0]["value"] = 2.0
    assert frozen.to_json()["measure"] == {"nested": [{"value": 1.0}]}
    with pytest.raises((TypeError, ValueError)):
        PocketEvidenceChannel(Availability.AVAILABLE, measure=object())
    with pytest.raises((TypeError, ValueError)):
        PocketEvidenceChannel(Availability.AVAILABLE, measure={"bad": float("nan")})


def test_evidence_channel_snapshots_to_json_objects_once_and_enforces_payload_limits() -> None:
    class MutableJSON:
        def __init__(self, payload: dict[str, object]) -> None:
            self.payload = payload
            self.calls = 0

        def to_json(self) -> dict[str, object]:
            self.calls += 1
            return self.payload

    measure = MutableJSON({"nested": {"value": 1.0}})
    provenance = MutableJSON({"run": {"id": "initial"}})
    channel = PocketEvidenceChannel(Availability.AVAILABLE, measure=measure, provenance=provenance)
    assert measure.calls == 1
    assert provenance.calls == 1
    measure.payload["nested"] = {"value": 2.0}
    provenance.payload["run"] = {"id": "tampered"}

    assert channel.to_json()["measure"] == {"nested": {"value": 1.0}}
    assert channel.to_json()["provenance"] == {"run": {"id": "initial"}}
    assert measure.calls == 1
    assert provenance.calls == 1
    assert json.loads(json.dumps(channel.to_json())) == channel.to_json()

    with pytest.raises(ValueError, match="limit"):
        PocketEvidenceChannel(Availability.AVAILABLE, measure={"items": tuple(range(100_001))})
    with pytest.raises(ValueError, match="string length"):
        PocketEvidenceChannel(Availability.AVAILABLE, measure="x" * 1_000_001)
    with pytest.raises(ValueError, match="string length"):
        PocketEvidenceChannel(Availability.AVAILABLE, measure={"x" * 1_000_001: True})


def test_evidence_channel_units_are_bounded_before_copy_and_bound_string_sizes() -> None:
    class OversizedUnits(dict[str, str]):
        def __len__(self) -> int:
            return evidence_builder_module.MAX_EVIDENCE_ITEMS + 1

        def items(self) -> object:
            raise AssertionError("oversized units must not be materialized")

    with pytest.raises(ValueError, match="units.*limit"):
        PocketEvidenceChannel(Availability.AVAILABLE, units=OversizedUnits())
    with pytest.raises(ValueError, match="units.*string length"):
        PocketEvidenceChannel(
            Availability.AVAILABLE,
            units={"x" * (evidence_builder_module.MAX_EVIDENCE_STRING_LENGTH + 1): "count"},
        )
    with pytest.raises(ValueError, match="units.*string length"):
        PocketEvidenceChannel(
            Availability.AVAILABLE,
            units={"count": "x" * (evidence_builder_module.MAX_EVIDENCE_STRING_LENGTH + 1)},
        )


def test_service_rejects_displacement_target_tampering_against_authoritative_map() -> None:
    reference = _lineaged_candidate("ref")
    target = _lineaged_candidate("target")
    reference_residue = ResidueId("ref", "1", "A", "10", None, "ALA")
    target_residue = ResidueId("target", "1", "A", "10", None, "ALA")
    row = ResidueCorrespondence(0, reference_residue, target_residue, "A", "A", "conserved")
    service_module._validate_correspondence_input({reference_residue: target_residue})
    vector = ResidueDisplacementVector(
        "A:10",
        reference_residue,
        target_residue,
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 1.0),
        (0.0, 0.0, 1.0),
        1.0,
    )
    service_module._validate_displacement_lineage((vector,), {reference_residue: target_residue})
    tampered = replace(vector, target_residue=ResidueId("target", "1", "B", "10", None, "ALA"))
    with pytest.raises(ValueError, match="authoritative correspondence"):
        PocketComparisonService().compare(reference, target, (row,), local_displacements=(tampered,))


def test_service_concordance_preserves_native_volume_sensitivity_and_surface_metadata() -> None:
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
    settings = PocketVolumeSettings()
    report = PocketComparisonService().compare(
        reference,
        target,
        (row,),
        reference_volumes={reference.candidate_id: _volume(reference, 8.0, settings)},
        target_volumes={target.candidate_id: _volume(target, 9.0, settings)},
        reference_surface_area_angstrom2=_surface(reference, 10.0),
        target_surface_area_angstrom2=_surface(target, 12.0),
    )
    assert report.concordance is not None
    volume_channel = report.concordance.volume_sensitivity
    assert volume_channel.availability is Availability.AVAILABLE
    assert volume_channel.measure is not None
    assert volume_channel.provenance is not None
    assert report.comparisons[0].surface_delta_angstrom2 == 2.0
    assert json.loads(json.dumps(report.to_json()))["concordance"]["volume_sensitivity"]["provenance"]
    geometry_json = report.concordance.geometry.to_json()
    assert geometry_json["provenance"]["surface_measurements"][0]["reference"]["method"] == "tests.surface"  # type: ignore[index]


def test_service_bounds_report_and_evidence_helpers_before_materialization() -> None:
    matching = PocketComparisonService().compare((), ()).matching
    with pytest.raises(TypeError):
        service_module.PocketComparisonReport(Availability.AVAILABLE, matching, comparisons=object())  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        service_module.PocketComparisonReport(Availability.AVAILABLE, matching, diagnostics=object())  # type: ignore[arg-type]

    class Oversized(Sequence[object]):
        def __len__(self) -> int:
            return service_module.MAX_INPUT_COLLECTION_SIZE + 1

        def __getitem__(self, index: int) -> object:
            raise AssertionError("oversized report input must not be materialized")

    with pytest.raises(ValueError):
        service_module.PocketComparisonReport(Availability.AVAILABLE, matching, comparisons=Oversized())  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        service_module.PocketComparisonReport(Availability.AVAILABLE, matching, diagnostics=Oversized())  # type: ignore[arg-type]

    reference = _lineaged_candidate("ref")
    valid = _volume(reference, 8.0, PocketVolumeSettings())
    with pytest.raises(TypeError):
        service_module._validate_volume_result(object(), reference)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="complete"):
        service_module._validate_volume_result(valid, _bare_candidate("bare"))
    assert service_module._volume_for(valid, None) is None
    assert service_module._volume_for({}, reference) is None
    service_module._validate_volume_input({"id": valid}, "volumes")
    with pytest.raises(ValueError):
        service_module._validate_mapping_identity({"id": valid}, (reference,), "volumes")


def test_service_propagates_comparison_invalid_and_numerical_states() -> None:
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
    error = Diagnostic("qc.error", DiagnosticSeverity.ERROR, "bad input")
    invalid = PocketComparisonService().compare(reference, target, (row,), qc=(error,))
    assert invalid.availability is Availability.INVALID_INPUT
    numerical = PocketComparisonService().compare(reference, target, (row,), qc=Availability.NUMERICAL_FAILURE)
    assert numerical.availability is Availability.NUMERICAL_FAILURE
