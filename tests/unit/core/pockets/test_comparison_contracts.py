from __future__ import annotations

import json

import pytest

import structlens.core.pockets.comparison_models as comparison_models_module
import structlens.core.provenance as provenance_module
from structlens.core.evidence import Availability
from structlens.core.interactions import (
    InteractionChange,
    InteractionDifference,
    InteractionRecord,
    InteractionType,
    ReferenceInteractionKey,
)
from structlens.core.models import ResidueCorrespondence, ResidueId
from structlens.core.pockets import AlphaSphere, PocketCandidate
from structlens.core.pockets.comparison import PocketComparison, PocketSurfaceMeasurement, compare_pocket_match
from structlens.core.pockets.matching import PocketMatch, PocketMatchState
from structlens.core.pockets.volume import compare_pocket_volumes, measure_pocket_volume
from structlens.core.pockets.volume_models import (
    PocketVolumeComparison,
    PocketVolumeResult,
    PocketVolumeSettings,
)
from structlens.core.provenance import MethodProvenance


def _residue(number: int, structure: str, name: str = "ALA") -> ResidueId:
    return ResidueId(structure, "1", "A", str(number), None, name)


def _candidate(structure: str, lining: tuple[int, ...]) -> PocketCandidate:
    atom_ids = tuple(f"{structure}-atom-{offset}" for offset in range(4))
    return PocketCandidate(
        (
            AlphaSphere(
                (0.0, 0.0, 0.0), 3.0, atom_ids, tuple(_residue(number, structure) for number in lining), atom_ids
            ),
        )
    )


def test_pocket_comparison_bounds_units_before_copy_and_rejects_oversized_strings() -> None:
    candidate = _candidate("ref", (10,))
    match = PocketMatch(candidate, None, PocketMatchState.UNMATCHED_REFERENCE)

    class OversizedUnits(dict[str, str]):
        def __len__(self) -> int:
            return provenance_module.MAX_UNIT_MAP_ITEMS + 1

        def items(self) -> object:
            raise AssertionError("oversized units must be rejected before iteration")

    with pytest.raises(ValueError, match="units.*item limit"):
        PocketComparison(Availability.NOT_APPLICABLE, match, units=OversizedUnits())

    huge = "x" * (provenance_module.MAX_UNIT_STRING_LENGTH + 1)
    for units in ({huge: "count"}, {"measure": huge}):
        with pytest.raises(ValueError, match="units.*string length"):
            PocketComparison(Availability.NOT_APPLICABLE, match, units=units)


def test_comparison_contract_rejects_invalid_nested_values_and_reports_unmatched() -> None:
    candidate = _candidate("ref", (10,))
    unmatched = PocketMatch(candidate, None, PocketMatchState.UNMATCHED_REFERENCE)
    result = compare_pocket_match(unmatched)
    assert result.availability is Availability.NOT_APPLICABLE
    assert result.lining_residue_losses == ()
    for kwargs in (
        {"availability": "bad", "match": unmatched},
        {"availability": Availability.AVAILABLE, "match": object()},
        {"availability": Availability.AVAILABLE, "match": unmatched, "volume": object()},
        {"availability": Availability.AVAILABLE, "match": unmatched, "surface_delta_angstrom2": float("nan")},
        {"availability": Availability.AVAILABLE, "match": unmatched, "lining_residue_gains": (object(),)},
        {"availability": Availability.AVAILABLE, "match": unmatched, "associated_mutations": (object(),)},
        {"availability": Availability.AVAILABLE, "match": unmatched, "interaction_changes": (object(),)},
        {"availability": Availability.AVAILABLE, "match": unmatched, "local_displacements": ((object(), 1.0),)},
        {"availability": Availability.AVAILABLE, "match": unmatched, "qc_diagnostics": (object(),)},
        {"availability": Availability.AVAILABLE, "match": unmatched, "diagnostics": (object(),)},
        {"availability": Availability.AVAILABLE, "match": unmatched, "units": {"distance": ""}},
    ):
        with pytest.raises((TypeError, ValueError)):
            PocketComparison(**kwargs)  # type: ignore[arg-type]


def test_volume_signature_excludes_candidate_geometry_for_real_measurements() -> None:
    reference = _candidate("ref", (10,))
    target = _candidate("target", (10,))
    target_sphere = AlphaSphere(
        (0.0, 0.0, 0.0),
        4.0,
        tuple(f"target-wide-atom-{index}" for index in range(4)),
        (_residue(10, "target"),),
        tuple(f"target-wide-atom-{index}" for index in range(4)),
    )
    target = PocketCandidate((target_sphere,))
    settings = PocketVolumeSettings(coarse_grid_spacing_angstrom=1.0, fine_grid_spacing_angstrom=0.5)
    reference_measure = measure_pocket_volume(reference, settings=settings)
    target_measure = measure_pocket_volume(target, settings=settings)
    result = compare_pocket_volumes(reference_measure, target_measure)
    assert result.availability is Availability.AVAILABLE
    assert result.delta_angstrom3 is not None
    assert "sphere_radii_angstrom" not in dict(reference_measure.compatibility_signature)
    assert reference_measure.provenance is not None
    assert "sphere_radii_angstrom" in reference_measure.provenance.parameters
    with pytest.raises(ValueError):
        PocketVolumeResult(
            Availability.NOT_APPLICABLE,
            compatibility_signature=(("sphere_radii_angstrom", (3.0,)),),
        )


def test_surface_measurement_keeps_method_provenance_and_rejects_incompatible_methods() -> None:
    reference = _candidate("ref", (10,))
    target = _candidate("target", (10,))
    row = ResidueCorrespondence(0, _residue(10, "ref"), _residue(10, "target"), "A", "A", "conserved")
    ref_surface = PocketSurfaceMeasurement(
        "ref-candidate", "source-ref", "selection-ref", 10.0, "mesh-v1", {"run": "r"}
    )
    target_surface = PocketSurfaceMeasurement(
        "target-candidate", "source-target", "selection-target", 12.0, "mesh-v2", {"run": "t"}
    )
    match = PocketMatch(reference, target, PocketMatchState.MATCHED)
    result = compare_pocket_match(
        match,
        correspondences=(row,),
        reference_surface_area_angstrom2=ref_surface,
        target_surface_area_angstrom2=target_surface,
    )
    assert result.surface_delta_angstrom2 is None
    assert result.reference_surface_method == "mesh-v1"
    assert result.target_surface_method == "mesh-v2"
    assert result.to_json()["reference_surface_provenance"] == {"run": "r"}
    assert any(item.code == "pocket.surface.incompatible_method" for item in result.diagnostics)


def test_interaction_filter_requires_complete_residue_identity() -> None:
    reference = _candidate("ref", (10,))
    target = _candidate("target", (10,))
    row = ResidueCorrespondence(0, _residue(10, "ref"), _residue(10, "target"), "A", "A", "conserved")
    key = ReferenceInteractionKey(InteractionType.HBOND_GEOMETRIC, "A:10")
    exact = InteractionDifference(
        key,
        InteractionChange.LOST,
        reference_record=InteractionRecord(
            "ref", InteractionType.HBOND_GEOMETRIC, _residue(10, "ref"), None, "N", None, 2.0
        ),
    )
    insertion = ResidueId("ref", "1", "A", "10", "A", "ALA")
    wrong_model = ResidueId("ref", "2", "A", "10", None, "ALA")
    wrong_insertion = InteractionDifference(
        key,
        InteractionChange.LOST,
        reference_record=InteractionRecord("ref", InteractionType.HBOND_GEOMETRIC, insertion, None, "O", None, 2.0),
    )
    wrong_model_change = InteractionDifference(
        key,
        InteractionChange.LOST,
        reference_record=InteractionRecord("ref", InteractionType.HBOND_GEOMETRIC, wrong_model, None, "S", None, 2.0),
    )
    result = compare_pocket_match(
        PocketMatch(reference, target, PocketMatchState.MATCHED),
        correspondences=(row,),
        interaction_differences=(wrong_model_change, wrong_insertion, exact),
    )
    assert result.interaction_changes == (exact,)


def test_duplicate_local_displacement_residue_is_rejected() -> None:
    reference = _candidate("ref", (10,))
    target = _candidate("target", (10,))
    row = ResidueCorrespondence(0, _residue(10, "ref"), _residue(10, "target"), "A", "A", "conserved")
    vectors = (
        (_residue(10, "ref"), 1.0),
        (_residue(10, "ref"), 2.0),
    )
    with pytest.raises(ValueError, match="duplicate"):
        compare_pocket_match(
            PocketMatch(reference, target, PocketMatchState.MATCHED),
            correspondences=(row,),
            local_displacements=vectors,  # type: ignore[arg-type]
        )


def test_surface_measurement_freezes_nested_provenance_and_serializes_to_plain_json() -> None:
    provenance = {"run": {"tags": ["surface", "validated"]}}
    measurement = PocketSurfaceMeasurement("candidate", "source", "selection", 1.0, "mesh-v1", provenance)
    provenance["run"]["tags"].append("tampered")
    assert measurement.to_json()["provenance"] == {"run": {"tags": ["surface", "validated"]}}
    assert json.loads(json.dumps(measurement.to_json()))["method"] == "mesh-v1"
    with pytest.raises((TypeError, ValueError)):
        PocketSurfaceMeasurement("candidate", "source", "selection", 1.0, "mesh-v1", object())  # type: ignore[arg-type]
    with pytest.raises((TypeError, ValueError)):
        PocketSurfaceMeasurement("candidate", "source", "selection", 1.0, "mesh-v1", {"value": float("nan")})


def test_surface_measurement_rejects_provenance_payload_over_collection_limit() -> None:
    with pytest.raises(ValueError, match="limit"):
        PocketSurfaceMeasurement(
            "candidate",
            "source",
            "selection",
            1.0,
            "mesh-v1",
            {"items": tuple(range(100_001))},
        )


def test_surface_measurement_units_are_bounded_before_copy_and_bound_string_sizes() -> None:
    class OversizedUnits(dict[str, str]):
        def __len__(self) -> int:
            return comparison_models_module.MAX_SURFACE_PROVENANCE_ITEMS + 1

        def items(self) -> object:
            raise AssertionError("oversized units must not be materialized")

    with pytest.raises(ValueError, match="units.*limit"):
        PocketSurfaceMeasurement(
            "candidate",
            "source",
            "selection",
            1.0,
            "mesh-v1",
            {"run": "one"},
            OversizedUnits(),
        )
    with pytest.raises(ValueError, match="units.*string length"):
        PocketSurfaceMeasurement(
            "candidate",
            "source",
            "selection",
            1.0,
            "mesh-v1",
            {"run": "one"},
            {"x" * (comparison_models_module.MAX_SURFACE_PROVENANCE_STRING_LENGTH + 1): "angstrom^2"},
        )
    with pytest.raises(ValueError, match="units.*string length"):
        PocketSurfaceMeasurement(
            "candidate",
            "source",
            "selection",
            1.0,
            "mesh-v1",
            {"run": "one"},
            {"surface_area": "x" * (comparison_models_module.MAX_SURFACE_PROVENANCE_STRING_LENGTH + 1)},
        )


def test_surface_measurement_snapshots_to_json_provenance_once() -> None:
    class MutableProvenance:
        def __init__(self) -> None:
            self.payload: dict[str, object] = {"run": {"id": "initial"}}
            self.calls = 0

        def to_json(self) -> dict[str, object]:
            self.calls += 1
            return self.payload

    source = MutableProvenance()
    measurement = PocketSurfaceMeasurement("candidate", "source", "selection", 1.0, "mesh-v1", source)  # type: ignore[arg-type]
    source.payload["run"] = {"id": "tampered"}

    assert source.calls == 1
    assert measurement.to_json()["provenance"] == {"run": {"id": "initial"}}
    assert json.loads(json.dumps(measurement.to_json()))["provenance"] == {"run": {"id": "initial"}}


def test_comparison_model_preserves_typed_surface_provenance_and_aliases() -> None:
    candidate = _candidate("ref", (10,))
    match = PocketMatch(candidate, candidate, PocketMatchState.MATCHED)
    method = MethodProvenance("tests.surface", "1", input_hashes={"run": "a" * 64})
    volume = PocketVolumeComparison(
        Availability.AVAILABLE,
        delta_angstrom3=1.0,
        relative_delta_fraction=0.5,
        reference_volume_angstrom3=2.0,
        target_volume_angstrom3=3.0,
    )
    result = PocketComparison(
        Availability.AVAILABLE,
        match,
        volume=volume,
        reference_surface_method="mesh-v1",
        target_surface_method="mesh-v1",
        reference_surface_provenance=method,
        target_surface_provenance={"run": "target"},
        reference_surface_units={"surface_area": "angstrom^2"},
        target_surface_units={"surface_area": "angstrom^2"},
    )
    assert result.status is Availability.AVAILABLE
    assert result.volume_delta_angstrom3 == 1.0
    assert result.reference_volume_angstrom3 == 2.0
    assert result.target_volume_angstrom3 == 3.0
    assert result.relative_volume_delta_fraction == 0.5
    assert result.lining_gains == result.lining_residue_gains
    assert result.lining_losses == result.lining_residue_losses
    assert result.mutation_associations == result.associated_mutations
    assert result.interaction_differences == result.interaction_changes
    assert result.interaction_change == result.interaction_changes
    assert result.ca_displacement is None
    assert result.sidechain_displacement is None
    assert result.to_json()["reference_surface_provenance"]["method_id"] == "tests.surface"  # type: ignore[index]
