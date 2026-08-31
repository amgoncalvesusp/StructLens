from __future__ import annotations

import math
from dataclasses import FrozenInstanceError

import pytest

from structlens.core.evidence import Availability
from structlens.core.models import AtomRecord, ComponentKind, ResidueId, StructureComponent
from structlens.core.pockets import POCKET_RADII_VERSION, AlphaSphere, PocketCandidate
from structlens.core.pockets.volume import (
    PocketVolumeComparison,
    PocketVolumeResult,
    PocketVolumeSensitivity,
    PocketVolumeSettings,
    compare_pocket_volumes,
    measure_pocket_volume,
    sensitivity_from_volumes,
)


def _residue(index: int) -> ResidueId:
    return ResidueId("synthetic", "1", "A", str(index), None, "ALA")


def _sphere(
    center: tuple[float, float, float],
    radius: float,
    index: int,
) -> AlphaSphere:
    atom_ids = tuple(f"sphere-{index}-atom-{offset}" for offset in range(4))
    return AlphaSphere(
        center_xyz=center,
        radius_angstrom=radius,
        touching_atom_ids=atom_ids,
        lining_residues=(_residue(index),),
        source_simplex_atom_ids=atom_ids,
    )


def _candidate(*spheres: AlphaSphere) -> PocketCandidate:
    return PocketCandidate(alpha_spheres=tuple(spheres))


def _atom(
    coordinate: tuple[float, float, float],
    *,
    element: str = "H",
    source_atom_id: str = "atom-1",
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


def _settings(
    *,
    coarse: float = 1.0,
    fine: float = 1.0,
    margin: float = 0.0,
    policy: str = "protein_only",
    max_voxels: int = 100_000,
    chunk_size: int = 256,
) -> PocketVolumeSettings:
    return PocketVolumeSettings(
        coarse_grid_spacing_angstrom=coarse,
        fine_grid_spacing_angstrom=fine,
        boundary_margin_angstrom=margin,
        component_exclusion_policy=policy,
        max_voxel_count=max_voxels,
        voxel_chunk_size=chunk_size,
    )


def _measure(
    candidate: PocketCandidate | None,
    *,
    settings: PocketVolumeSettings | None = None,
    protein_atoms: tuple[AtomRecord, ...] = (),
    components: tuple[StructureComponent, ...] = (),
) -> PocketVolumeResult:
    return measure_pocket_volume(
        candidate,
        protein_atoms=protein_atoms,
        retained_components=components,
        settings=settings or _settings(),
    )


def test_one_sphere_uses_cell_centers_and_has_exact_known_grid_count() -> None:
    # A radius-1 sphere has a 2 x 2 x 2 bounding grid at spacing 1.0. Every
    # cell center is at +/-0.5 and is inside; vertex sampling would report 0.
    candidate = _candidate(_sphere((0.0, 0.0, 0.0), 1.0, 1))

    result = _measure(candidate)

    assert result.availability is Availability.AVAILABLE
    assert result.coarse_voxel_count == 8
    assert result.fine_voxel_count == 8
    assert result.coarse_volume_angstrom3 == pytest.approx(8.0)
    assert result.fine_volume_angstrom3 == pytest.approx(8.0)
    assert result.grid_shape == (2, 2, 2)
    assert result.grid_origin_xyz == pytest.approx((-1.0, -1.0, -1.0))
    assert result.grid_phase == "cell_center"
    assert result.units["volume"] == "angstrom^3"


def test_overlapping_spheres_count_the_union_once() -> None:
    # Each radius-1 sphere contributes two x-planes and four y/z cells. The
    # shared x-plane must be counted once: 8 + 8 - 4 = 12 cells.
    candidate = _candidate(
        _sphere((-0.5, 0.0, 0.0), 1.0, 1),
        _sphere((0.5, 0.0, 0.0), 1.0, 2),
    )

    result = _measure(candidate)

    assert result.coarse_voxel_count == 12
    assert result.coarse_volume_angstrom3 == pytest.approx(12.0)
    assert result.grid_shape == (3, 2, 2)


def test_grid_origin_phase_and_serialization_are_explicit_and_deterministic() -> None:
    candidate = _candidate(_sphere((2.0, -3.0, 4.0), 1.0, 1))

    first = _measure(candidate, settings=_settings(margin=0.25))
    second = _measure(candidate, settings=_settings(margin=0.25))

    assert first == second
    assert first.grid_phase == "cell_center"
    assert first.to_json()["grid"]["origin_xyz"] == pytest.approx((0.75, -4.25, 2.75))
    assert first.to_json()["grid"]["phase"] == "cell_center"
    assert first.to_json()["grid"]["sampling"] == "origin + (index + 0.5) * spacing"


def test_protein_vdw_exclusion_removes_occupied_cells_without_zero_coercion() -> None:
    candidate = _candidate(_sphere((0.0, 0.0, 0.0), 2.0, 1))
    protein_atom = _atom((0.0, 0.0, 0.0), element="H", source_atom_id="protein-h")

    geometric_void = _measure(candidate)
    protein_excluded = _measure(candidate, protein_atoms=(protein_atom,))

    assert geometric_void.coarse_voxel_count == 32
    assert protein_excluded.coarse_voxel_count == 24
    assert protein_excluded.coarse_volume_angstrom3 == pytest.approx(24.0)
    assert protein_excluded.coarse_volume_angstrom3 < geometric_void.coarse_volume_angstrom3
    assert protein_excluded.coarse_volume_angstrom3 != 0.0


def test_protein_only_keeps_holo_ligand_void_but_unoccupied_excludes_selected_ligand() -> None:
    candidate = _candidate(_sphere((0.0, 0.0, 0.0), 2.0, 1))
    ligand = _component("ligand-1", ComponentKind.LIGAND, _atom((0.0, 0.0, 0.0), source_atom_id="ligand-h"))
    water = _component("water-1", ComponentKind.WATER, _atom((0.0, 0.0, 0.0), source_atom_id="water-h"))

    protein_only = _measure(candidate, components=(ligand, water))
    unoccupied = _measure(
        candidate,
        settings=_settings(policy="unoccupied"),
        components=(ligand, water),
    )
    unoccupied_without_water = _measure(
        candidate,
        settings=_settings(policy="unoccupied"),
        components=(ligand,),
    )

    assert protein_only.coarse_voxel_count == 32
    assert unoccupied.coarse_voxel_count == 24
    assert unoccupied.coarse_voxel_count == unoccupied_without_water.coarse_voxel_count
    assert unoccupied.fine_volume_angstrom3 < protein_only.fine_volume_angstrom3
    assert protein_only.provenance_parameters["component_exclusion_policy"] == "protein_only"
    assert unoccupied.provenance_parameters["component_exclusion_policy"] == "unoccupied"


def test_unoccupied_ignores_polymer_and_components_without_explicit_selection() -> None:
    candidate = _candidate(_sphere((0.0, 0.0, 0.0), 2.0, 1))
    atom = _atom((0.0, 0.0, 0.0), source_atom_id="retained-h")
    polymer = _component("polymer-1", ComponentKind.POLYMER_RESIDUE, atom)
    unspecified_ligand = StructureComponent(
        "unspecified-ligand",
        ComponentKind.LIGAND,
        (atom,),
        model_id="1",
        author_chain_id="A",
        residue_name="LIG",
        auth_seq_id="2",
    )

    result = _measure(
        candidate,
        settings=_settings(policy="unoccupied"),
        components=(polymer, unspecified_ligand),
    )

    assert result.coarse_voxel_count == 32


def test_coarse_fine_sensitivity_is_absolute_and_relative_to_positive_fine_volume() -> None:
    candidate = _candidate(_sphere((0.0, 0.0, 0.0), 1.0, 1))

    result = _measure(candidate, settings=_settings(coarse=1.0, fine=0.5))

    assert result.coarse_voxel_count == 8
    assert result.fine_voxel_count == 32
    assert result.coarse_volume_angstrom3 == pytest.approx(8.0)
    assert result.fine_volume_angstrom3 == pytest.approx(32 * 0.5**3)
    assert result.absolute_sensitivity_angstrom3 == pytest.approx(4.0)
    assert result.relative_sensitivity == pytest.approx(1.0)


@pytest.mark.parametrize(
    ("coarse", "fine", "expected_absolute", "expected_relative"),
    (
        (0.0, 2.0, 2.0, 1.0),
        (2.0, 0.0, 2.0, None),
        (0.0, 0.0, 0.0, None),
    ),
)
def test_sensitivity_handles_zero_and_nonzero_combinations_without_nan(
    coarse: float,
    fine: float,
    expected_absolute: float,
    expected_relative: float | None,
) -> None:
    sensitivity = sensitivity_from_volumes(coarse, fine)

    assert sensitivity.absolute_angstrom3 == pytest.approx(expected_absolute)
    assert sensitivity.relative_fraction == expected_relative
    assert sensitivity.relative_fraction is None or math.isfinite(sensitivity.relative_fraction)


def test_numerical_zero_is_available_and_zero_zero_sensitivity_stays_unavailable_relative() -> None:
    # Both grid phases miss this tiny sphere. That is a measured zero, not a
    # missing candidate and not an invented NaN-valued ratio.
    candidate = _candidate(_sphere((0.5, 0.5, 0.5), 0.01, 1))

    result = _measure(candidate, settings=_settings(coarse=1.0, fine=0.5))

    assert result.availability is Availability.AVAILABLE
    assert result.coarse_voxel_count == 0
    assert result.fine_voxel_count == 0
    assert result.coarse_volume_angstrom3 == 0.0
    assert result.fine_volume_angstrom3 == 0.0
    assert result.absolute_sensitivity_angstrom3 == 0.0
    assert result.relative_sensitivity is None
    assert "NaN" not in str(result.to_json())


def test_missing_candidate_is_not_detected_and_never_masquerades_as_zero_volume() -> None:
    result = _measure(None)

    assert result.availability is Availability.NOT_APPLICABLE
    assert result.coarse_voxel_count is None
    assert result.fine_voxel_count is None
    assert result.coarse_volume_angstrom3 is None
    assert result.fine_volume_angstrom3 is None
    assert result.absolute_sensitivity_angstrom3 is None
    assert result.relative_sensitivity is None
    assert any(item.code == "pocket.volume.no_candidate" for item in result.diagnostics)


def test_voxel_cap_returns_typed_numerical_failure_with_remediation() -> None:
    candidate = _candidate(_sphere((0.0, 0.0, 0.0), 2.0, 1))

    result = _measure(candidate, settings=_settings(max_voxels=4))

    assert result.availability is Availability.NUMERICAL_FAILURE
    assert result.coarse_volume_angstrom3 is None
    assert result.fine_volume_angstrom3 is None
    diagnostic = next(item for item in result.diagnostics if item.code == "pocket.volume.resource_limit")
    assert diagnostic.remediation
    assert "voxel" in diagnostic.message.lower()


def test_chunked_and_single_chunk_measurements_are_identical() -> None:
    candidate = _candidate(
        _sphere((-0.5, 0.0, 0.0), 1.0, 1),
        _sphere((0.5, 0.0, 0.0), 1.0, 2),
    )

    one_by_one = _measure(candidate, settings=_settings(chunk_size=1))
    large_chunk = _measure(candidate, settings=_settings(chunk_size=10_000))

    assert one_by_one == large_chunk
    assert one_by_one.coarse_voxel_count == 12


def test_translation_preserves_counts_when_origin_translates_with_region() -> None:
    candidate = _candidate(_sphere((0.0, 0.0, 0.0), 1.0, 1))
    translated = _candidate(_sphere((13.25, -4.75, 9.5), 1.0, 1))

    first = _measure(candidate)
    second = _measure(translated)

    assert first.coarse_voxel_count == second.coarse_voxel_count == 8
    assert first.fine_voxel_count == second.fine_voxel_count == 8
    assert first.coarse_volume_angstrom3 == second.coarse_volume_angstrom3
    assert second.grid_origin_xyz == pytest.approx((12.25, -5.75, 8.5))


def test_rotation_is_reported_as_bounded_discretization_sensitivity_not_exact_invariance() -> None:
    candidate = _candidate(
        _sphere((0.0, 0.0, 0.0), 1.0, 1),
        _sphere((1.4, 0.0, 0.0), 0.8, 2),
    )
    rotated = _candidate(
        _sphere((0.0, 0.0, 0.0), 1.0, 1),
        _sphere((0.7 * math.sqrt(2), 0.7 * math.sqrt(2), 0.0), 0.8, 2),
    )

    first = _measure(candidate, settings=_settings(coarse=0.5, fine=0.5))
    second = _measure(rotated, settings=_settings(coarse=0.5, fine=0.5))

    assert first.rotation_error_bound_angstrom3 > 0.0
    assert second.rotation_error_bound_angstrom3 > 0.0
    assert abs(first.fine_volume_angstrom3 - second.fine_volume_angstrom3) <= max(
        first.rotation_error_bound_angstrom3,
        second.rotation_error_bound_angstrom3,
    )


def test_rotation_bound_is_the_declared_half_diagonal_boundary_shell_bound() -> None:
    result = _measure(
        _candidate(_sphere((0.0, 0.0, 0.0), 1.0, 1)),
        settings=_settings(coarse=0.5, fine=0.5),
    )
    half_diagonal = math.sqrt(3.0) * 0.5 / 2.0
    shell = (4.0 / 3.0) * math.pi * (
        (1.0 + half_diagonal) ** 3 - (1.0 - half_diagonal) ** 3
    )

    assert result.rotation_error_bound_angstrom3 == pytest.approx(2.0 * shell)


@pytest.mark.parametrize(
    "difference",
    ("radii", "spheres", "grid", "policy"),
)
def test_volume_comparison_disables_delta_when_method_defining_inputs_differ(difference: str) -> None:
    baseline_candidate = _candidate(_sphere((0.0, 0.0, 0.0), 1.0, 1))
    baseline_settings = _settings()
    baseline = _measure(baseline_candidate, settings=baseline_settings)

    if difference == "radii":
        other = _measure(_candidate(_sphere((0.0, 0.0, 0.0), 1.1, 1)), settings=baseline_settings)
    elif difference == "spheres":
        other = _measure(
            _candidate(
                _sphere((0.0, 0.0, 0.0), 1.0, 1),
                _sphere((3.0, 0.0, 0.0), 0.5, 2),
            ),
            settings=baseline_settings,
        )
    elif difference == "grid":
        other = _measure(baseline_candidate, settings=_settings(coarse=0.5, fine=0.5))
    else:
        other = _measure(baseline_candidate, settings=_settings(policy="unoccupied"))

    comparison = compare_pocket_volumes(baseline, other)

    assert comparison.availability is Availability.NOT_APPLICABLE
    assert comparison.delta_angstrom3 is None
    assert any(item.code == "pocket.volume.incompatible_settings" for item in comparison.diagnostics)


def test_volume_result_is_immutable_and_serializes_radii_version() -> None:
    result = _measure(_candidate(_sphere((0.0, 0.0, 0.0), 1.0, 1)))

    assert result.provenance_parameters["radii_version"] == POCKET_RADII_VERSION
    assert result.to_json()["units"]["volume"] == "angstrom^3"
    with pytest.raises(FrozenInstanceError):
        result.coarse_volume_angstrom3 = 99.0  # type: ignore[misc]


def test_volume_settings_reject_a_radii_version_not_used_by_the_calculation() -> None:
    with pytest.raises(ValueError, match="radii_version"):
        PocketVolumeSettings(radii_version="unavailable-table")


@pytest.mark.parametrize(
    ("kwargs", "message"),
    (
        ({"coarse_grid_spacing_angstrom": "bad"}, "coarse_grid_spacing_angstrom"),
        ({"fine_grid_spacing_angstrom": math.inf}, "fine_grid_spacing_angstrom"),
        ({"boundary_margin_angstrom": -0.1}, "boundary_margin_angstrom"),
        ({"max_voxel_count": 0}, "max_voxel_count"),
        ({"voxel_chunk_size": True}, "voxel_chunk_size"),
    ),
)
def test_volume_settings_validate_numeric_bounds_and_integer_controls(
    kwargs: dict[str, object],
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        PocketVolumeSettings(**kwargs)


def test_volume_settings_expose_deterministic_compatibility_metadata() -> None:
    settings = _settings(coarse=1.5, fine=0.75, margin=0.25, policy="unoccupied")

    assert settings.compatibility_parameters["coarse_grid_spacing_angstrom"] == 1.5
    assert settings.compatibility_parameters["component_exclusion_policy"] == "unoccupied"
    assert settings.compatibility_signature == tuple(settings.compatibility_parameters.items())
    assert settings.to_json()["voxel_chunk_size"] == 256


def test_volume_sensitivity_validates_relative_requirements_and_aliases() -> None:
    with pytest.raises(ValueError, match="relative_fraction"):
        PocketVolumeSensitivity(relative_fraction=0.5)

    sensitivity = PocketVolumeSensitivity(absolute_angstrom3=2.0, relative_fraction=0.25)

    assert sensitivity.absolute_sensitivity_angstrom3 == pytest.approx(2.0)
    assert sensitivity.relative_sensitivity == pytest.approx(0.25)
    assert sensitivity.to_json()["relative_fraction"] == pytest.approx(0.25)


def test_volume_result_rejects_a_grid_origin_without_three_coordinates() -> None:
    with pytest.raises(ValueError, match="grid_origin_xyz"):
        PocketVolumeResult(
            availability=Availability.AVAILABLE,
            coarse_voxel_count=1,
            fine_voxel_count=1,
            coarse_volume_angstrom3=1.0,
            fine_volume_angstrom3=1.0,
            grid_shape=(1, 1, 1),
            coarse_grid_shape=(1, 1, 1),
            fine_grid_shape=(1, 1, 1),
            grid_origin_xyz=(0.0, 0.0),
        )


def test_volume_result_rejects_non_volume_settings_objects() -> None:
    with pytest.raises(TypeError, match="settings"):
        PocketVolumeResult(
            availability=Availability.AVAILABLE,
            coarse_voxel_count=1,
            fine_voxel_count=1,
            coarse_volume_angstrom3=1.0,
            fine_volume_angstrom3=1.0,
            grid_shape=(1, 1, 1),
            coarse_grid_shape=(1, 1, 1),
            fine_grid_shape=(1, 1, 1),
            grid_origin_xyz=(0.0, 0.0, 0.0),
            settings=object(),
        )


def test_volume_result_coerces_string_availability_and_exposes_result_aliases() -> None:
    result = PocketVolumeResult(
        availability="available",
        coarse_voxel_count=1,
        fine_voxel_count=8,
        coarse_volume_angstrom3=1.0,
        fine_volume_angstrom3=8.0,
        grid_shape=(1, 1, 1),
        coarse_grid_shape=(1, 1, 1),
        fine_grid_shape=(2, 2, 2),
        grid_origin_xyz=(0.0, 0.0, 0.0),
        sensitivity=PocketVolumeSensitivity(absolute_angstrom3=7.0, relative_fraction=0.875),
        settings=PocketVolumeSettings(),
        candidate_id="candidate-1",
    )

    assert result.availability is Availability.AVAILABLE
    assert result.status is Availability.AVAILABLE
    assert result.absolute_sensitivity_angstrom3 == pytest.approx(7.0)
    assert result.relative_sensitivity == pytest.approx(0.875)
    assert result.to_json()["candidate_id"] == "candidate-1"


@pytest.mark.parametrize(
    ("kwargs", "error_type", "message"),
    (
        ({"grid_shape": (1, 1)}, ValueError, "grid shape"),
        ({"grid_phase": "vertex"}, ValueError, "grid_phase"),
        ({"sensitivity": object()}, TypeError, "sensitivity"),
        ({"provenance": object()}, TypeError, "provenance"),
        ({"candidate_id": "  "}, ValueError, "candidate_id"),
        ({"units": {"": "angstrom^3"}}, ValueError, "mapping keys"),
        ({"units": {"volume": math.inf}}, ValueError, "finite"),
        ({"provenance_parameters": {"selection": object()}}, TypeError, "unsupported JSON value type"),
    ),
)
def test_volume_result_rejects_invalid_public_contract_values(
    kwargs: dict[str, object],
    error_type: type[Exception],
    message: str,
) -> None:
    payload: dict[str, object] = {
        "availability": Availability.AVAILABLE,
        "coarse_voxel_count": 1,
        "fine_voxel_count": 1,
        "coarse_volume_angstrom3": 1.0,
        "fine_volume_angstrom3": 1.0,
        "grid_shape": (1, 1, 1),
        "coarse_grid_shape": (1, 1, 1),
        "fine_grid_shape": (1, 1, 1),
        "grid_origin_xyz": (0.0, 0.0, 0.0),
        "settings": PocketVolumeSettings(),
    }
    payload.update(kwargs)
    with pytest.raises(error_type, match=message):
        PocketVolumeResult(**payload)


def test_volume_comparison_coerces_string_availability_and_serializes_aliases() -> None:
    comparison = PocketVolumeComparison(
        availability="available",
        delta_angstrom3=-2.0,
        relative_delta_fraction=-0.25,
        reference_volume_angstrom3=8.0,
        target_volume_angstrom3=6.0,
    )

    assert comparison.availability is Availability.AVAILABLE
    assert comparison.status is Availability.AVAILABLE
    assert comparison.relative_delta == pytest.approx(-0.25)
    assert comparison.to_json()["relative_delta"] == pytest.approx(-0.25)
    assert comparison.to_json()["target_volume_angstrom3"] == pytest.approx(6.0)


def test_compatible_volume_comparison_reports_target_minus_reference() -> None:
    candidate = _candidate(_sphere((0.0, 0.0, 0.0), 2.0, 1))
    reference = _measure(candidate)
    target = _measure(
        candidate,
        protein_atoms=(_atom((0.0, 0.0, 0.0), source_atom_id="protein-h"),),
    )

    comparison = compare_pocket_volumes(reference, target)

    assert comparison.availability is Availability.AVAILABLE
    assert comparison.delta_angstrom3 == pytest.approx(
        target.fine_volume_angstrom3 - reference.fine_volume_angstrom3
    )
    assert comparison.reference_volume_angstrom3 == reference.fine_volume_angstrom3
    assert comparison.target_volume_angstrom3 == target.fine_volume_angstrom3
