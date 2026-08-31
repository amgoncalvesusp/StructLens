from __future__ import annotations

import math

import numpy as np
import pytest

from structlens.core.models import AtomRecord
from structlens.core.quality.clashes import (
    VDW_RADII,
    VDW_RADII_VERSION,
    ClashAtom,
    ClashDiagnostic,
    ClashScreeningResult,
    screen_heavy_atom_overlaps,
)


def atom(
    atom_id: str,
    element: str,
    coordinate: tuple[float, float, float],
    *,
    atom_name: str | None = None,
    residue_id: str | None = None,
    component_id: str | None = None,
    chain_id: str | None = None,
    residue_number: int | str | None = None,
    residue_name: str | None = None,
) -> ClashAtom:
    return ClashAtom(
        atom_id,
        element,
        coordinate,
        atom_name=atom_name,
        residue_id=residue_id,
        component_id=component_id,
        chain_id=chain_id,
        residue_number=residue_number,
        residue_name=residue_name,
    )


def test_explicit_versioned_radii_and_stable_observation() -> None:
    assert VDW_RADII_VERSION
    assert VDW_RADII["C"] == pytest.approx(1.70)
    result = screen_heavy_atom_overlaps(
        [atom("z", "C", (0.0, 0.0, 0.0)), atom("a", "O", (2.0, 0.0, 0.0))],
        overlap_tolerance_angstrom=0.0,
    )

    assert isinstance(result, ClashScreeningResult)
    assert result.method_name == "heavy-atom overlap screening"
    assert result.vdw_radii_version == VDW_RADII_VERSION
    observation = result.observations[0]
    assert (observation.atom_a_id, observation.atom_b_id) == ("a", "z")
    assert (observation.element_a, observation.element_b) == ("O", "C")
    assert observation.distance_angstrom == pytest.approx(2.0)
    assert observation.overlap_angstrom == pytest.approx(1.22)
    assert (observation.vdw_radius_a_angstrom, observation.vdw_radius_b_angstrom) == pytest.approx((1.52, 1.70))


def test_cutoff_boundary_is_reported_and_overlap_is_not_reduced_by_threshold() -> None:
    boundary = screen_heavy_atom_overlaps(
        [atom("a", "C", (0.0, 0.0, 0.0)), atom("b", "C", (3.0, 0.0, 0.0))]
    )
    inside = screen_heavy_atom_overlaps(
        [atom("a", "C", (0.0, 0.0, 0.0)), atom("b", "C", (2.999999, 0.0, 0.0))]
    )

    assert boundary.observations[0].overlap_angstrom == pytest.approx(0.4)
    assert inside.observations[0].overlap_angstrom == pytest.approx(0.400001)


def test_heavy_atom_filter_same_residue_component_and_peptide_bond() -> None:
    atoms = [
        atom("same-residue-a", "C", (0.0, 0.0, 0.0), residue_id="r1", component_id="r1"),
        atom("same-residue-b", "N", (0.1, 0.0, 0.0), residue_id="r1", component_id="r1"),
        atom("same-component-a", "C", (10.0, 0.0, 0.0), residue_id="r2", component_id="cmp"),
        atom("same-component-b", "N", (10.1, 0.0, 0.0), residue_id="r3", component_id="cmp"),
        atom("peptide-c", "C", (20.0, 0.0, 0.0), atom_name="C", chain_id="A", residue_number=1),
        atom("peptide-n", "N", (21.33, 0.0, 0.0), atom_name="N", chain_id="A", residue_number=2),
        atom("adjacent-ca-a", "C", (25.0, 0.0, 0.0), atom_name="CA", chain_id="A", residue_number=1),
        atom("adjacent-ca-b", "C", (25.1, 0.0, 0.0), atom_name="CA", chain_id="A", residue_number=2),
        atom("valid-a", "C", (30.0, 0.0, 0.0), residue_id="r4", component_id="r4"),
        atom("valid-b", "O", (30.1, 0.0, 0.0), residue_id="r5", component_id="r5"),
        atom("hydrogen", "H", (30.2, 0.0, 0.0), residue_id="r5", component_id="r5"),
    ]

    result = screen_heavy_atom_overlaps(atoms, overlap_tolerance_angstrom=0.0)

    assert [(item.atom_a_id, item.atom_b_id) for item in result.observations] == [("valid-a", "valid-b")]
    assert result.hydrogen_count == 1
    assert result.screened_atom_count == 10


def test_peptide_exclusion_requires_same_chain_and_known_adjacent_numbers() -> None:
    atoms = [
        atom("c", "C", (0.0, 0.0, 0.0), atom_name="C", chain_id="A", residue_number=1),
        atom("peptide-n", "N", (1.33, 0.0, 0.0), atom_name="N", chain_id="A", residue_number=2),
        atom("n-other-chain", "N", (0.1, 0.0, 0.0), atom_name="N", chain_id="B", residue_number=2),
        atom("n-gap", "N", (10.0, 0.0, 0.0), atom_name="N", chain_id="A", residue_number=4),
        atom("n-unknown", "N", (20.0, 0.0, 0.0), atom_name="N", chain_id="A"),
    ]
    result = screen_heavy_atom_overlaps(atoms, overlap_tolerance_angstrom=0.0)

    pairs = [(item.atom_a_id, item.atom_b_id) for item in result.observations]
    assert ("c", "n-other-chain") in pairs
    assert ("c", "peptide-n") not in pairs


def test_implausibly_short_peptide_contact_is_not_hidden() -> None:
    result = screen_heavy_atom_overlaps(
        [
            atom("c", "C", (0.0, 0.0, 0.0), atom_name="C", chain_id="A", residue_number=1),
            atom("n", "N", (0.5, 0.0, 0.0), atom_name="N", chain_id="A", residue_number=2),
        ]
    )

    assert [(item.atom_a_id, item.atom_b_id) for item in result.observations] == [("c", "n")]


def test_probable_disulfide_and_explicit_source_connections_are_excluded() -> None:
    disulfide = screen_heavy_atom_overlaps(
        [
            atom("sg-a", "S", (0.0, 0.0, 0.0), atom_name="SG", residue_id="A:10", residue_name="CYS"),
            atom("sg-b", "S", (2.03, 0.0, 0.0), atom_name="SG", residue_id="A:40", residue_name="CYS"),
        ]
    )
    connected = screen_heavy_atom_overlaps(
        [atom("a", "C", (0.0, 0.0, 0.0)), atom("b", "O", (1.2, 0.0, 0.0))],
        excluded_atom_pairs=(("b", "a"),),
    )

    assert disulfide.observations == ()
    assert connected.observations == ()


def test_unknown_heavy_element_is_diagnostic_and_skipped() -> None:
    result = screen_heavy_atom_overlaps(
        [atom("known", "C", (0.0, 0.0, 0.0)), atom("mystery", "Xx", (0.1, 0.0, 0.0))]
    )

    assert result.observations == ()
    assert result.unknown_atom_count == 1
    assert result.diagnostics == (
        ClashDiagnostic(
            code="UNKNOWN_VDW_RADIUS",
            severity="warning",
            message="No validated van der Waals radius is available; atom was skipped.",
            atom_id="mystery",
            element="XX",
        ),
    )


def test_result_is_immutable_and_rigid_body_invariant() -> None:
    source = [
        atom("b", "N", (0.0, 0.0, 0.0)),
        atom("a", "O", (2.0, 0.0, 0.0)),
        atom("c", "S", (10.0, 0.0, 0.0)),
    ]
    rotation = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    transformed = [
        atom(item.atom_id, item.element, tuple(rotation @ np.asarray(item.coordinate) + (4.0, -2.0, 3.0)))
        for item in source
    ]

    first = screen_heavy_atom_overlaps(source, overlap_tolerance_angstrom=0.0)
    second = screen_heavy_atom_overlaps(transformed, overlap_tolerance_angstrom=0.0)

    assert first.observations == second.observations
    with pytest.raises((AttributeError, TypeError)):
        first.observations += ()  # type: ignore[misc]


def test_invalid_coordinates_tolerance_and_duplicate_ids_are_explicit() -> None:
    with pytest.raises(ValueError, match="finite"):
        atom("bad", "C", (math.nan, 0.0, 0.0))
    with pytest.raises(ValueError, match="tolerance"):
        screen_heavy_atom_overlaps([], overlap_tolerance_angstrom=-1.0)

    result = screen_heavy_atom_overlaps(
        [atom("same", "C", (0.0, 0.0, 0.0)), atom("same", "O", (0.1, 0.0, 0.0))]
    )
    assert result.observations == ()
    assert result.diagnostics[0].code == "DUPLICATE_ATOM_ID"


def test_large_sparse_input_uses_indexed_screen_and_returns_no_false_pairs() -> None:
    atoms = [atom(f"a-{i}", "C", (i * 10.0, 0.0, 0.0)) for i in range(2000)]

    result = screen_heavy_atom_overlaps(atoms)

    assert result.observations == ()
    assert result.screened_atom_count == 2000


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"atom_id": "", "element": "C"}, "atom_id"),
        ({"atom_id": "a", "element": ""}, "element"),
        ({"atom_id": "a", "element": "C", "coordinate": ("no", 0, 0)}, "numeric"),
        ({"atom_id": "a", "element": "C", "coordinate": (0, 0)}, "three"),
        ({"atom_id": "a", "element": "C", "residue_index": True}, "residue_index"),
        ({"atom_id": "a", "element": "C", "residue_number": True}, "residue_number"),
        ({"atom_id": "a", "element": "C", "residue_number": ""}, "residue_number"),
        ({"atom_id": "a", "element": "C", "residue_index": 1, "residue_number": "2"}, "disagree"),
    ],
)
def test_clash_atom_rejects_malformed_context(kwargs: dict[str, object], message: str) -> None:
    defaults: dict[str, object] = {"coordinate": (0.0, 0.0, 0.0)}
    defaults.update(kwargs)
    with pytest.raises(ValueError, match=message):
        ClashAtom(**defaults)  # type: ignore[arg-type]


def test_context_normalization_and_atom_record_adapter() -> None:
    value = ClashAtom(
        " a ",
        " c ",
        (0, 1, 2),
        atom_name=" c ",
        residue_id=" r ",
        component_id=" cmp ",
        chain_id=" A ",
        residue_number=" 12A ",
    )
    assert (value.atom_id, value.element, value.atom_name, value.residue_number) == ("a", "C", "C", "12A")
    record = AtomRecord("CA", "C", (2.0, 0.0, 0.0), source_atom_id="record-atom")
    result = screen_heavy_atom_overlaps([record, value], overlap_tolerance_angstrom=0.0)
    assert result.input_atom_count == 2
    assert result.observations[0].atom_a_id == "a"


def test_invalid_inputs_and_empty_known_set_are_reported_without_tree() -> None:
    with pytest.raises(TypeError, match="ClashAtom or AtomRecord"):
        screen_heavy_atom_overlaps([object()])  # type: ignore[list-item]
    with pytest.raises(ValueError, match="finite"):
        screen_heavy_atom_overlaps([], overlap_tolerance_angstrom=float("inf"))
    result = screen_heavy_atom_overlaps([atom("h", "H", (0, 0, 0))], overlap_tolerance_angstrom=10.0)
    assert result.screened_atom_count == 0
    assert result.hydrogen_count == 1
    assert not result.has_warnings


def test_local_value_objects_validate_and_keep_diagnostics_immutable() -> None:
    diagnostic = ClashDiagnostic(" code ", " WARNING ", " message ", context=(("source", 1),))
    assert (diagnostic.code, diagnostic.severity, diagnostic.context) == ("CODE", "warning", (("source", "1"),))
    with pytest.raises(ValueError, match="required"):
        ClashDiagnostic("", "warning", "message")
    with pytest.raises(ValueError, match="context"):
        ClashDiagnostic("code", "warning", "message", context=(("", "bad"),))
    with pytest.raises(ValueError, match="non-negative"):
        ClashScreeningResult(input_atom_count=-1)
    with pytest.raises(TypeError, match="observations"):
        ClashScreeningResult(observations=(object(),))  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="finite"):
        ClashScreeningResult(overlap_tolerance_angstrom=float("nan"))
