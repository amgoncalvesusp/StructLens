"""Tests for the atomic coordinate quality gate."""

from __future__ import annotations

import math
from decimal import Decimal

import pytest

from structlens.core.evidence import Availability, Diagnostic, DiagnosticSeverity
from structlens.core.models import (
    AtomRecord,
    ProteinChain,
    ProteinStructure,
    ResidueId,
    ResidueNumbering,
    ResidueRecord,
)
from structlens.core.provenance import MethodProvenance
from structlens.core.quality.coordinates import (
    check_coordinate_quality,
    check_residue_quality,
    diagnostics_for_atoms,
    vdw_radius_angstrom,
)
from structlens.core.quality.models import (
    CoordinateAtom,
    CoordinateQCSettings,
    CoordinateResidue,
    StructureQualityReport,
)


def atom(name: str = "CA", element: str = "C", coordinate: object = (0.0, 0.0, 0.0), **kwargs: object) -> CoordinateAtom:
    return CoordinateAtom(name=name, element=element, coordinate=coordinate, **kwargs)


def residue(*atoms: CoordinateAtom) -> CoordinateResidue:
    residue_id = ResidueId("protein", "1", "A", "10", None, "ALA")
    return CoordinateResidue(residue_id=residue_id, atoms=atoms)


def codes(diagnostics: tuple[object, ...]) -> list[str]:
    return [diagnostic.code for diagnostic in diagnostics]  # type: ignore[union-attr]


def test_lenient_coordinate_snapshot_preserves_invalid_values_and_is_immutable() -> None:
    invalid = atom(coordinate=(math.nan, 1.0, 2.0), occupancy=1.2, b_factor=-2.0)

    assert math.isnan(float(invalid.coordinate[0]))
    assert invalid.occupancy == 1.2
    assert invalid.b_factor == -2.0
    with pytest.raises((AttributeError, TypeError)):
        invalid.coordinate[0] = 0.0  # type: ignore[index]


def test_atom_screen_reports_atomic_failures_with_context_and_stable_order() -> None:
    diagnostics = diagnostics_for_atoms(
        (
            atom("CA", "C", (math.nan, 0.0, 0.0), occupancy=1.2, b_factor=-1.0),
            atom("CB", "C", (0.0, 0.0, math.inf), occupancy=float("inf"), b_factor=float("nan")),
            atom("XX", "Q", (0.0, 0.0, 0.0)),
        ),
        source_id="input-a",
    )

    assert codes(diagnostics) == [
        "coordinate.non_finite",
        "coordinate.non_finite",
        "occupancy.invalid",
        "occupancy.invalid",
        "b_factor.invalid",
        "b_factor.invalid",
        "element.unknown_radius",
    ]
    assert diagnostics[0].severity is DiagnosticSeverity.ERROR  # type: ignore[union-attr]
    assert diagnostics[0].source_id == "input-a"  # type: ignore[union-attr]
    assert diagnostics[-1].atom_id == "XX"  # type: ignore[union-attr]


def test_residue_screen_reports_altloc_sum_duplicate_identity_and_missing_backbone() -> None:
    diagnostics = check_residue_quality(
        residue(
            atom("CA", altloc="A", occupancy=0.7),
            atom("CA", altloc="B", occupancy=0.7),
            atom("N", altloc="A", occupancy=0.7),
            atom("N", altloc="A", occupancy=0.7),
        )
    )

    assert codes(diagnostics) == [
        "altloc.occupancy_sum",
        "atom.duplicate_identity",
        "backbone.missing_atoms",
    ]
    assert diagnostics[0].residue_id.endswith(":ALA")  # type: ignore[union-attr]
    assert diagnostics[-1].atom_id is None  # type: ignore[union-attr]


def test_unknown_radius_is_explicit_and_known_radii_are_positive() -> None:
    assert vdw_radius_angstrom("C") > 0.0
    assert vdw_radius_angstrom("cl") > 0.0
    assert vdw_radius_angstrom("Q") is None


def test_quality_report_contains_typed_availability_counts_and_optional_provenance() -> None:
    report = check_coordinate_quality(
        (residue(atom("N", element="N"), atom("CA"), atom("C"), atom("O", element="O")),),
        source_id="input-a",
    )

    assert isinstance(report, StructureQualityReport)
    assert report.availability is Availability.AVAILABLE
    assert report.counts["residues"] == 1
    assert report.counts["atoms"] == 4
    assert report.provenance is None
    assert report.settings == CoordinateQCSettings()


def test_quality_report_is_unavailable_when_an_atomic_error_is_present() -> None:
    report = check_coordinate_quality((residue(atom(coordinate=(math.inf, 0.0, 0.0))),), source_id="input-a")

    assert report.availability is Availability.INVALID_INPUT
    assert "coordinate.non_finite" in codes(report.diagnostics)


def test_strict_atom_records_are_accepted_by_the_same_qc_boundary() -> None:
    strict = AtomRecord("CA", "C", (1.0, 2.0, 3.0), occupancy=0.5, b_factor=10.0)

    diagnostics = diagnostics_for_atoms((strict,))

    assert diagnostics == ()


def test_lenient_snapshot_handles_scalars_strings_and_rejects_arbitrary_objects() -> None:
    scalar = CoordinateAtom("CA", "C", Decimal("1.5"))
    text = CoordinateAtom("CA", "C", "xyz")
    assert scalar.coordinate == (1.5,)
    assert text.coordinate == ("x", "y", "z")
    with pytest.raises(TypeError):
        CoordinateAtom("CA", "C", (object(), 0.0, 0.0))
    with pytest.raises(TypeError):
        CoordinateAtom("CA", "C", (complex(1, 2), 0.0, 0.0))
    normalized = CoordinateAtom(" CA ", " c ", (0.0, 0.0, 0.0), altloc=" ", source_atom_id=" ", source_serial=" 4 ", formal_charge=1)
    assert normalized.name == "CA"
    assert normalized.element == "C"
    assert normalized.altloc is None
    assert normalized.source_atom_id is None
    assert normalized.source_serial == "4"
    assert normalized.formal_charge == 1


def test_coordinate_shape_and_non_numeric_values_are_separate_diagnostics() -> None:
    diagnostics = diagnostics_for_atoms(
        (atom("CA", coordinate=(0.0, 1.0)), atom("CB", coordinate=("bad", 1.0, 2.0)))
    )

    assert codes(diagnostics) == ["coordinate.invalid", "coordinate.invalid"]
    assert "exactly three" in diagnostics[0].message
    assert "non-numeric" in diagnostics[1].message


def test_duplicate_source_atom_id_is_screened_even_when_names_differ() -> None:
    diagnostics = diagnostics_for_atoms(
        (
            atom("CA", source_atom_id="site-1"),
            atom("CB", source_atom_id="site-1"),
        )
    )

    assert [item.code for item in diagnostics] == ["atom.duplicate_identity"]
    assert diagnostics[0].atom_id == "site-1"


def test_altloc_sum_uses_pdb_rounding_tolerance_and_missing_occupancy_is_not_zero() -> None:
    base = residue(
        atom("CA", altloc="A", occupancy=0.5),
        atom("CA", altloc="B", occupancy=0.5),
    )
    near = residue(
        atom("CA", altloc="A", occupancy=0.51),
        atom("CA", altloc="B", occupancy=0.51),
    )
    missing = residue(
        atom("CA", altloc="A"),
        atom("CA", altloc="B", occupancy=0.7),
    )

    assert "altloc.occupancy_sum" not in codes(check_residue_quality(base))
    assert "altloc.occupancy_sum" in codes(check_residue_quality(near))
    assert "altloc.occupancy_sum" not in codes(check_residue_quality(missing))


def test_structure_inputs_count_legacy_residues_and_keep_deterministic_order() -> None:
    residue_id = ResidueId("protein", "1", "A", "10", "A", "ALA")
    record = ResidueRecord(
        residue_id,
        ResidueNumbering("10", None, "A"),
        "ALA",
        "A",
        (AtomRecord("CA", "C", (0.0, 0.0, 0.0)),),
    )
    rich_chain = ProteinChain("protein", "1", "A", residue_records=(record,))
    legacy_chain = ProteinChain("protein", "1", "B", residues=(ResidueId("protein", "1", "B", "11", None, "GLY"),))
    report = check_coordinate_quality(ProteinStructure("protein", (rich_chain, legacy_chain)))

    assert report.counts["chains"] == 2
    assert report.counts["residues"] == 2
    assert "backbone.missing_atoms" in codes(report.diagnostics)
    assert any(item.residue_id is not None and item.residue_id.endswith(":GLY") for item in report.diagnostics)


def test_boundary_rejects_text_mixed_and_non_iterable_inputs() -> None:
    with pytest.raises(TypeError):
        check_coordinate_quality("CA")
    with pytest.raises(TypeError):
        check_coordinate_quality(1)
    with pytest.raises(TypeError):
        check_coordinate_quality((atom(), "not an atom"))


def test_empty_input_is_not_applicable_and_report_serializes_typed_provenance() -> None:
    provenance = MethodProvenance(
        method_id="coordinate_qc",
        method_version="1.0",
        input_hashes={"source": "a" * 64},
    )
    report = check_coordinate_quality((), provenance=provenance)

    assert report.availability is Availability.NOT_APPLICABLE
    payload = report.to_json()
    assert payload["availability"] == "not_applicable"
    assert payload["provenance"]["artifact_id"] == provenance.artifact_id  # type: ignore[index]


def test_model_contracts_validate_settings_and_snapshot_metadata() -> None:
    with pytest.raises(ValueError):
        CoordinateQCSettings(altloc_occupancy_tolerance=-1.0)
    with pytest.raises(ValueError):
        CoordinateQCSettings(occupancy_minimum=0.9, occupancy_maximum=0.1)


def test_custom_scalar_limits_are_reported_as_declared_not_as_hardcoded_defaults() -> None:
    diagnostics = diagnostics_for_atoms(
        (atom(occupancy=0.2, b_factor=2.0),),
        settings=CoordinateQCSettings(occupancy_minimum=0.3, b_factor_minimum=3.0),
    )

    assert "[0.3, 1]" in diagnostics[0].message
    assert "minimum 3" in diagnostics[1].message
    with pytest.raises(ValueError):
        StructureQualityReport("invalid")  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        StructureQualityReport(Availability.AVAILABLE, settings=object())  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        StructureQualityReport(Availability.AVAILABLE, provenance=object())  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        StructureQualityReport(Availability.AVAILABLE, counts={"atoms": -1})
    with pytest.raises(TypeError):
        CoordinateResidue("bad", ())  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        CoordinateResidue(ResidueId("protein", "1", "A", "10", None, "ALA"), (AtomRecord("CA", "C", (0, 0, 0)),))  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        CoordinateAtom.from_atom_record(object())  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        CoordinateResidue.from_residue_record(object())  # type: ignore[arg-type]


def test_report_contract_exposes_counts_and_validates_all_boundary_types() -> None:
    diagnostic = Diagnostic("test.warning", DiagnosticSeverity.WARNING, "warning")
    report = StructureQualityReport(
        "available",
        diagnostics=(diagnostic,),
        counts={"atoms": 2},
    )
    assert report.status is Availability.AVAILABLE
    assert report.error_count == 0
    assert report.warning_count == 1
    with pytest.raises(TypeError):
        StructureQualityReport(Availability.AVAILABLE, counts=[])  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        StructureQualityReport(Availability.AVAILABLE, counts={"": 1})
    with pytest.raises(ValueError):
        StructureQualityReport(Availability.AVAILABLE, counts={"atoms": True})
    with pytest.raises(TypeError):
        StructureQualityReport(Availability.AVAILABLE, diagnostics=(object(),))  # type: ignore[arg-type]
