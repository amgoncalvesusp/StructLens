"""Tests for deterministic peptide and Cα geometry screening."""

from __future__ import annotations

import math

import pytest

from structlens.core.evidence import DiagnosticSeverity
from structlens.core.models import AtomRecord, ProteinChain, ResidueId, ResidueNumbering, ResidueRecord
from structlens.core.quality.geometry import (
    screen_ca_pseudo_geometry,
    screen_chain_geometry,
    screen_peptide_continuity,
)


def _residue(
    number: str,
    *,
    insertion_code: str | None = None,
    c: tuple[float, float, float] | None = (0.0, 0.0, 0.0),
    n: tuple[float, float, float] | None = (1.32, 0.0, 0.0),
    ca: tuple[float, float, float] | None = (0.0, 3.8, 0.0),
    chain_id: str = "A",
    model_id: str = "1",
) -> ResidueRecord:
    residue_id = ResidueId("structure", model_id, chain_id, number, insertion_code, "ALA")
    atoms = tuple(
        atom
        for atom in (
            AtomRecord("C", "C", c) if c is not None else None,
            AtomRecord("N", "N", n) if n is not None else None,
            AtomRecord("CA", "C", ca) if ca is not None else None,
        )
        if atom is not None
    )
    return ResidueRecord(
        residue_id=residue_id,
        numbering=ResidueNumbering(number, number, insertion_code),
        residue_name="ALA",
        one_letter="A",
        atoms=atoms,
    )


def _chain(*residues: ResidueRecord) -> ProteinChain:
    return ProteinChain(
        "structure",
        "1",
        "A",
        residues=tuple(residue.residue_id for residue in residues),
        sequence="A" * len(residues),
        residue_records=residues,
    )


def test_boundaries_are_inclusive_and_good_pairs_are_silent() -> None:
    first = _residue("10", c=(0.0, 0.0, 0.0), ca=(0.0, 0.0, 0.0))
    second = _residue("11", n=(2.0, 0.0, 0.0), ca=(2.5, 0.0, 0.0))

    assert screen_peptide_continuity(_chain(first, second)) == ()
    assert screen_ca_pseudo_geometry(_chain(first, second)) == ()


def test_outliers_have_stable_codes_context_and_deterministic_order() -> None:
    first = _residue("10", c=(0.0, 0.0, 0.0), ca=(0.0, 0.0, 0.0))
    second = _residue("11", n=(2.01, 0.0, 0.0), ca=(4.51, 0.0, 0.0))
    chain = _chain(first, second)

    peptide = screen_peptide_continuity(chain)
    alpha = screen_ca_pseudo_geometry(chain)

    assert [item.code for item in peptide] == ["quality.geometry.peptide_cn_gap"]
    assert peptide[0].severity is DiagnosticSeverity.WARNING
    assert peptide[0].residue_id == "structure:1:A:10:ALA"
    assert "observed-coordinate discontinuity" in peptide[0].message
    assert [item.code for item in alpha] == ["quality.geometry.ca_distance_outlier"]
    assert alpha[0].severity is DiagnosticSeverity.WARNING
    assert alpha[0].residue_id == "structure:1:A:10:ALA"


def test_missing_atoms_are_reported_without_fabricating_a_distance() -> None:
    first = _residue("10", c=None)
    second = _residue("11", ca=None, n=None)
    chain = _chain(first, second)

    peptide = screen_peptide_continuity(chain)
    alpha = screen_ca_pseudo_geometry(chain)

    assert [item.code for item in peptide] == ["quality.geometry.peptide_cn_unavailable"]
    assert [item.code for item in alpha] == ["quality.geometry.ca_distance_unavailable"]
    assert all(item.severity is DiagnosticSeverity.WARNING for item in (*peptide, *alpha))
    assert all("unavailable" in item.message.lower() for item in (*peptide, *alpha))


def test_numbering_gaps_and_reverse_order_are_not_treated_as_peptide_segments() -> None:
    first = _residue("10", c=(0.0, 0.0, 0.0), ca=(0.0, 0.0, 0.0))
    gap = _residue("12", n=(9.0, 0.0, 0.0), ca=(9.0, 0.0, 0.0))
    reverse = _residue("11", n=(9.0, 0.0, 0.0), ca=(9.0, 0.0, 0.0))

    assert screen_chain_geometry(_chain(first, gap)) == ()
    assert screen_chain_geometry(_chain(gap, reverse)) == ()


def test_insertion_codes_follow_source_order_without_crossing_a_gap() -> None:
    first = _residue("10", c=(0.0, 0.0, 0.0), ca=(0.0, 0.0, 0.0))
    inserted = _residue("10", insertion_code="A", n=(2.1, 0.0, 0.0), ca=(5.0, 0.0, 0.0))
    next_residue = _residue("11", n=(1.32, 0.0, 0.0), ca=(3.8, 0.0, 0.0))
    diagnostics = screen_chain_geometry(_chain(first, inserted, next_residue))

    assert [item.code for item in diagnostics] == [
        "quality.geometry.peptide_cn_gap",
        "quality.geometry.ca_distance_outlier",
        "quality.geometry.ca_distance_outlier",
    ]
    assert [item.residue_id for item in diagnostics] == [
        "structure:1:A:10:ALA",
        "structure:1:A:10:ALA",
        "structure:1:A:10A:ALA",
    ]


def test_multiple_chains_are_screened_independently_and_not_joined() -> None:
    chain_a = _chain(_residue("10"), _residue("11", n=(2.1, 0.0, 0.0)))
    chain_b = ProteinChain(
        "structure",
        "1",
        "B",
        residues=(_residue("20", chain_id="B").residue_id,),
        sequence="A",
        residue_records=(_residue("20", chain_id="B"),),
    )

    diagnostics = screen_peptide_continuity((chain_a, chain_b))

    assert len(diagnostics) == 1
    assert diagnostics[0].residue_id == "structure:1:A:10:ALA"


def test_distances_are_rigid_transform_invariant() -> None:
    first = _residue("10", c=(2.0, 3.0, 4.0), ca=(2.0, 3.0, 4.0))
    second = _residue("11", n=(3.32, 3.0, 4.0), ca=(5.8, 3.0, 4.0))
    original = _chain(first, second)

    angle = math.pi / 3.0

    def transform(point: tuple[float, float, float]) -> tuple[float, float, float]:
        x, y, z = point
        return (x * math.cos(angle) - y * math.sin(angle) + 11.0, x * math.sin(angle) + y * math.cos(angle) - 7.0, z + 2.0)

    transformed_records = []
    for record in (first, second):
        transformed_records.append(
            _residue(
                record.residue_id.auth_seq_id,
                c=transform(next(atom.coordinate for atom in record.atoms if atom.name == "C")),
                n=transform(next(atom.coordinate for atom in record.atoms if atom.name == "N")),
                ca=transform(next(atom.coordinate for atom in record.atoms if atom.name == "CA")),
            )
        )
    transformed = _chain(*transformed_records)

    assert screen_chain_geometry(original) == screen_chain_geometry(transformed)


def test_limits_are_explicit_and_invalid_ranges_are_rejected() -> None:
    chain = _chain(_residue("10"), _residue("11"))

    with pytest.raises(ValueError, match="finite and non-negative"):
        screen_peptide_continuity(chain, min_cn_distance_angstrom=float("nan"))
    with pytest.raises(ValueError, match="less than or equal"):
        screen_ca_pseudo_geometry(chain, min_ca_distance_angstrom=4.0, max_ca_distance_angstrom=3.0)
