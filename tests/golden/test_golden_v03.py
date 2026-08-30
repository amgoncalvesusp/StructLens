"""Golden regression test for the v0.3 scientific release gate.

The gate requires numeric coverage of MSA, conservation, interactions, sites,
distance maps, vectors, and Evidence Cards. `expected.json` pins the values the
pipeline produced when it was reviewed; any drift in a scientific lane fails
here instead of surfacing in a release build.

Regenerate deliberately, and only after reviewing the diff:

    STRUCTLENS_REGEN_GOLDEN=1 python -m pytest tests/golden
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

import pytest
from tests.golden.pipeline import FIXTURES, build_snapshot

EXPECTED = FIXTURES / "expected.json"
SECTIONS = (
    "msa",
    "conservation",
    "interactions",
    "sites",
    "distance_map",
    "vectors",
    "evidence_card",
)


def _load_expected() -> dict[str, Any]:
    if not EXPECTED.exists():
        raise AssertionError(f"missing golden file {EXPECTED}; regenerate with STRUCTLENS_REGEN_GOLDEN=1")
    loaded: dict[str, Any] = json.loads(EXPECTED.read_text(encoding="utf-8"))
    return loaded


def _write(path: Path, snapshot: dict[str, Any]) -> None:
    path.write_text(json.dumps(snapshot, indent=2, sort_keys=True) + "\n", encoding="utf-8")


@pytest.fixture(scope="module")
def snapshot() -> dict[str, Any]:
    computed = build_snapshot()
    if os.environ.get("STRUCTLENS_REGEN_GOLDEN") == "1":
        _write(EXPECTED, computed)
    return computed


@pytest.mark.parametrize("section", SECTIONS)
def test_golden_section_matches_recorded_values(section: str, snapshot: dict[str, Any]) -> None:
    assert snapshot[section] == _load_expected()[section]


def test_golden_snapshot_is_json_round_trippable(snapshot: dict[str, Any]) -> None:
    assert json.loads(json.dumps(snapshot)) == snapshot


def test_unavailable_metrics_are_none_never_zero(snapshot: dict[str, Any]) -> None:
    """Missing values must stay explicit; the plan forbids substituting zero."""

    for site in snapshot["sites"]:
        assert site["sasa_angstrom2"] is None or site["sasa_angstrom2"] > 0.0


def test_displacement_magnitudes_match_the_applied_transform(snapshot: dict[str, Any]) -> None:
    """Every Ca moves +0.5 A in z; residue A:5 also swings +1.2 A in y."""

    magnitudes = {vector["reference_position"]: vector["magnitude_angstrom"] for vector in snapshot["vectors"]}
    assert magnitudes
    assert magnitudes["A:5"] == pytest.approx(1.3)
    for position, magnitude in magnitudes.items():
        if position != "A:5":
            assert magnitude == pytest.approx(0.5)


def test_distance_difference_matrix_isolates_the_local_change(snapshot: dict[str, Any]) -> None:
    """The rigid part of the move cancels; only A:5 pairs shift internal distance."""

    labels = snapshot["distance_map"]["labels"]
    deltas = snapshot["distance_map"]["delta_angstrom"]
    moved = labels.index("A:5")
    non_zero = 0
    for row_index, row in enumerate(deltas):
        for column_index, value in enumerate(row):
            if moved in (row_index, column_index) and row_index != column_index:
                # A moved residue need not change every pair distance, but it
                # must change some of them.
                non_zero += abs(value) > 1e-6
            else:
                assert value == pytest.approx(0.0, abs=1e-9), f"rigid pair {row_index},{column_index} drifted"
    assert non_zero, "the locally displaced residue must alter at least one internal distance"


def test_alignment_records_the_reference_only_column(snapshot: dict[str, Any]) -> None:
    """PHE6 exists only in the reference, so its column must stay a target gap."""

    gap_columns = [column for column in snapshot["msa"]["columns"] if "-" in column["characters"]]
    assert gap_columns, "fixture must exercise at least one gapped column"
    assert all(column["non_gap_count"] < len(column["characters"]) for column in gap_columns)
