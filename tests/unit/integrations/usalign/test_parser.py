"""Contract tests for parsing the stable parts of US-align text output."""

from __future__ import annotations

from pathlib import Path

from structlens.integrations.usalign.parser import parse_usalign_output

USALIGN_OUTPUT = (Path(__file__).parent / "fixtures" / "usalign_basic.out").read_text()


def test_parser_extracts_primary_tm_score_version_and_metadata() -> None:
    parsed = parse_usalign_output(USALIGN_OUTPUT)

    assert parsed.tm_score == 0.9
    assert parsed.version == "20240101"
    assert parsed.metadata["structure_1"] == "reference.pdb"


def test_parser_turns_gapped_alignment_into_explicit_aligned_pairs() -> None:
    parsed = parse_usalign_output(USALIGN_OUTPUT)

    indices = [
        (pair.reference_index, pair.target_index) for pair in parsed.aligned_pairs
    ]
    assert indices == [
        (0, 0),
        (1, 1),
        (None, 2),
        (2, 3),
    ]


def test_parser_extracts_optional_rigid_transform() -> None:
    parsed = parse_usalign_output(USALIGN_OUTPUT)

    assert parsed.transform is not None
    assert parsed.transform.translation == (1.0, 2.0, 3.0)
    assert parsed.transform.rotation == (
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
    )
