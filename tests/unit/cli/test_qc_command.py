from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from structlens.cli.commands import run_qc
from structlens.cli.inputs import capture_input
from structlens.cli.main import main


@pytest.mark.parametrize(
    ("suffix", "expected_format", "source"),
    [
        (".pdb", "pdb", "numbering_altloc.pdb"),
        (".ent", "pdb", "numbering_altloc.pdb"),
        (".cif", "mmcif", "numbering.mmcif"),
        (".mmcif", "mmcif", "numbering.mmcif"),
        (".pdb.gz", "pdb", "numbering_altloc.pdb"),
        (".ent.gz", "pdb", "numbering_altloc.pdb"),
        (".cif.gz", "mmcif", "numbering.mmcif"),
        (".mmcif.gz", "mmcif", "numbering.mmcif"),
    ],
)
def test_capture_input_supports_structure_aliases_and_gzip(
    tmp_path: Path,
    suffix: str,
    expected_format: str,
    source: str,
) -> None:
    source_bytes = (Path("tests/fixtures/parsing") / source).read_bytes()
    path = tmp_path / f"structure{suffix}"
    path.write_bytes(gzip.compress(source_bytes) if suffix.endswith(".gz") else source_bytes)
    captured = capture_input(path)
    assert captured.snapshot.logical_format == expected_format
    assert captured.selection.display_name == path.name


def test_capture_input_rejects_unsupported_format(tmp_path: Path) -> None:
    path = tmp_path / "structure.xyz"
    path.write_text("not a supported structure", encoding="utf-8")
    with pytest.raises(ValueError, match="format|supported"):
        capture_input(path)


def test_qc_renders_counts_diagnostics_and_selection() -> None:
    result = run_qc(Path("tests/fixtures/parsing/numbering.mmcif"))
    assert result.status == "available"
    assert result.payload["selection"]["model_id"] == "1"
    assert result.payload["quality"]["counts"]["chains"] == 1
    assert "diagnostics" in result.payload["quality"]


def test_qc_json_is_deterministic(tmp_path: Path) -> None:
    path = Path("tests/fixtures/parsing/numbering.mmcif")
    first = run_qc(path).json_bytes()
    second = run_qc(path).json_bytes()
    assert first == second
    destination = tmp_path / "qc.json"
    assert main(["qc", str(path), "--json", str(destination)]) == 0
    assert json.loads(destination.read_text(encoding="utf-8"))["source"]["content_id"]


def test_qc_invalid_path_is_boundary_error() -> None:
    result = run_qc(Path("missing-structure.pdb"))
    assert result.exit_code != 0
    assert result.error is not None


def test_input_boundary_rejects_empty_or_duplicate_selection_values(tmp_path: Path) -> None:
    fixture = Path("tests/fixtures/parsing/numbering.mmcif")
    for kwargs, message in (
        ({"model_id": ""}, "model"),
        ({"author_chain_ids": ("A", "A")}, "repeat"),
    ):
        try:
            capture_input(fixture, **kwargs)
        except ValueError as error:
            assert message in str(error)
        else:  # pragma: no cover - defensive assertion
            raise AssertionError("invalid selection was accepted")
    try:
        capture_input(tmp_path, model_id="1")
    except ValueError as error:
        assert "regular file" in str(error)
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("directory was accepted as a structure")
