from __future__ import annotations

from pathlib import Path

from structlens.cli.main import main


def test_cli_help_returns_success(capsys) -> None:
    assert main(["--help"]) == 0
    assert "compare" in capsys.readouterr().out


def test_cli_compare_emits_reproducible_metrics(capsys) -> None:
    fixture = Path("tests/fixtures/parsing/numbering_altloc.pdb")

    assert main(["compare", str(fixture), str(fixture), "--mode", "sequence"]) == 0
    output = capsys.readouterr().out
    assert "Sequence identity" in output
    assert "Strict Cα RMSD" in output


def test_cli_invalid_input_returns_nonzero(capsys) -> None:
    assert main(["compare", "missing.pdb", "other.pdb"]) != 0
    assert "Error:" in capsys.readouterr().err
