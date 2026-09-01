from __future__ import annotations

import argparse
from pathlib import Path

import pytest

import structlens.cli.main as main_module
from structlens.cli.commands import WorkflowResult, run_compare
from structlens.cli.main import main


def test_cli_help_returns_success(capsys) -> None:
    assert main(["--help"]) == 0
    output = capsys.readouterr().out
    assert "compare" in output
    assert "qc" in output
    assert "pockets" in output
    assert "pocket-compare" in output


def test_cli_compare_emits_reproducible_metrics(capsys) -> None:
    fixture = Path("tests/fixtures/parsing/numbering_altloc.pdb")

    assert main(["compare", str(fixture), str(fixture), "--mode", "sequence"]) == 0
    output = capsys.readouterr().out
    assert "Sequence identity" in output
    assert "Strict Cα RMSD" in output


def test_cli_compare_supports_canonical_report_and_tsv_outputs(tmp_path: Path, capsys) -> None:
    fixture = Path("tests/fixtures/parsing/numbering.mmcif")
    report = tmp_path / "report.json"
    table = tmp_path / "report.tsv"

    assert (
        main(
            [
                "compare",
                str(fixture),
                str(fixture),
                "--mode",
                "sequence",
                "--report",
                str(report),
                "--tsv",
                str(table),
            ]
        )
        == 0
    )
    assert report.exists()
    assert table.read_text(encoding="utf-8").startswith("section\trole\tmetric")
    assert "Availability:" in capsys.readouterr().out


def test_cli_rejects_unknown_command_and_conflicting_outputs(capsys, tmp_path: Path) -> None:
    assert main(["unknown"]) != 0
    assert "invalid choice" in capsys.readouterr().err
    fixture = Path("tests/fixtures/parsing/numbering.mmcif")
    destination = tmp_path / "same.json"
    assert main(["compare", str(fixture), str(fixture), "--json", str(destination), "--report", str(destination)]) != 0
    assert "output" in capsys.readouterr().err.lower()


def test_cli_invalid_input_returns_nonzero(capsys) -> None:
    assert main(["compare", "missing.pdb", "other.pdb"]) != 0
    assert "Error:" in capsys.readouterr().err


def test_cli_rejects_invalid_selection_and_manual_mode(capsys) -> None:
    fixture = Path("tests/fixtures/parsing/numbering.mmcif")
    assert main(["compare", str(fixture), str(fixture), "--reference-model", ""]) != 0
    assert "model" in capsys.readouterr().err.lower()
    assert main(["compare", str(fixture), str(fixture), "--mode", "manual"]) != 0
    assert "manual" in capsys.readouterr().err.lower()


def test_compare_workflow_rejects_manual_mode_and_output_collision() -> None:
    fixture = Path("tests/fixtures/parsing/numbering.mmcif")
    assert run_compare(fixture, fixture, mode="manual").exit_code == 2
    result = run_compare(fixture, fixture, mode="sequence", csv_path=fixture, tsv_path=fixture)
    assert result.exit_code == 2


def test_compare_workflow_rejects_json_report_aliases_at_function_boundary(tmp_path: Path) -> None:
    fixture = Path("tests/fixtures/parsing/numbering.mmcif")
    result = run_compare(
        fixture,
        fixture,
        mode="sequence",
        json_path=tmp_path / "json.json",
        report_path=tmp_path / "report.json",
    )
    assert result.exit_code == 2
    assert "either" in (result.error or "").lower()


def test_nonzero_workflow_summary_is_sent_to_stderr(monkeypatch: pytest.MonkeyPatch, capsys) -> None:
    monkeypatch.setattr(
        main_module,
        "_dispatch",
        lambda args: WorkflowResult("numerical_failure", summary="Status: numerical_failure", exit_code=2),
    )
    assert main(["pockets", "fixture.pdb"]) == 2
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "Status: numerical_failure" in captured.err


def test_dispatch_programming_errors_escape_main(monkeypatch: pytest.MonkeyPatch, capsys) -> None:
    def raise_error(error: Exception):
        def dispatch(args):
            raise error

        return dispatch

    for error in (TypeError("programming type fault"), KeyError("programming key fault")):
        monkeypatch.setattr(main_module, "_dispatch", raise_error(error))
        with pytest.raises(type(error), match=str(error).strip("'")):
            main(["qc", "fixture.pdb"])
        captured = capsys.readouterr()
        assert captured.out == ""
        assert captured.err == ""


@pytest.mark.parametrize(
    ("command", "argv", "runner"),
    [
        ("qc", ["qc", "fixture.pdb"], "run_qc"),
        ("pockets", ["pockets", "fixture.pdb"], "run_pockets"),
        ("pocket-compare", ["pocket-compare", "reference.pdb", "target.pdb"], "run_pocket_compare"),
    ],
)
def test_dispatches_each_non_compare_command(
    monkeypatch: pytest.MonkeyPatch,
    capsys,
    command: str,
    argv: list[str],
    runner: str,
) -> None:
    monkeypatch.setattr(main_module, runner, lambda *args, **kwargs: WorkflowResult("available", summary="ok"))
    assert main(argv) == 0
    assert capsys.readouterr().out.strip() == "ok"


def test_dispatch_rejects_unknown_namespace_command() -> None:
    with pytest.raises(ValueError, match="unknown command"):
        main_module._dispatch(argparse.Namespace(command="unsupported"))
