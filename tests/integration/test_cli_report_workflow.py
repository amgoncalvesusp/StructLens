from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from structlens.application.report_serialization import deserialize_report
from structlens.cli.main import main
from structlens.core.parsing import capture_snapshot


@pytest.mark.parametrize("command", ["qc", "pockets"])
def test_cli_single_structure_commands_are_real_fixture_workflows(command: str, capsys) -> None:
    fixture = Path("tests/fixtures/parsing/numbering.mmcif")
    assert main([command, str(fixture)]) in {0, 2}
    captured = capsys.readouterr()
    assert "Status:" in captured.out + captured.err


def test_cli_compare_exports_canonical_json_tables_and_xlsx(tmp_path: Path, capsys) -> None:
    fixture = Path("tests/fixtures/parsing/numbering.mmcif")
    report_path = tmp_path / "analysis.json"
    csv_path = tmp_path / "analysis.csv"
    tsv_path = tmp_path / "analysis.tsv"
    xlsx_path = tmp_path / "analysis.xlsx"
    assert (
        main(
            [
                "compare",
                str(fixture),
                str(fixture),
                "--report",
                str(report_path),
                "--csv",
                str(csv_path),
                "--tsv",
                str(tsv_path),
                "--output",
                str(xlsx_path),
            ]
        )
        == 0
    )
    payload = json.loads(report_path.read_text(encoding="utf-8"))
    snapshots = (capture_snapshot(fixture), capture_snapshot(fixture))
    report = deserialize_report(payload, snapshots=snapshots)
    assert report.report_id == payload["report_id"]
    assert report.reference_selection.content_id == payload["reference_selection"]["content_id"]
    assert payload["provenance"]["input_hashes"]["reference_raw"] == snapshots[0].raw_sha256
    assert payload["provenance"]["input_hashes"]["reference_content"] == snapshots[0].content_id
    assert payload["provenance"]["input_hashes"]["target_raw"] == snapshots[1].raw_sha256
    assert payload["provenance"]["input_hashes"]["target_content"] == snapshots[1].content_id
    assert csv_path.read_text(encoding="utf-8").startswith("section,role,metric")
    assert tsv_path.read_text(encoding="utf-8").startswith("section\trole\tmetric")
    assert xlsx_path.stat().st_size > 0
    assert "Sequence identity:" in capsys.readouterr().out


def test_cli_compare_pdb_altloc_exports_nullable_ambiguous_msa_metrics(tmp_path: Path, capsys) -> None:
    fixture = Path("tests/fixtures/parsing/numbering_altloc.pdb")
    report_path = tmp_path / "altloc-report.json"
    assert (
        main(
            [
                "compare",
                str(fixture),
                str(fixture),
                "--mode",
                "sequence",
                "--report",
                str(report_path),
            ]
        )
        == 0
    )

    payload = json.loads(report_path.read_text(encoding="utf-8"))
    reference_snapshot = capture_snapshot(fixture)
    target_snapshot = capture_snapshot(fixture)
    report = deserialize_report(payload, snapshots=(reference_snapshot, target_snapshot))
    assert report.reference_selection.content_id == reference_snapshot.content_id
    assert report.target_selection.content_id == target_snapshot.content_id
    assert payload["provenance"]["input_hashes"] == {
        "reference_content": reference_snapshot.content_id,
        "reference_raw": reference_snapshot.raw_sha256,
        "target_content": target_snapshot.content_id,
        "target_raw": target_snapshot.raw_sha256,
    }
    ambiguous = [column for column in payload["msa"]["columns"] if column["ambiguous_fraction"] == 1.0]
    assert ambiguous
    assert all(column["conservation_score"] is None for column in ambiguous)
    assert all(column["entropy_bits"] is None for column in ambiguous)
    assert "Sequence identity:" in capsys.readouterr().out


def test_cli_help_module_smoke(capsys) -> None:
    assert main(["--help"]) == 0
    assert "usage:" in capsys.readouterr().out.lower()


def test_cli_module_help_subprocess_smoke() -> None:
    repository = Path(__file__).resolve().parents[2]
    environment = os.environ.copy()
    source_root = str(repository / "src")
    environment["PYTHONPATH"] = os.pathsep.join(filter(None, (source_root, environment.get("PYTHONPATH", ""))))
    completed = subprocess.run(
        [sys.executable, "-m", "structlens.cli.main", "--help"],
        cwd=repository,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0
    assert "pocket-compare" in completed.stdout


def test_installed_cli_help_subprocess_smoke() -> None:
    executable = shutil.which("structlens")
    if executable is None:
        pytest.skip("structlens console entrypoint is not installed")
    completed = subprocess.run([executable, "--help"], capture_output=True, text=True, check=False)
    assert completed.returncode == 0
    assert "pocket-compare" in completed.stdout


def test_cli_dispatches_pocket_compare_and_empty_command(capsys) -> None:
    fixture = Path("tests/fixtures/parsing/numbering.mmcif")
    assert main(["pocket-compare", str(fixture), str(fixture)]) in {0, 2}
    captured = capsys.readouterr()
    assert "Status:" in captured.out + captured.err or "Error:" in captured.out + captured.err
    assert main([]) == 0
    assert "usage:" in capsys.readouterr().out.lower()


def test_real_pocket_compare_writes_complete_invalid_input_evidence(tmp_path: Path, capsys) -> None:
    fixture = Path("tests/fixtures/parsing/numbering.mmcif")
    output = tmp_path / "pocket-compare.json"
    reference_snapshot = capture_snapshot(fixture)
    target_snapshot = capture_snapshot(fixture)

    assert main(["pocket-compare", str(fixture), str(fixture), "--json", str(output)]) == 2
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "Status: invalid_input" in captured.err
    assert "pocket.detect.insufficient_atoms" in captured.err

    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["reference_selection"]["content_id"] == reference_snapshot.content_id
    assert payload["target_selection"]["content_id"] == target_snapshot.content_id
    assert payload["reference_source"]["raw_sha256"] == reference_snapshot.raw_sha256
    assert payload["target_source"]["raw_sha256"] == target_snapshot.raw_sha256
    assert payload["reference_source"]["content_id"] == reference_snapshot.content_id
    assert payload["target_source"]["content_id"] == target_snapshot.content_id

    report = payload["analysis_report"]
    for key in ("input_quality", "availability", "diagnostics", "interactions", "displacement_vectors", "provenance"):
        assert key in report
    assert report["provenance"]["input_hashes"] == {
        "reference_content": reference_snapshot.content_id,
        "reference_raw": reference_snapshot.raw_sha256,
        "target_content": target_snapshot.content_id,
        "target_raw": target_snapshot.raw_sha256,
    }
    assert report["analysis"] is not None
    assert "mutations" in report["analysis"]
    assert "correspondences" in report["analysis"]
    assert "transform" in report["analysis"]

    pockets = payload["pockets"]
    assert pockets["availability"] == "invalid_input"
    assert pockets["units"] == {"length": "angstrom", "volume": "angstrom^3"}
    assert pockets["reference"]["detection"]["availability"] == "invalid_input"
    assert pockets["target"]["detection"]["availability"] == "invalid_input"
    assert pockets["reference"]["detection"]["candidates"] == []
    assert pockets["target"]["detection"]["candidates"] == []
    assert pockets["reference"]["volumes"] == []
    assert pockets["target"]["volumes"] == []
    assert pockets["comparison_report"]["matching"] == pockets["matching"]
    assert pockets["comparison_report"]["comparisons"] == pockets["comparisons"]
    assert pockets["comparison_report"]["concordance"] == pockets["concordance"]
    assert "diagnostics" in pockets["comparison_report"]
    assert pockets["diagnostics"]
    assert any(item["code"] == "pocket.detect.insufficient_atoms" for item in pockets["diagnostics"])
    assert (
        pockets["reference"]["detection"]["provenance"]["input_hashes"]["raw_source"] == reference_snapshot.raw_sha256
    )
    assert pockets["target"]["detection"]["provenance"]["input_hashes"]["raw_source"] == target_snapshot.raw_sha256
    assert (
        pockets["reference"]["detection"]["provenance"]["input_hashes"]["logical_content"]
        == reference_snapshot.content_id
    )
    assert pockets["target"]["detection"]["provenance"]["input_hashes"]["logical_content"] == target_snapshot.content_id
