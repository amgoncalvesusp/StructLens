"""Tests for safe US-align invocation and domain correspondence conversion."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from structlens.core.models import (
    AtomRecord,
    ProteinChain,
    ResidueId,
    ResidueNumbering,
    ResidueRecord,
    StructuralAlignmentSettings,
)
from structlens.integrations.usalign.adapter import USAlignAdapter
from structlens.integrations.usalign.executable import USAlignExecutionError

USALIGN_OUTPUT = """\
US-align (Version 20240101)
TM-score= 0.90000 (if normalized by length of Structure_1)
AC-D
:: :
ATGD
"""


def _chain(structure_id: str, sequence: str) -> ProteinChain:
    residues = tuple(
        ResidueId(structure_id, "1", "A", str(index), None, "ALA")
        for index, _ in enumerate(sequence, start=1)
    )
    records = tuple(
        ResidueRecord(
            residue,
            ResidueNumbering(residue.auth_seq_id, residue.auth_seq_id, None),
            "ALA",
            sequence[index - 1],
            (AtomRecord("CA", "C", (float(index), float(index % 2), 0.0), source_atom_id=f"{index}-CA"),),
        )
        for index, residue in enumerate(residues, start=1)
    )
    return ProteinChain(structure_id, "1", "A", residues, sequence, records)


def test_adapter_uses_argument_list_and_converts_pairs_to_correspondences(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    executable_path = tmp_path / "USalign"
    executable_path.touch()
    captured: dict[str, object] = {}

    def fake_run(*args: object, **kwargs: object) -> subprocess.CompletedProcess[str]:
        captured["args"] = args
        captured["kwargs"] = kwargs
        paths = tuple(Path(item) for item in args[0][1:])  # type: ignore[index]
        captured["paths"] = paths
        for path in paths:
            atom_lines = [line for line in path.read_text("ascii").splitlines() if line.startswith("ATOM")]
            assert atom_lines
            assert all(line[12:16].strip() == "CA" for line in atom_lines)
            assert all(line[21] == "A" for line in atom_lines)
        return subprocess.CompletedProcess(args[0], 0, USALIGN_OUTPUT, "")

    monkeypatch.setattr(
        "structlens.integrations.usalign.adapter.subprocess.run", fake_run
    )
    adapter = USAlignAdapter(executable=executable_path)

    result = adapter.align(
        _chain("ref", "ACD"), _chain("target", "ATGD"), StructuralAlignmentSettings()
    )

    command = captured["args"][0]  # type: ignore[index]
    assert command[0] == str(executable_path)
    assert tuple(command[1:]) == tuple(str(path) for path in captured["paths"])  # type: ignore[arg-type]
    assert captured["kwargs"] == {
        "capture_output": True,
        "check": False,
        "shell": False,
        "text": True,
        "timeout": 120.0,
    }
    assert [item.status.value for item in result.correspondences] == [
        "conserved",
        "substitution",
        "insertion",
        "conserved",
    ]
    assert result.tm_score == 0.9
    assert result.executable_version == "20240101"
    assert all(not path.exists() for path in captured["paths"])  # type: ignore[union-attr]


def test_adapter_raises_typed_error_when_usalign_fails(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    executable_path = tmp_path / "USalign"
    executable_path.touch()
    monkeypatch.setattr(
        "structlens.integrations.usalign.adapter.subprocess.run",
        lambda *_args, **_kwargs: subprocess.CompletedProcess([], 2, "", "bad input"),
    )
    adapter = USAlignAdapter(executable=executable_path)

    with pytest.raises(USAlignExecutionError, match="bad input"):
        adapter.align(
            _chain("ref", "A"), _chain("target", "A"), StructuralAlignmentSettings()
        )
