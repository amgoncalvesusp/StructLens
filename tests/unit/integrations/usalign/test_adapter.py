"""Tests for safe US-align invocation and domain correspondence conversion."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from structlens.core.models import ProteinChain, ResidueId, StructuralAlignmentSettings
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
    return ProteinChain(structure_id, "1", "A", residues, sequence)


def test_adapter_uses_argument_list_and_converts_pairs_to_correspondences(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    reference_path = tmp_path / "reference.pdb"
    target_path = tmp_path / "target.pdb"
    executable_path = tmp_path / "USalign"
    reference_path.touch()
    target_path.touch()
    executable_path.touch()
    captured: dict[str, object] = {}

    def fake_run(*args: object, **kwargs: object) -> subprocess.CompletedProcess[str]:
        captured["args"] = args
        captured["kwargs"] = kwargs
        return subprocess.CompletedProcess(args[0], 0, USALIGN_OUTPUT, "")

    monkeypatch.setattr(
        "structlens.integrations.usalign.adapter.subprocess.run", fake_run
    )
    adapter = USAlignAdapter(
        executable=executable_path,
        structure_paths={"ref": reference_path, "target": target_path},
    )

    result = adapter.align(
        _chain("ref", "ACD"), _chain("target", "ATGD"), StructuralAlignmentSettings()
    )

    assert captured["args"] == (
        [str(executable_path), str(reference_path), str(target_path)],
    )
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


def test_adapter_raises_typed_error_when_usalign_fails(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    path = tmp_path / "input.pdb"
    executable_path = tmp_path / "USalign"
    path.touch()
    executable_path.touch()
    monkeypatch.setattr(
        "structlens.integrations.usalign.adapter.subprocess.run",
        lambda *_args, **_kwargs: subprocess.CompletedProcess([], 2, "", "bad input"),
    )
    adapter = USAlignAdapter(
        executable=executable_path, structure_paths={"ref": path, "target": path}
    )

    with pytest.raises(USAlignExecutionError, match="bad input"):
        adapter.align(
            _chain("ref", "A"), _chain("target", "A"), StructuralAlignmentSettings()
        )
