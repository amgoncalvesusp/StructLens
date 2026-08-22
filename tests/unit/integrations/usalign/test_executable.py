"""Discovery tests for the locally installed US-align executable."""

from __future__ import annotations

from pathlib import Path

import pytest

from structlens.integrations.usalign.executable import (
    USAlignNotFoundError,
    discover_executable,
)


def test_discovery_prefers_an_explicit_existing_executable(tmp_path: Path) -> None:
    executable = tmp_path / "USalign.exe"
    executable.touch()

    assert discover_executable(executable) == executable


def test_discovery_uses_path_when_no_explicit_executable(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        "structlens.integrations.usalign.executable.shutil.which",
        lambda name: "/tools/USalign" if name == "USalign" else None,
    )

    assert discover_executable() == Path("/tools/USalign")


def test_discovery_raises_actionable_error_for_missing_executable(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        "structlens.integrations.usalign.executable.shutil.which",
        lambda _name: None,
    )

    with pytest.raises(USAlignNotFoundError, match="PATH|executable"):
        discover_executable()
