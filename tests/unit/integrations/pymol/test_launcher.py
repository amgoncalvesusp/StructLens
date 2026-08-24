from pathlib import Path

import pytest

from structlens.core.errors import PyMOLNotConfiguredError
from structlens.integrations.pymol.launcher import PyMOLLauncher


def test_launcher_prefers_explicit_executable(tmp_path: Path) -> None:
    executable = tmp_path / "pymol"
    executable.write_text("", encoding="utf-8")
    assert PyMOLLauncher(executable).locate() == executable


def test_launcher_reports_actionable_missing_configuration(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("STRUCTLENS_PYMOL_EXECUTABLE", raising=False)
    monkeypatch.setattr("structlens.integrations.pymol.launcher.shutil.which", lambda _name: None)
    with pytest.raises(PyMOLNotConfiguredError, match="Export|configured"):
        PyMOLLauncher().locate()
