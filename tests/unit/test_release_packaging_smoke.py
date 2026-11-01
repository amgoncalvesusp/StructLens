"""Run the release entry point and keep package job inputs and metadata consistent."""

import tomllib
from io import BytesIO
from pathlib import Path

import pytest

from structlens import __version__
from structlens.cli.main import main as cli_main
from structlens.gui.main import main as gui_main

ROOT = Path(__file__).resolve().parents[2]


def test_tag_metadata_installer_notes_and_package_version_agree():
    version = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]["version"]
    assert __version__ == version
    assert f"version: {version}\n" in (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    iss = (ROOT / "packaging/windows/structlens.iss").read_text(encoding="utf-8")
    assert f'#define MyAppVersion "{version}"' in iss
    assert "RELEASE_NOTES_v{#MyAppVersion}.md" in iss
    assert "RELEASE_NOTES_v${VERSION}.md" in (ROOT / "packaging/linux/build_packages.sh").read_text(encoding="utf-8")


def test_release_build_checks_installed_gui_and_checksums():
    workflow = (ROOT / ".github/workflows/release.yml").read_text(encoding="utf-8")
    assert workflow.count('--hidden-import matplotlib.backends.backend_qtagg') == 2
    assert workflow.count('".[gui,charts]"') == 2
    for required in ("--cov-fail-under=80", "/VERYSILENT", "dpkg-deb -x", "--appimage-extract",
                     "smoke_release_package.py", "--verify-tag", "--notes-file", "sha256sum --check"):
        assert required in workflow, required


def test_backend_fetch_rejects_wrong_hash_without_replacing_existing_asset(monkeypatch, tmp_path):
    from scripts import fetch_backends

    dest = tmp_path / "muscle"
    dest.write_bytes(b"verified previous bytes")
    monkeypatch.setattr(fetch_backends, "urlopen", lambda *_args, **_kwargs: BytesIO(b"wrong contents"))
    with pytest.raises(ValueError, match="SHA-256"):
        fetch_backends.download(fetch_backends.MUSCLE_URLS["Linux"], dest)
    assert dest.read_bytes() == b"verified previous bytes"


def test_gui_release_probe_opens_scientific_report_and_refuses_invalid_input(monkeypatch, tmp_path):
    monkeypatch.setenv("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication

    monkeypatch.setattr(QApplication, "exec", lambda *_args: 92)
    structure = ROOT / "tests/fixtures/parsing/enriched.pdb"
    report = tmp_path / "report.json"
    assert cli_main(["compare", str(structure), str(structure), "--mode", "sequence", "--report", str(report)]) == 0
    assert gui_main(["StructLens", "--smoke-report", str(report), str(structure)]) == 0
    report.write_text("{}", encoding="utf-8")
    assert gui_main(["StructLens", "--smoke-report", str(report), str(structure)]) != 0
