from pathlib import Path

import pytest

from structlens.integrations.usalign.executable import platform_key, resolve_backend


def test_platform_key_normalizes_release_architectures() -> None:
    assert platform_key("Windows", "AMD64") == "windows-x64"
    assert platform_key("Linux", "aarch64") == "linux-arm64"
    assert platform_key("Darwin", "x86_64") == "macos-x64"


def test_custom_backend_reports_provenance(tmp_path: Path) -> None:
    executable = tmp_path / "USalign"
    executable.write_text("#!/bin/sh\n", encoding="utf-8")
    backend = resolve_backend(executable)
    assert backend.path == executable
    assert backend.source == "custom"


def test_unknown_platform_is_actionable() -> None:
    with pytest.raises(Exception, match="platform"):
        platform_key("Plan9", "x86_64")


def test_unknown_architecture_is_actionable() -> None:
    with pytest.raises(Exception, match="architecture"):
        platform_key("Linux", "riscv64")
