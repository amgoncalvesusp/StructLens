"""Resolve the platform-pinned MUSCLE executable."""

from __future__ import annotations

import platform
from pathlib import Path


def bundled_executable(*, system: str | None = None, machine: str | None = None) -> Path | None:
    system_name = (system or platform.system()).lower()
    machine_name = (machine or platform.machine()).lower()
    if system_name.startswith("win"):
        key, filename = "windows-x64", "muscle.exe"
    elif system_name == "linux" and machine_name in {"x86_64", "amd64", "x64"}:
        key, filename = "linux-x64", "muscle"
    else:
        return None
    candidate = Path(__file__).parents[2] / "resources" / "muscle" / key / filename
    return candidate if candidate.is_file() else None


__all__ = ["bundled_executable"]
