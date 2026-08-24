"""Safe local discovery and typed failures for the US-align executable."""

from __future__ import annotations

import platform
import shutil
from dataclasses import dataclass
from pathlib import Path

from structlens.core.errors import (
    BundledBackendUnavailableError,
)
from structlens.core.errors import (
    USAlignExecutionError as CoreUSAlignExecutionError,
)
from structlens.core.errors import (
    USAlignNotFoundError as CoreUSAlignNotFoundError,
)


class USAlignError(Exception):
    """Base class for US-align adapter failures."""


class USAlignNotFoundError(CoreUSAlignNotFoundError, USAlignError):
    """Raised when US-align is neither configured nor available on PATH."""


class USAlignExecutionError(CoreUSAlignExecutionError, USAlignError):
    """Raised when a US-align process cannot produce a valid alignment."""


class USAlignOutputError(USAlignExecutionError):
    """Raised when successful US-align output lacks required result data."""


@dataclass(frozen=True, slots=True)
class USAlignBackend:
    """Resolved executable plus provenance used in analysis records."""

    path: Path
    version: str
    source: str
    platform: str


def platform_key(
    system: str | None = None,
    machine: str | None = None,
) -> str:
    """Return the resource directory key for a supported release platform."""

    system_name = (system or platform.system()).lower()
    machine_name = (machine or platform.machine()).lower().replace("amd64", "x86_64")
    if machine_name in {"aarch64", "arm64"}:
        architecture = "arm64"
    elif machine_name in {"x86_64", "x86-64", "x64"}:
        architecture = "x64"
    else:
        raise BundledBackendUnavailableError(
            f"The bundled US-align backend is not available for architecture "
            f"{machine or platform.machine()} on {system or platform.system()}."
        )
    if system_name.startswith("win"):
        return f"windows-{architecture}"
    if system_name == "darwin":
        return f"macos-{architecture}"
    if system_name == "linux":
        return f"linux-{architecture}"
    raise BundledBackendUnavailableError(
        f"The bundled US-align backend is not available for platform {system or platform.system()} {machine or platform.machine()}."
    )


def bundled_executable(
    *,
    system: str | None = None,
    machine: str | None = None,
) -> Path | None:
    """Locate a shipped binary, returning ``None`` for source-only installs."""

    key = platform_key(system, machine)
    filename = "USalign.exe" if key.startswith("windows-") else "USalign"
    candidate = Path(__file__).parents[2] / "resources" / "usalign" / key / filename
    return candidate if candidate.is_file() else None


def resolve_backend(
    explicit: str | Path | None = None,
    *,
    allow_path_fallback: bool = True,
) -> USAlignBackend:
    """Resolve custom, bundled, then diagnostic PATH backends in that order."""

    if explicit is not None:
        path = _resolve_explicit(explicit)
        return USAlignBackend(path, _read_version(path), "custom", platform_key())
    bundled = bundled_executable()
    if bundled is not None:
        _ensure_executable(bundled)
        return USAlignBackend(bundled, _read_version(bundled), "bundled", platform_key())
    if allow_path_fallback:
        path_on_path = _path_executable()
        if path_on_path is not None:
            return USAlignBackend(path_on_path, _read_version(path_on_path), "PATH", platform_key())
    raise USAlignNotFoundError(
        "The bundled US-align backend is unavailable for this installation. "
        "Open Diagnostics to verify the StructLens installation or choose a "
        "custom US-align executable in Advanced Settings."
    )


def discover_executable(explicit: str | Path | None = None) -> Path:
    """Find a custom, bundled, or diagnostic PATH executable."""

    return resolve_backend(explicit).path


def _resolve_explicit(explicit: str | Path) -> Path:

    configured = Path(explicit)
    if configured.is_file():
        return configured
    located = shutil.which(str(explicit))
    if located is not None:
        return Path(located)
    raise USAlignNotFoundError(
        f"US-align executable '{explicit}' was not found. Configure an existing "
        "executable path or add USalign to PATH."
    )


def _path_executable() -> Path | None:
    for name in ("USalign", "US-align", "USalign.exe", "US-align.exe"):
        located = shutil.which(name)
        if located is not None:
            return Path(located)
    return None


def _read_version(path: Path) -> str:
    version_file = path.parent.parent / "VERSION"
    if version_file.is_file():
        text = version_file.read_text(encoding="utf-8").strip()
        if text:
            return text
    return "unknown"


def _ensure_executable(path: Path) -> None:
    if platform.system() != "Windows":
        path.chmod(path.stat().st_mode | 0o111)


__all__ = [
    "USAlignError",
    "USAlignExecutionError",
    "USAlignNotFoundError",
    "USAlignOutputError",
    "USAlignBackend",
    "bundled_executable",
    "discover_executable",
    "platform_key",
    "resolve_backend",
]
