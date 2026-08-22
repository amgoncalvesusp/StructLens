"""Safe local discovery and typed failures for the US-align executable."""

from __future__ import annotations

import shutil
from pathlib import Path

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


def discover_executable(explicit: str | Path | None = None) -> Path:
    """Find a pre-installed US-align executable without downloading anything."""

    if explicit is not None:
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

    for name in ("USalign", "US-align", "USalign.exe", "US-align.exe"):
        located = shutil.which(name)
        if located is not None:
            return Path(located)
    raise USAlignNotFoundError(
        "US-align executable was not found on PATH. Install US-align locally or "
        "configure its executable path."
    )


__all__ = [
    "USAlignError",
    "USAlignExecutionError",
    "USAlignNotFoundError",
    "USAlignOutputError",
    "discover_executable",
]
