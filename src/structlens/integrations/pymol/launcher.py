"""Optional external PyMOL launch handoff without shell interpolation."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from structlens.core.errors import PyMOLNotConfiguredError, PyMOLPluginUnavailableError


@dataclass(frozen=True, slots=True)
class LaunchResult:
    executable: Path
    bundle_path: Path
    process_id: int | None


class PyMOLLauncher:
    """Locate PyMOL and launch a bundle using a safe argument vector."""

    def __init__(self, executable: str | Path | None = None) -> None:
        self.executable = executable

    def locate(self) -> Path:
        configured = self.executable or os.environ.get("STRUCTLENS_PYMOL_EXECUTABLE")
        if configured:
            path = Path(configured)
            if path.is_file():
                return path
            located = shutil.which(str(configured))
            if located:
                return Path(located)
            raise PyMOLNotConfiguredError(f"Configured PyMOL executable was not found: {configured}")
        for name in ("pymol", "pymol.exe", "PyMOL.exe"):
            located = shutil.which(name)
            if located:
                return Path(located)
        raise PyMOLNotConfiguredError(
            "PyMOL is not configured. Export the validated bundle and open it from the PyMOL plugin."
        )

    def launch_bundle(self, bundle_path: str | Path, *, extra_args: Sequence[str] = ()) -> LaunchResult:
        path = Path(bundle_path)
        if not path.is_file():
            raise PyMOLNotConfiguredError(f"PyMOL bundle does not exist: {path}")
        executable = self.locate()
        startup = f"structlens_open {json.dumps(str(path))}"
        command = [str(executable), "-d", startup, *extra_args]
        try:
            process = subprocess.Popen(command, shell=False)
        except OSError as error:
            raise PyMOLPluginUnavailableError(f"Could not start PyMOL: {error}") from error
        return LaunchResult(executable, path, process.pid)


__all__ = ["LaunchResult", "PyMOLLauncher"]
