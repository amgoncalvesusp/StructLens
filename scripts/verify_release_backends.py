"""Fail-fast smoke checks for bundled scientific backends."""

from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path


def check(executable: Path, *arguments: str) -> None:
    if not executable.is_file():
        raise SystemExit(f"Missing bundled backend: {executable}")
    if os.name != "nt":
        executable.chmod(executable.stat().st_mode | 0o111)
    result = subprocess.run([str(executable), *arguments], capture_output=True, text=True, shell=False, timeout=30)
    if result.returncode != 0:
        raise SystemExit(f"Backend failed smoke check: {executable}\n{result.stderr}")


def check_muscle_compliance(root: Path) -> None:
    """Refuse to ship the GPL-3 MUSCLE binary without its licence and source.

    Redistributing the executable obliges us to carry the complete licence text
    and the exact corresponding source alongside it. See
    docs/distribution-licensing.md.
    """

    muscle_root = root / "muscle"
    required = {
        "GPL-3 licence text": muscle_root / "LICENSE",
        "corresponding source archive": muscle_root / "muscle-5.3-source.tar.gz",
    }
    missing = [f"{label} ({path})" for label, path in required.items() if not path.is_file()]
    if missing:
        raise SystemExit(
            "MUSCLE compliance files are missing, so this build must not redistribute it:\n  "
            + "\n  ".join(missing)
            + "\nRun scripts/fetch_backends.py, or omit MUSCLE from the artifact."
        )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("src/structlens/resources"))
    parser.add_argument("--platform", choices=("windows", "linux"), required=True)
    args = parser.parse_args()
    check_muscle_compliance(args.root)
    if args.platform == "windows":
        check(args.root / "muscle" / "windows-x64" / "muscle.exe", "-version")
        check(args.root / "usalign" / "windows-x64" / "USalign.exe", "-h")
    else:
        check(args.root / "muscle" / "linux-x64" / "muscle", "-version")
        check(args.root / "usalign" / "linux-x64" / "USalign", "-h")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
