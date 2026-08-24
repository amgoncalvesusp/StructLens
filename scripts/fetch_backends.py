"""Fetch pinned offline backend assets for a release build.

The script is intentionally explicit: it never follows a floating latest URL
and it records the downloaded asset's SHA-256 next to the executable.
"""

from __future__ import annotations

import argparse
import hashlib
import platform
from pathlib import Path
from urllib.request import urlopen

MUSCLE_VERSION = "5.3"
MUSCLE_URLS = {
    "Windows": "https://github.com/rcedgar/muscle/releases/download/v5.3/muscle-win64.v5.3.exe",
    "Linux": "https://github.com/rcedgar/muscle/releases/download/v5.3/muscle-linux-x86.v5.3",
}


def download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with urlopen(url, timeout=60) as response:  # noqa: S310 - pinned HTTPS URL above
        destination.write_bytes(response.read())
    digest = hashlib.sha256(destination.read_bytes()).hexdigest()
    destination.with_suffix(destination.suffix + ".sha256").write_text(f"{digest}  {destination.name}\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("src/structlens/resources"))
    parser.add_argument("--platform", choices=("Windows", "Linux"), default=platform.system())
    args = parser.parse_args()
    if args.platform not in MUSCLE_URLS:
        raise SystemExit(f"Unsupported backend platform: {args.platform}")
    suffix = "muscle.exe" if args.platform == "Windows" else "muscle"
    target = args.output / "muscle" / ("windows-x64" if args.platform == "Windows" else "linux-x64") / suffix
    download(MUSCLE_URLS[args.platform], target)
    (args.output / "muscle" / "VERSION").write_text(MUSCLE_VERSION + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
