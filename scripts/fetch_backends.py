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
# MUSCLE is GPL-3. Redistributing the binary obliges us to ship its licence text
# and the exact corresponding source, so both are fetched from the same pinned
# tag as the executable. See docs/distribution-licensing.md.
MUSCLE_LICENSE_URL = "https://raw.githubusercontent.com/rcedgar/muscle/v5.3/LICENSE"
MUSCLE_SOURCE_URL = "https://github.com/rcedgar/muscle/archive/refs/tags/v5.3.tar.gz"
MUSCLE_SOURCE_ARCHIVE = "muscle-5.3-source.tar.gz"
MUSCLE_SHA256 = {
    MUSCLE_URLS["Linux"]: "318abeb951d786a3e2532714cc81ad3b3d8f79a2b517dc31316eeb5b694db2bc",
    MUSCLE_URLS["Windows"]: "a86637db99deb3efa1cb6762c4e3160ca33644608270eb167e4214f1f9409711",
    MUSCLE_LICENSE_URL: "3972dc9744f6499f0f9b2dbf76696f2ae7ad8af9b23dde66d6af86c9dfb36986",
    MUSCLE_SOURCE_URL: "74b22a94e630b16015c2bd9ae83aa2be2c2048d3e41f560b2d4a954725c81968",
}


def download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with urlopen(url, timeout=60) as response:  # noqa: S310 - pinned HTTPS URL above
        data = response.read(64 * 1024 * 1024 + 1)
    digest = hashlib.sha256(data).hexdigest()
    if digest != MUSCLE_SHA256[url]:
        raise ValueError(f"SHA-256 mismatch for pinned MUSCLE asset: {destination.name}")
    destination.write_bytes(data)
    destination.with_suffix(destination.suffix + ".sha256").write_text(
        f"{digest}  {destination.name}\n", encoding="utf-8"
    )


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
    muscle_root = args.output / "muscle"
    (muscle_root / "VERSION").write_text(f"{MUSCLE_VERSION}\n", encoding="utf-8")
    # Without these two files the artifact would redistribute a GPL-3 binary
    # with neither its licence nor its corresponding source.
    download(MUSCLE_LICENSE_URL, muscle_root / "LICENSE")
    download(MUSCLE_SOURCE_URL, muscle_root / MUSCLE_SOURCE_ARCHIVE)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
