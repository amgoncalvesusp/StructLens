"""Build the multi-resolution Windows icon from the approved StructLens mark."""

from __future__ import annotations

from pathlib import Path

from PIL import Image


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    source = root / "src" / "structlens" / "plugin" / "assets" / "structlens_icon.png"
    target = root / "packaging" / "windows" / "structlens.ico"
    image = Image.open(source).convert("RGBA")
    image.save(
        target,
        format="ICO",
        sizes=[(256, 256), (128, 128), (64, 64), (48, 48), (32, 32), (16, 16)],
    )


if __name__ == "__main__":
    main()
