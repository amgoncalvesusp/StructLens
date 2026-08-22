from __future__ import annotations

from pathlib import Path

import pytest
from PIL import Image

from structlens.application.export_service import export_publication_image


def test_publication_image_export_sets_dimensions_and_dpi(tmp_path: Path) -> None:
    image = Image.new("RGB", (40, 40), "white")
    output = tmp_path / "figure.tiff"

    export_publication_image(image, output, width_mm=25.4, dpi=600)
    saved = Image.open(output)

    assert saved.size == (600, 600)
    assert saved.info["dpi"] == pytest.approx((600, 600), abs=1)


def test_jpeg_export_rejects_transparent_input(tmp_path: Path) -> None:
    image = Image.new("RGBA", (20, 20), (255, 255, 255, 0))

    with pytest.raises(ValueError, match="transparency"):
        export_publication_image(image, tmp_path / "figure.jpg", dpi=300)
