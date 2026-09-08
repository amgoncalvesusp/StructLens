"""Run the GUI from this checkout, without installing over another version."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from structlens.gui.main import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
