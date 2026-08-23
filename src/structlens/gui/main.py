"""Launch the standalone StructLens Evidence Bench desktop application."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any


def main(argv: list[str] | None = None) -> int:
    """Start the file-based GUI without requiring a PyMOL host."""

    from structlens.plugin.gui.main_panel import build_qt_panel
    from structlens.plugin.gui.qt_compat import load_qt

    qt = load_qt()
    if qt is None:
        print(
            "StructLens GUI requires PySide6 or PyQt5. "
            "Install with: python -m pip install 'structlens[gui]'",
            file=sys.stderr,
        )
        return 2
    arguments = argv if argv is not None else sys.argv
    application = qt.widgets.QApplication.instance()
    if application is None:
        application = qt.widgets.QApplication(arguments)
    application.setApplicationName("StructLens")
    application.setApplicationDisplayName("StructLens")
    icon_path = Path(__file__).parents[1] / "plugin" / "assets" / "structlens_icon.png"
    if icon_path.exists():
        application.setWindowIcon(qt.gui.QIcon(str(icon_path)))
    panel: Any = build_qt_panel(command=None)
    panel.resize(1280, 820)
    panel.show()
    return int(application.exec())


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = ["main"]
