"""Launch the standalone StructLens Evidence Bench desktop application."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any


def _fit_window_to_screen(panel: Any) -> None:
    """Use Qt logical pixels and leave space for window decorations."""

    screen = panel.screen()
    if screen is None:
        panel.resize(1024, 680)
        return
    available = screen.availableGeometry()
    width = min(1280, max(1, available.width() - 32))
    height = min(820, max(1, available.height() - 64))
    panel.setMinimumSize(min(panel.minimumWidth(), width), min(panel.minimumHeight(), height))
    panel.resize(width, height)
    panel.move(available.x() + (available.width() - width) // 2, available.y() + (available.height() - height) // 2)


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
    _fit_window_to_screen(panel)
    panel.show()
    return int(application.exec())


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = ["main"]
