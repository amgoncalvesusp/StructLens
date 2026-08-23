"""PyMOL plugin entry point with lazy host imports."""

from __future__ import annotations

from .gui.main_panel import build_qt_panel


def __init_plugin__() -> None:
    """Register the StructLens menu item when called by PyMOL."""

    try:
        from pymol import cmd, plugins  # type: ignore[import-not-found]
    except ImportError as exc:
        raise RuntimeError(
            "StructLens plugin must be loaded from a PyMOL environment"
        ) from exc

    plugins.addmenuitemqt("StructLens", lambda: build_qt_panel(command=cmd))


__all__ = ["__init_plugin__"]
