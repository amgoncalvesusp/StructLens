"""Lazy Qt binding discovery for the optional PyMOL plugin surface."""

from __future__ import annotations

import os
from dataclasses import dataclass
from types import ModuleType


@dataclass(frozen=True, slots=True)
class QtBindings:
    """The small, binding-neutral module bundle used by the panel factory."""

    widgets: ModuleType
    core: ModuleType
    gui: ModuleType
    binding_name: str


def load_qt() -> QtBindings | None:
    """Return the requested Qt binding, preferring PySide6 when available."""

    requested_binding = os.environ.get("STRUCTLENS_QT_BINDING", "").lower()
    if requested_binding not in {"pyqt5", "qt5"}:
        try:
            from PySide6 import QtCore as pyside_core
            from PySide6 import QtGui as pyside_gui
            from PySide6 import QtWidgets as pyside_widgets

            return QtBindings(
                widgets=pyside_widgets,
                core=pyside_core,
                gui=pyside_gui,
                binding_name="PySide6",
            )
        except ImportError:
            pass

    if requested_binding not in {"pyside6", "qt6"}:
        try:
            from PyQt5 import QtCore as pyqt_core
            from PyQt5 import QtGui as pyqt_gui
            from PyQt5 import QtWidgets as pyqt_widgets

            return QtBindings(
                widgets=pyqt_widgets,
                core=pyqt_core,
                gui=pyqt_gui,
                binding_name="PyQt5",
            )
        except ImportError:
            pass
    return None


def signal_type(core: ModuleType) -> object:
    """Resolve the signal factory shared by PySide and PyQt."""

    signal = getattr(core, "Signal", None)
    if signal is not None:
        return signal
    return core.pyqtSignal


def user_role(core: ModuleType) -> object:
    """Resolve ``Qt.UserRole`` across Qt5 and Qt6 enum layouts."""

    qt = core.Qt
    legacy_role = getattr(qt, "UserRole", None)
    if legacy_role is not None:
        return legacy_role
    return qt.ItemDataRole.UserRole


__all__ = ["QtBindings", "load_qt", "signal_type", "user_role"]
