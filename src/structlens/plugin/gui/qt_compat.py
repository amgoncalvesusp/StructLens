"""Lazy Qt binding discovery for the optional PyMOL plugin surface."""

from __future__ import annotations

import importlib
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
    """Return PySide6 first, then PyQt5, without importing either at module load."""

    for binding in ("PySide6", "PyQt5"):
        try:
            return QtBindings(
                widgets=importlib.import_module(f"{binding}.QtWidgets"),
                core=importlib.import_module(f"{binding}.QtCore"),
                gui=importlib.import_module(f"{binding}.QtGui"),
                binding_name=binding,
            )
        except ImportError:
            continue
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
    return getattr(qt, "UserRole", qt.ItemDataRole.UserRole)


__all__ = ["QtBindings", "load_qt", "signal_type", "user_role"]
