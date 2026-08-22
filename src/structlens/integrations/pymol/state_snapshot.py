"""Minimal reversible snapshot for StructLens-owned PyMOL names."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class PyMOLStateSnapshot:
    owned_selections: tuple[str, ...] = ()
    owned_objects: tuple[str, ...] = ()

    def restore(self, command: object) -> None:
        for name in self.owned_selections + self.owned_objects:
            delete = getattr(command, "delete", None)
            if delete is not None:
                delete(name)


__all__ = ["PyMOLStateSnapshot"]
