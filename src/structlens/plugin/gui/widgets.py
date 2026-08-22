"""Small GUI-independent descriptors used by optional Qt pages."""

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class PageDescriptor:
    name: str
    purpose: str


__all__ = ["PageDescriptor"]
