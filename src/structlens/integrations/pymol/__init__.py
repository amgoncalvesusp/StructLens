"""Optional PyMOL integration; importing this module does not require PyMOL."""

from .selections import residue_selection, selection_name

__all__ = ["residue_selection", "selection_name"]
