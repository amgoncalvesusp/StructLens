"""Optional PyMOL integration; importing this module does not require PyMOL."""

from .launcher import LaunchResult, PyMOLLauncher

__all__ = ["LaunchResult", "PyMOLLauncher"]

from .selections import residue_selection, selection_name

__all__ = ["residue_selection", "selection_name"]
