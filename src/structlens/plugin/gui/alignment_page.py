from .widgets import PageDescriptor

# Compatibility module: the integrated page is now named Structures.
PAGE = PageDescriptor("Structures", "Choose how StructLens should determine equivalent residues before superposition.")

__all__ = ["PAGE"]
