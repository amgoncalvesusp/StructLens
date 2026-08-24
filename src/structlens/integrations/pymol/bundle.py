"""Compatibility import for the StructLens-PyMOL interchange contract."""

from structlens.integrations.pymol_bundle import (
    BUNDLE_FORMAT,
    BUNDLE_SCHEMA_VERSION,
    SUPPORTED_BUNDLE_SCHEMA_MAJOR,
    validate_pymol_bundle,
    write_pymol_bundle,
)

__all__ = [
    "BUNDLE_FORMAT",
    "BUNDLE_SCHEMA_VERSION",
    "SUPPORTED_BUNDLE_SCHEMA_MAJOR",
    "validate_pymol_bundle",
    "write_pymol_bundle",
]
