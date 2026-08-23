"""Typed failures exposed by StructLens at application boundaries."""


class StructLensError(Exception):
    """Base class for expected, user-actionable failures."""


class InputFormatError(StructLensError):
    """Input file is unreadable or unsupported."""


class ChainNotFoundError(StructLensError):
    """Requested chain does not exist in the normalized structure."""


class MappingError(StructLensError):
    """Residue mapping cannot satisfy the requested workflow."""


class AnalysisCancelledError(StructLensError):
    """An explicit user cancellation stopped an in-progress comparison."""


class InsufficientAtomsError(StructLensError):
    """Required atoms are absent for a requested geometry calculation."""


class USAlignNotFoundError(StructLensError):
    """Configured US-align executable cannot be found."""


class USAlignExecutionError(StructLensError):
    """US-align returned an execution failure."""


class ProjectSchemaError(StructLensError):
    """Project JSON does not match a supported schema."""


__all__ = [
    "AnalysisCancelledError",
    "ChainNotFoundError",
    "InputFormatError",
    "InsufficientAtomsError",
    "MappingError",
    "ProjectSchemaError",
    "StructLensError",
    "USAlignExecutionError",
    "USAlignNotFoundError",
]
