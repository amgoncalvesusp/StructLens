"""Typed failures exposed by StructLens at application boundaries."""


class StructLensError(Exception):
    """Base class for expected, user-actionable failures."""


class InputFormatError(StructLensError):
    pass


class ChainNotFoundError(StructLensError):
    pass


class MappingError(StructLensError):
    pass


class InsufficientAtomsError(StructLensError):
    pass


class USAlignNotFoundError(StructLensError):
    pass


class USAlignExecutionError(StructLensError):
    pass


class ProjectSchemaError(StructLensError):
    pass


__all__ = [
    "ChainNotFoundError",
    "InputFormatError",
    "InsufficientAtomsError",
    "MappingError",
    "ProjectSchemaError",
    "StructLensError",
    "USAlignExecutionError",
    "USAlignNotFoundError",
]
