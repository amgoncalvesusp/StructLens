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


class BundledBackendUnavailableError(StructLensError):
    """The packaged structural backend is unavailable for this platform."""


class UnsupportedPlatformError(StructLensError):
    """The requested operation is not supported on the current platform."""


class PyMOLNotConfiguredError(StructLensError):
    """PyMOL was not configured for an optional launch handoff."""


class PyMOLPluginUnavailableError(StructLensError):
    """The StructLens-PyMOL plugin is not available for an optional launch."""


class BundleValidationError(StructLensError):
    """A .structlens-pymol bundle is malformed or incomplete."""


class BundleCompatibilityError(BundleValidationError):
    """A .structlens-pymol bundle uses an unsupported schema major version."""


class UnsafeBundleError(BundleValidationError):
    """A bundle contains unsafe paths or executable content."""


class MultiStructureAlignmentError(StructLensError):
    """A multiple-structure analysis cannot be completed."""


__all__ = [
    "AnalysisCancelledError",
    "BundledBackendUnavailableError",
    "BundleCompatibilityError",
    "BundleValidationError",
    "ChainNotFoundError",
    "InputFormatError",
    "InsufficientAtomsError",
    "MappingError",
    "ProjectSchemaError",
    "PyMOLNotConfiguredError",
    "PyMOLPluginUnavailableError",
    "StructLensError",
    "MultiStructureAlignmentError",
    "UnsafeBundleError",
    "UnsupportedPlatformError",
    "USAlignExecutionError",
    "USAlignNotFoundError",
]
