"""Input parsing and normalization for structure and sequence files."""

from .coordinate_audit import (
    CoordinateAudit,
    CoordinateQualityError,
    audit_mmcif_coordinates,
    audit_mmcif_source,
    audit_pdb_coordinates,
    audit_pdb_source,
)
from .fasta import ProteinSequence, load_fasta
from .limits import ParseLimits, SnapshotError, SnapshotLimitError, StructureParseError
from .mmcif import load_mmcif
from .models import (
    AltlocPolicy,
    AssemblyScope,
    ChainLocator,
    InputSelection,
    ParsedStructure,
    StructureFormat,
    StructureMetadata,
)
from .normalize import load_structure, load_structure_evidence, load_structure_legacy
from .pdb import load_pdb
from .snapshot import SourceSnapshot, capture_snapshot

__all__ = [
    "ParseLimits",
    "CoordinateAudit",
    "CoordinateQualityError",
    "ProteinSequence",
    "AltlocPolicy",
    "AssemblyScope",
    "ChainLocator",
    "InputSelection",
    "ParsedStructure",
    "SnapshotError",
    "SnapshotLimitError",
    "StructureParseError",
    "SourceSnapshot",
    "StructureFormat",
    "StructureMetadata",
    "capture_snapshot",
    "audit_mmcif_coordinates",
    "audit_mmcif_source",
    "audit_pdb_coordinates",
    "audit_pdb_source",
    "load_fasta",
    "load_mmcif",
    "load_pdb",
    "load_structure",
    "load_structure_evidence",
    "load_structure_legacy",
]
