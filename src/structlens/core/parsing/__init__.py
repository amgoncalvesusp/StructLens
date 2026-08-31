"""Input parsing and normalization for structure and sequence files."""

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
    "load_fasta",
    "load_mmcif",
    "load_pdb",
    "load_structure",
    "load_structure_evidence",
    "load_structure_legacy",
]
