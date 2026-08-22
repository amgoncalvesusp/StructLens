"""Construction and validation of authoritative residue correspondences."""

from .sequence_mapper import SequenceResidueMapper

__all__ = ["SequenceResidueMapper"]
from .validator import MappingDecision, choose_auto_mapping, validate_correspondence

__all__ = [
    "ManualResidueMapper",
    "MappingDecision",
    "choose_auto_mapping",
    "validate_correspondence",
]
from .manual_mapper import ManualResidueMapper
