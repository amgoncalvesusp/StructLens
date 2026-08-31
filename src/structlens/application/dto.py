"""Stable request/summary DTOs for CLI and GUI adapters."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

from structlens.core.interactions import InteractionThresholds
from structlens.core.models import AnalysisResult, AnalysisSettings, ResidueId
from structlens.core.msa import MSASettings
from structlens.core.parsing import InputSelection, SourceSnapshot
from structlens.core.sites import SiteDefinition


@dataclass(frozen=True, slots=True)
class AnalysisRequest:
    reference_path: Path
    target_paths: tuple[Path, ...]
    settings: AnalysisSettings = AnalysisSettings()


@dataclass(frozen=True, slots=True)
class AnalysisReportRequest:
    """Content-addressed pairwise request for the canonical report workflow."""

    reference_snapshot: SourceSnapshot
    target_snapshot: SourceSnapshot
    reference_selection: InputSelection
    target_selection: InputSelection
    analysis_settings: AnalysisSettings = AnalysisSettings()
    msa_settings: MSASettings = MSASettings()
    interaction_thresholds: InteractionThresholds = InteractionThresholds()
    site_definitions: tuple[SiteDefinition, ...] = ()
    manual_pairs: tuple[tuple[ResidueId, ResidueId], ...] = ()
    minimum_vector_magnitude_angstrom: float = 0.5
    maximum_vectors: int = 100

    def __post_init__(self) -> None:
        if not isinstance(self.reference_snapshot, SourceSnapshot):
            raise TypeError("reference_snapshot must be a SourceSnapshot")
        if not isinstance(self.target_snapshot, SourceSnapshot):
            raise TypeError("target_snapshot must be a SourceSnapshot")
        if not isinstance(self.reference_selection, InputSelection):
            raise TypeError("reference_selection must be an InputSelection")
        if not isinstance(self.target_selection, InputSelection):
            raise TypeError("target_selection must be an InputSelection")
        if self.reference_selection.content_id != self.reference_snapshot.content_id:
            raise ValueError("reference selection content does not match its source snapshot")
        if self.target_selection.content_id != self.target_snapshot.content_id:
            raise ValueError("target selection content does not match its source snapshot")
        if self.reference_selection.format.value != self.reference_snapshot.logical_format:
            raise ValueError("reference selection format does not match its source snapshot")
        if self.target_selection.format.value != self.target_snapshot.logical_format:
            raise ValueError("target selection format does not match its source snapshot")
        if not isinstance(self.analysis_settings, AnalysisSettings):
            raise TypeError("analysis_settings must be AnalysisSettings")
        if not isinstance(self.msa_settings, MSASettings):
            raise TypeError("msa_settings must be MSASettings")
        if not isinstance(self.interaction_thresholds, InteractionThresholds):
            raise TypeError("interaction_thresholds must be InteractionThresholds")
        sites = tuple(self.site_definitions)
        if any(not isinstance(item, SiteDefinition) for item in sites):
            raise TypeError("site_definitions must contain SiteDefinition values")
        pairs = tuple(tuple(pair) for pair in self.manual_pairs)
        if any(
            len(pair) != 2 or not isinstance(pair[0], ResidueId) or not isinstance(pair[1], ResidueId)
            for pair in pairs
        ):
            raise TypeError("manual_pairs must contain reference/target ResidueId pairs")
        if len({pair[0] for pair in pairs}) != len(pairs) or len({pair[1] for pair in pairs}) != len(pairs):
            raise ValueError("manual_pairs must map each reference and target residue at most once")
        if self.analysis_settings.alignment_mode.value == "manual" and not pairs:
            raise ValueError("manual alignment mode requires manual_pairs")
        if self.analysis_settings.alignment_mode.value != "manual" and pairs:
            raise ValueError("manual_pairs require manual alignment mode")
        minimum = float(self.minimum_vector_magnitude_angstrom)
        if not math.isfinite(minimum) or minimum < 0.0:
            raise ValueError("minimum_vector_magnitude_angstrom must be finite and non-negative")
        if isinstance(self.maximum_vectors, bool) or not isinstance(self.maximum_vectors, int) or self.maximum_vectors < 1:
            raise ValueError("maximum_vectors must be a positive integer")
        object.__setattr__(self, "site_definitions", sites)
        object.__setattr__(self, "manual_pairs", pairs)
        object.__setattr__(self, "minimum_vector_magnitude_angstrom", minimum)


@dataclass(frozen=True, slots=True)
class TargetSummary:
    target_id: str
    sequence_identity: float
    sequence_coverage: float
    strict_rmsd_angstrom: float | None
    refined_rmsd_angstrom: float | None
    mutation_count: int

    @classmethod
    def from_result(cls, result: AnalysisResult) -> TargetSummary:
        return cls(
            result.target_id,
            result.sequence_identity,
            result.sequence_coverage,
            result.strict_rmsd_angstrom,
            result.refined_rmsd_angstrom,
            result.mutation_count,
        )


__all__ = ["AnalysisReportRequest", "AnalysisRequest", "TargetSummary"]
