"""Stable request/summary DTOs for CLI and GUI adapters."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from structlens.core.models import AnalysisResult, AnalysisSettings


@dataclass(frozen=True, slots=True)
class AnalysisRequest:
    reference_path: Path
    target_paths: tuple[Path, ...]
    settings: AnalysisSettings = AnalysisSettings()


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


__all__ = ["AnalysisRequest", "TargetSummary"]
