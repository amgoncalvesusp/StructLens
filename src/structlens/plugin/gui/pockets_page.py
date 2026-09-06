"""Immutable presentation records for Sites & Pockets evidence."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType

from structlens.core.evidence import Availability, Diagnostic
from structlens.core.pockets.matching import PocketMatchState
from structlens.core.pockets.volume_models import PocketVolumeResult
from structlens.core.reports.pockets import PocketReportSnapshot

from .method_help import pocket_volume_explanation


@dataclass(frozen=True, slots=True)
class PocketCandidateRow:
    candidate_id: str
    role: str
    sphere_count: int
    lining_residue_count: int
    centroid_xyz: tuple[float, float, float]
    state: PocketMatchState | None = None


@dataclass(frozen=True, slots=True)
class PocketMatchRow:
    state: PocketMatchState
    reference_candidate_id: str | None = None
    target_candidate_id: str | None = None
    lining_jaccard: float | None = None
    centroid_distance_angstrom: float | None = None
    score: float | None = None
    diagnostics: tuple[Diagnostic, ...] = ()


@dataclass(frozen=True, slots=True)
class PocketDetail:
    candidate_id: str
    role: str
    units: Mapping[str, str]
    coarse_volume_angstrom3: float | None = None
    fine_volume_angstrom3: float | None = None
    absolute_sensitivity_angstrom3: float | None = None
    relative_sensitivity: float | None = None
    coarse_voxel_count: int | None = None
    fine_voxel_count: int | None = None
    method_explanation: str = ""

    def __post_init__(self) -> None:
        object.__setattr__(self, "units", MappingProxyType(dict(self.units)))


@dataclass(frozen=True, slots=True)
class PocketPresentation:
    availability: Availability
    candidate_rows: tuple[PocketCandidateRow, ...] = ()
    match_rows: tuple[PocketMatchRow, ...] = ()
    details: tuple[PocketDetail, ...] = ()
    diagnostics: tuple[Diagnostic, ...] = ()
    units: Mapping[str, str] = field(default_factory=dict)
    provenance: object | None = None
    can_detect: bool = False
    can_measure: bool = False
    can_compare: bool = False
    disabled_reasons: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "candidate_rows", tuple(self.candidate_rows))
        object.__setattr__(self, "match_rows", tuple(self.match_rows))
        object.__setattr__(self, "details", tuple(self.details))
        object.__setattr__(self, "diagnostics", tuple(self.diagnostics))
        object.__setattr__(self, "units", MappingProxyType(dict(self.units)))

    def detail_for(self, candidate_id: str) -> PocketDetail | None:
        return next((item for item in self.details if item.candidate_id == candidate_id), None)


class PocketPresenter:
    """Pure projection of a canonical pocket report; no science is performed."""

    @staticmethod
    def present(snapshot: PocketReportSnapshot) -> PocketPresentation:
        if not isinstance(snapshot, PocketReportSnapshot):
            raise TypeError("snapshot must be a PocketReportSnapshot")
        matches = () if snapshot.matching is None else tuple(snapshot.matching.matches)
        match_by_id: dict[tuple[str, str], PocketMatchState] = {}
        match_rows = tuple(
            PocketMatchRow(
                state=match.state,
                reference_candidate_id=(match.reference_candidate.candidate_id if match.reference_candidate else None),
                target_candidate_id=(match.target_candidate.candidate_id if match.target_candidate else None),
                lining_jaccard=match.lining_jaccard,
                centroid_distance_angstrom=match.centroid_distance_angstrom,
                score=match.score,
                diagnostics=match.diagnostics,
            )
            for match in matches
        )
        for match in matches:
            if match.reference_candidate:
                match_by_id[("reference", match.reference_candidate.candidate_id)] = match.state
            if match.target_candidate:
                match_by_id[("target", match.target_candidate.candidate_id)] = match.state

        candidates = tuple(
            PocketCandidateRow(
                candidate_id=candidate.candidate_id,
                role=detection.role,
                sphere_count=len(candidate.alpha_spheres),
                lining_residue_count=len(candidate.lining_residues),
                centroid_xyz=candidate.centroid_xyz,
                state=match_by_id.get((detection.role, candidate.candidate_id)),
            )
            for detection in snapshot.detections
            for candidate in detection.candidates
        )
        volumes = {(item.role, item.result.candidate_id): item.result for item in snapshot.volumes if item.result and item.result.candidate_id}
        details = tuple(
            _detail_for(candidate, volumes.get((candidate.role, candidate.candidate_id)))
            for candidate in candidates
        )
        unavailable = snapshot.availability in {
            Availability.INVALID_INPUT,
            Availability.DEPENDENCY_UNAVAILABLE,
            Availability.NUMERICAL_FAILURE,
        }
        disabled: list[str] = []
        if unavailable:
            disabled.append(_capability_reason(snapshot))
        if not candidates:
            disabled.append("No detected pocket candidates are available for this report.")
        if not matches:
            disabled.append("Pocket comparison is unavailable until both structures provide candidates.")
        return PocketPresentation(
            availability=snapshot.availability or Availability.NOT_APPLICABLE,
            candidate_rows=candidates,
            match_rows=match_rows,
            details=details,
            diagnostics=snapshot.diagnostics,
            units=snapshot.units,
            provenance=snapshot.provenance,
            can_detect=not unavailable,
            can_measure=bool(candidates) and not unavailable,
            can_compare=bool(matches) and not unavailable,
            disabled_reasons=tuple(dict.fromkeys(disabled)),
        )


def _detail_for(row: PocketCandidateRow, result: PocketVolumeResult | None) -> PocketDetail:
    if result is None:
        return PocketDetail(row.candidate_id, row.role, {"volume": "Å³", "length": "Å"}, method_explanation=pocket_volume_explanation())
    return PocketDetail(
        candidate_id=row.candidate_id,
        role=row.role,
        units={"volume": "Å³", "length": "Å"},
        coarse_volume_angstrom3=result.coarse_volume_angstrom3,
        fine_volume_angstrom3=result.fine_volume_angstrom3,
        absolute_sensitivity_angstrom3=result.absolute_sensitivity_angstrom3,
        relative_sensitivity=result.relative_sensitivity,
        coarse_voxel_count=result.coarse_voxel_count,
        fine_voxel_count=result.fine_voxel_count,
        method_explanation=pocket_volume_explanation(),
    )


def _capability_reason(snapshot: PocketReportSnapshot) -> str:
    first = next(iter(snapshot.diagnostics), None)
    availability = snapshot.availability or Availability.NOT_APPLICABLE
    return first.message if first is not None else f"Pocket evidence is {availability.value}."


__all__ = ["PocketCandidateRow", "PocketDetail", "PocketMatchRow", "PocketPresentation", "PocketPresenter"]
