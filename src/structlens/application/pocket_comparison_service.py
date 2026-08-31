"""Application orchestration for mutation-associated pocket comparison."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, TypeAlias, cast

from structlens.core.difference_maps import ResidueDisplacementVector
from structlens.core.evidence import (
    Availability,
    Diagnostic,
    DiagnosticSeverity,
    PocketEvidenceChannel,
    PocketEvidenceSections,
)
from structlens.core.interactions import InteractionDifference
from structlens.core.models import MutationEvent, ResidueCorrespondence, ResidueId, StructuralTransform
from structlens.core.pockets.comparison import PocketComparison, PocketSurfaceMeasurement, compare_pocket_match
from structlens.core.pockets.focused import FocusedPocketSelection
from structlens.core.pockets.matching import (
    MatchInput,
    PocketMatch,
    PocketMatchingResult,
    PocketMatchingSettings,
    match_pocket_candidates,
)
from structlens.core.pockets.models import PocketCandidate
from structlens.core.pockets.volume_models import PocketVolumeResult, validate_pocket_volume_result

CandidateInput: TypeAlias = Sequence[PocketCandidate] | PocketCandidate | FocusedPocketSelection
VolumeInput: TypeAlias = PocketVolumeResult
SurfaceInput: TypeAlias = PocketSurfaceMeasurement

MAX_INPUT_COLLECTION_SIZE = 100_000


@dataclass(frozen=True, slots=True)
class PocketComparisonReport:
    """Immutable collection returned by :class:`PocketComparisonService`."""

    availability: Availability
    matching: PocketMatchingResult
    comparisons: tuple[PocketComparison, ...] = field(default_factory=tuple)
    diagnostics: tuple[Diagnostic, ...] = field(default_factory=tuple)
    concordance: PocketEvidenceSections | None = None

    def __post_init__(self) -> None:
        availability = self.availability
        if not isinstance(availability, Availability):
            availability = Availability(availability)
            object.__setattr__(self, "availability", availability)
        if not isinstance(self.matching, PocketMatchingResult):
            raise TypeError("matching must be a PocketMatchingResult")
        if not isinstance(self.comparisons, Sequence) or isinstance(self.comparisons, (str, bytes)):
            raise TypeError("comparisons must be a bounded Sequence")
        if len(self.comparisons) > MAX_INPUT_COLLECTION_SIZE:
            raise ValueError("comparisons exceed the bounded application input limit")
        comparisons = tuple(self.comparisons)
        if any(not isinstance(item, PocketComparison) for item in comparisons):
            raise TypeError("comparisons must contain PocketComparison values")
        if not isinstance(self.diagnostics, Sequence) or isinstance(self.diagnostics, (str, bytes)):
            raise TypeError("diagnostics must be a bounded Sequence")
        if len(self.diagnostics) > MAX_INPUT_COLLECTION_SIZE:
            raise ValueError("diagnostics exceed the bounded application input limit")
        diagnostics = tuple(self.diagnostics)
        if any(not isinstance(item, Diagnostic) for item in diagnostics):
            raise TypeError("diagnostics must contain Diagnostic values")
        object.__setattr__(self, "comparisons", comparisons)
        object.__setattr__(self, "diagnostics", diagnostics)
        if self.concordance is not None and not isinstance(self.concordance, PocketEvidenceSections):
            raise TypeError("concordance must be PocketEvidenceSections or None")

    @property
    def status(self) -> Availability:
        return self.availability

    @property
    def matches(self) -> tuple[PocketMatch, ...]:
        return self.matching.matches

    @property
    def pocket_comparisons(self) -> tuple[PocketComparison, ...]:
        return self.comparisons

    @property
    def pocket_concordance(self) -> PocketEvidenceSections | None:
        return self.concordance

    def to_json(self) -> dict[str, object]:
        return {
            "availability": self.availability.value,
            "matching": self.matching.to_json(),
            "comparisons": [item.to_json() for item in self.comparisons],
            "diagnostics": [item.to_json() for item in self.diagnostics],
            "concordance": self.concordance.to_json() if self.concordance is not None else None,
        }


def _candidates(
    value: CandidateInput,
    *,
    maximum_candidate_count: int | None = None,
) -> tuple[tuple[PocketCandidate, ...], tuple[Diagnostic, ...], Availability]:
    if isinstance(value, FocusedPocketSelection):
        if value.availability is Availability.AVAILABLE and value.candidate is not None:
            return (value.candidate,), value.diagnostics, Availability.AVAILABLE
        return (), value.diagnostics, value.availability
    if isinstance(value, PocketCandidate):
        return (value,), (), Availability.AVAILABLE
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise TypeError("candidate inputs must be PocketCandidate or bounded Sequence")
    candidate_limit = (
        min(MAX_INPUT_COLLECTION_SIZE, maximum_candidate_count)
        if maximum_candidate_count is not None
        else MAX_INPUT_COLLECTION_SIZE
    )
    over_candidate_limit = len(value) > candidate_limit
    # Built-in collections are already materialized and can still be inspected
    # for candidate collisions; custom oversized Sequences fail closed without
    # invoking iteration.  The hard application cap remains absolute.
    if over_candidate_limit and (not isinstance(value, (tuple, list)) or len(value) > MAX_INPUT_COLLECTION_SIZE):
        return (
            (),
            (
                Diagnostic(
                    "pocket.compare.resource_limit",
                    DiagnosticSeverity.ERROR,
                    "Candidate collection exceeds the bounded application input limit.",
                ),
            ),
            Availability.NUMERICAL_FAILURE,
        )
    values = tuple(value)
    if any(not isinstance(item, PocketCandidate) for item in values):
        raise TypeError("candidate inputs must contain PocketCandidate values")
    return values, (), Availability.AVAILABLE


def _validate_lineage(candidates: Sequence[PocketCandidate], side: str) -> None:
    if any(item.source_content_id is None or item.selection_id is None for item in candidates):
        raise ValueError(f"{side} candidates require complete source/selection lineage")
    if len({item.candidate_id for item in candidates}) != len(candidates):
        raise ValueError(f"{side} candidate_id collision makes ordering ambiguous")
    if len({item.selection_id for item in candidates}) > 1:
        raise ValueError(f"{side} candidates must have homogeneous selection lineage")
    if len({item.source_content_id for item in candidates}) > 1:
        raise ValueError(f"{side} candidates must have homogeneous source lineage")


def _validate_volume_result(value: PocketVolumeResult, candidate: PocketCandidate) -> None:
    if not isinstance(value, PocketVolumeResult) or not isinstance(candidate, PocketCandidate):
        raise TypeError("volume result and candidate must use typed contracts")
    if candidate.source_content_id is None or candidate.selection_id is None:
        raise ValueError("candidate requires complete source/selection lineage")
    if value.candidate_id != candidate.candidate_id:
        raise ValueError("pocket volume result candidate identity does not match candidate")
    if value.provenance is None:
        raise ValueError("pocket volume result requires provenance at application boundary")
    inputs = value.provenance.input_hashes
    if inputs.get("candidate") != candidate.candidate_id:
        raise ValueError("pocket volume provenance does not identify candidate")
    if inputs.get("source_content") != candidate.source_content_id:
        raise ValueError("pocket volume provenance source hash does not match candidate lineage")
    parameters = value.provenance.parameters
    if parameters.get("candidate_id") != candidate.candidate_id:
        raise ValueError("pocket volume provenance candidate_id does not match candidate")
    lineage = parameters.get("candidate_lineage")
    if not isinstance(lineage, Mapping):
        raise ValueError("pocket volume provenance requires complete candidate_lineage")
    expected_lineage = {
        "candidate_id": candidate.candidate_id,
        "source_content_id": candidate.source_content_id,
        "selection_id": candidate.selection_id,
    }
    if dict(lineage) != expected_lineage:
        raise ValueError("pocket volume provenance candidate_lineage does not match candidate")
    for name in ("source_content_id", "selection_id"):
        if parameters.get(name) != getattr(candidate, name):
            raise ValueError("pocket volume provenance does not match candidate lineage")
    result_parameters = value.provenance_parameters
    if not isinstance(result_parameters, Mapping):
        raise ValueError("pocket volume result requires complete provenance parameters")
    if (value.provenance.method_id, value.provenance.method_version) == ("structlens.pocket.volume", "0.4.0"):
        if result_parameters.get("raw_source_hash") != inputs.get("raw_source"):
            raise ValueError("pocket volume raw-source snapshot does not match method provenance")
        if (
            result_parameters.get("logical_content") != inputs.get("logical_content")
            or result_parameters.get("logical_content") != candidate.source_content_id
        ):
            raise ValueError("pocket volume logical-content hash snapshot does not match candidate provenance")
    for name, expected in expected_lineage.items():
        if result_parameters.get(name) != expected:
            raise ValueError("pocket volume result provenance does not match candidate lineage")
    expected_sphere_metadata: dict[str, object] = {
        "sphere_ids": tuple(sphere.sphere_id for sphere in candidate.alpha_spheres),
        "sphere_radii_angstrom": tuple(sphere.radius_angstrom for sphere in candidate.alpha_spheres),
        "sphere_count": len(candidate.alpha_spheres),
    }
    for label, candidate_parameters in (("result", result_parameters), ("method", parameters)):
        for sphere_name, sphere_expected in expected_sphere_metadata.items():
            if candidate_parameters.get(sphere_name) != sphere_expected:
                raise ValueError(f"pocket volume {label} {sphere_name} does not match candidate sphere geometry")
    validate_pocket_volume_result(value)


def _volume_for(
    value: Mapping[str, VolumeInput] | VolumeInput | None, candidate: PocketCandidate | None
) -> VolumeInput | None:
    if value is None or candidate is None:
        return None
    if isinstance(value, PocketVolumeResult):
        _validate_volume_result(value, candidate)
        return value
    if not isinstance(value, Mapping):
        raise TypeError("volume inputs must be typed volume result or candidate-id mapping")
    selected = value.get(candidate.candidate_id)
    if selected is not None and not isinstance(selected, PocketVolumeResult):
        raise TypeError("volume mappings must contain pocket volume values")
    if isinstance(selected, PocketVolumeResult):
        _validate_volume_result(selected, candidate)
    return selected


def _surface_for(
    value: Mapping[str, SurfaceInput] | SurfaceInput | None, candidate: PocketCandidate | None, name: str
) -> PocketSurfaceMeasurement | None:
    if value is None or candidate is None:
        return None
    selected: object | None
    if isinstance(value, PocketSurfaceMeasurement):
        selected = value
    elif isinstance(value, Mapping):
        selected = value.get(candidate.candidate_id)
    else:
        raise TypeError(f"{name} must contain typed surface measurements")
    if selected is None:
        return None
    if not isinstance(selected, PocketSurfaceMeasurement):
        raise TypeError(f"{name} must contain typed surface measurements")
    if (
        selected.candidate_id != candidate.candidate_id
        or selected.source_content_id != candidate.source_content_id
        or selected.selection_id != candidate.selection_id
    ):
        raise ValueError(f"{name} measurement identity does not match candidate lineage")
    return selected


def _validate_volume_input(value: object, name: str) -> None:
    if value is None or isinstance(value, PocketVolumeResult):
        return
    if not isinstance(value, Mapping):
        raise TypeError(f"{name} must contain PocketVolumeResult typed measurements")
    if len(value) > MAX_INPUT_COLLECTION_SIZE:
        raise ValueError(f"{name} exceeds the bounded application input limit")
    for key, item in value.items():
        if not isinstance(key, str) or not isinstance(item, PocketVolumeResult):
            raise TypeError(f"{name} must map candidate IDs to typed volume measurements")


def _validate_surface_input(value: object, name: str) -> None:
    if value is None or isinstance(value, PocketSurfaceMeasurement):
        return
    if not isinstance(value, Mapping):
        raise TypeError(f"{name} must contain typed surface measurements")
    if len(value) > MAX_INPUT_COLLECTION_SIZE:
        raise ValueError(f"{name} exceeds the bounded application input limit")
    for key, item in value.items():
        if not isinstance(key, str) or not isinstance(item, PocketSurfaceMeasurement):
            raise TypeError(f"{name} must map candidate IDs to typed surface measurements")


def _validate_displacement_input(value: object) -> None:
    if value is None:
        return
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise TypeError("local_displacements must contain typed displacement vectors")
    if len(value) > MAX_INPUT_COLLECTION_SIZE:
        raise ValueError("local_displacements exceeds the bounded application input limit")
    if any(not isinstance(item, ResidueDisplacementVector) for item in value):
        raise TypeError("local_displacements must contain typed displacement vectors")
    residues = tuple(item.reference_residue for item in value)
    if len(set(residues)) != len(residues):
        raise ValueError("local_displacements must not contain duplicate residues")


def _validate_displacement_lineage(value: object, correspondences: object) -> None:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        return
    if isinstance(correspondences, Mapping):
        forward = dict(correspondences)
    else:
        if not isinstance(correspondences, Sequence) or isinstance(correspondences, (str, bytes)):
            return
        forward = {
            item.reference: item.target
            for item in correspondences
            if isinstance(item, ResidueCorrespondence) and item.reference is not None and item.target is not None
        }
    for item in value:
        if isinstance(item, ResidueDisplacementVector) and item.reference_residue in forward:
            if item.target_residue != forward[item.reference_residue]:
                raise ValueError("local displacement target does not match authoritative correspondence")


def _validate_qc_input(value: object, name: str) -> None:
    if value is None or isinstance(value, Availability):
        return
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        if len(value) > MAX_INPUT_COLLECTION_SIZE:
            raise ValueError(f"{name} exceeds the bounded application input limit")
        if any(not isinstance(item, Diagnostic) for item in value):
            raise TypeError(f"{name} diagnostics must contain Diagnostic values")
        return
    if not hasattr(value, "availability"):
        raise TypeError(f"{name} must expose typed availability")
    try:
        Availability(cast(Any, value).availability)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{name} must expose typed availability") from exc
    diagnostics = cast(Any, value).diagnostics if hasattr(value, "diagnostics") else ()
    if not isinstance(diagnostics, Sequence) or isinstance(diagnostics, (str, bytes)):
        raise TypeError(f"{name} diagnostics must be a bounded Sequence")
    if len(diagnostics) > MAX_INPUT_COLLECTION_SIZE:
        raise ValueError(f"{name} diagnostics exceed the bounded application input limit")
    if any(not isinstance(item, Diagnostic) for item in diagnostics):
        raise TypeError(f"{name} diagnostics must contain Diagnostic values")


def _validate_diagnostics_input(value: object, name: str) -> None:
    if value is None:
        return
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise TypeError(f"{name} must be a bounded Sequence")
    if len(value) > MAX_INPUT_COLLECTION_SIZE:
        raise ValueError(f"{name} exceeds the bounded application input limit")
    if any(not isinstance(item, Diagnostic) for item in value):
        raise TypeError(f"{name} must contain Diagnostic values")


def _validate_mapping_identity(value: object, candidates: Sequence[PocketCandidate], name: str) -> None:
    if not isinstance(value, Mapping):
        return
    candidate_ids = {item.candidate_id for item in candidates}
    unknown = set(value) - candidate_ids
    if unknown:
        raise ValueError(f"{name} contains evidence for an unknown candidate identity")


def _validate_correspondence_input(value: object) -> None:
    if isinstance(value, Mapping):
        if len(value) > MAX_INPUT_COLLECTION_SIZE:
            raise ValueError("correspondences exceeds the bounded application input limit")
        pairs = tuple(value.items())
        if any(not isinstance(left, ResidueId) or not isinstance(right, ResidueId) for left, right in pairs):
            raise TypeError("correspondences must contain ResidueId pairs")
        if len({left for left, _ in pairs}) != len(pairs) or len({right for _, right in pairs}) != len(pairs):
            raise ValueError("correspondences must be bijective")
        return
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise TypeError("correspondences must be a bounded Sequence or Mapping")
    if len(value) > MAX_INPUT_COLLECTION_SIZE:
        raise ValueError("correspondences exceeds the bounded application input limit")
    rows = tuple(value)
    if any(not isinstance(item, ResidueCorrespondence) for item in rows):
        raise TypeError("correspondences must contain ResidueCorrespondence values")
    references = tuple(item.reference for item in rows if item.reference is not None)
    targets = tuple(item.target for item in rows if item.target is not None)
    if len(set(references)) != len(references) or len(set(targets)) != len(targets):
        raise ValueError("correspondences must be bijective")


def _validate_sequence_contract(value: object, name: str, item_type: type[object]) -> None:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise TypeError(f"{name} must be a bounded Sequence")
    if len(value) > MAX_INPUT_COLLECTION_SIZE:
        raise ValueError(f"{name} exceeds the bounded application input limit")
    if any(not isinstance(item, item_type) for item in value):
        raise TypeError(f"{name} must contain {item_type.__name__} values")


def _validate_bound_volume(value: object, candidates: Sequence[PocketCandidate], name: str) -> None:
    if value is not None and not isinstance(value, (PocketVolumeResult, Mapping)):
        raise TypeError(f"{name} must use PocketVolumeResult evidence bound to a candidate")
    if isinstance(value, PocketVolumeResult):
        if len(candidates) != 1:
            raise ValueError(f"{name} direct volume requires one unambiguous candidate")
        _validate_volume_result(value, candidates[0])
    elif isinstance(value, Mapping):
        for key, item in value.items():
            candidate = next((candidate for candidate in candidates if candidate.candidate_id == key), None)
            _validate_volume_result(item, candidate)  # type: ignore[arg-type]


def _validate_bound_surface(value: object, candidates: Sequence[PocketCandidate], name: str) -> None:
    if isinstance(value, PocketSurfaceMeasurement):
        if len(candidates) != 1:
            raise ValueError(f"{name} direct surface requires one unambiguous candidate")
        if (
            value.candidate_id != candidates[0].candidate_id
            or value.source_content_id != candidates[0].source_content_id
            or value.selection_id != candidates[0].selection_id
        ):
            raise ValueError(f"{name} measurement identity does not match candidate lineage")
    elif isinstance(value, Mapping):
        for key, item in value.items():
            candidate = next((candidate for candidate in candidates if candidate.candidate_id == key), None)
            _surface_for(item, candidate, name)


class PocketComparisonService:
    """Coordinate matching and comparison of selected pocket evidence."""

    def __init__(self, settings: PocketMatchingSettings | None = None) -> None:
        if settings is not None and not isinstance(settings, PocketMatchingSettings):
            raise TypeError("settings must be PocketMatchingSettings or None")
        self._settings = settings or PocketMatchingSettings()

    @property
    def settings(self) -> PocketMatchingSettings:
        return self._settings

    def compare(
        self,
        reference_candidates: CandidateInput,
        target_candidates: CandidateInput,
        correspondences: MatchInput = (),
        *,
        transform: StructuralTransform | None = None,
        settings: PocketMatchingSettings | None = None,
        reference_volumes: Mapping[str, VolumeInput] | VolumeInput | None = None,
        target_volumes: Mapping[str, VolumeInput] | VolumeInput | None = None,
        reference_volume: VolumeInput | None = None,
        target_volume: VolumeInput | None = None,
        reference_surface_areas_angstrom2: Mapping[str, SurfaceInput] | SurfaceInput | None = None,
        target_surface_areas_angstrom2: Mapping[str, SurfaceInput] | SurfaceInput | None = None,
        reference_surface_area_angstrom2: SurfaceInput | None = None,
        target_surface_area_angstrom2: SurfaceInput | None = None,
        mutations: Sequence[MutationEvent] = (),
        interaction_differences: Sequence[InteractionDifference] = (),
        local_displacements: Mapping[ResidueId, float] | Sequence[object] | None = None,
        qc: Availability | Sequence[Diagnostic] | object | None = None,
        reference_qc: Availability | Sequence[Diagnostic] | object | None = None,
        target_qc: Availability | Sequence[Diagnostic] | object | None = None,
        qc_diagnostics: Sequence[Diagnostic] | None = None,
    ) -> PocketComparisonReport:
        """Compare candidates while preserving unmatched and unavailable states."""

        run_settings = settings if settings is not None else self._settings
        if not isinstance(run_settings, PocketMatchingSettings):
            raise TypeError("settings must be PocketMatchingSettings or None")
        if transform is not None and not isinstance(transform, StructuralTransform):
            raise TypeError("transform must be StructuralTransform or None")
        if reference_volumes is not None and reference_volume is not None:
            raise ValueError("provide either reference_volumes or reference_volume, not both")
        if target_volumes is not None and target_volume is not None:
            raise ValueError("provide either target_volumes or target_volume, not both")
        if reference_surface_areas_angstrom2 is not None and reference_surface_area_angstrom2 is not None:
            raise ValueError(
                "provide either reference_surface_areas_angstrom2 or reference_surface_area_angstrom2, not both"
            )
        if target_surface_areas_angstrom2 is not None and target_surface_area_angstrom2 is not None:
            raise ValueError("provide either target_surface_areas_angstrom2 or target_surface_area_angstrom2, not both")
        _validate_correspondence_input(correspondences)
        _validate_sequence_contract(mutations, "mutations", MutationEvent)
        _validate_sequence_contract(interaction_differences, "interaction_differences", InteractionDifference)
        _validate_displacement_input(local_displacements)
        _validate_displacement_lineage(local_displacements, correspondences)
        _validate_qc_input(qc, "qc")
        _validate_qc_input(reference_qc, "reference_qc")
        _validate_qc_input(target_qc, "target_qc")
        _validate_diagnostics_input(qc_diagnostics, "qc_diagnostics")
        _validate_volume_input(
            reference_volumes if reference_volumes is not None else reference_volume, "reference volumes"
        )
        _validate_volume_input(target_volumes if target_volumes is not None else target_volume, "target volumes")
        _validate_surface_input(
            reference_surface_areas_angstrom2
            if reference_surface_areas_angstrom2 is not None
            else reference_surface_area_angstrom2,
            "reference surfaces",
        )
        _validate_surface_input(
            target_surface_areas_angstrom2
            if target_surface_areas_angstrom2 is not None
            else target_surface_area_angstrom2,
            "target surfaces",
        )
        references, reference_diagnostics, reference_state = _candidates(
            reference_candidates, maximum_candidate_count=run_settings.maximum_candidate_count
        )
        targets, target_diagnostics, target_state = _candidates(
            target_candidates, maximum_candidate_count=run_settings.maximum_candidate_count
        )
        if reference_state is Availability.AVAILABLE:
            _validate_lineage(references, "reference")
        if target_state is Availability.AVAILABLE:
            _validate_lineage(targets, "target")
        _validate_mapping_identity(
            reference_volumes if reference_volumes is not None else reference_volume,
            references,
            "reference volumes",
        )
        _validate_mapping_identity(
            target_volumes if target_volumes is not None else target_volume,
            targets,
            "target volumes",
        )
        _validate_mapping_identity(
            reference_surface_areas_angstrom2
            if reference_surface_areas_angstrom2 is not None
            else reference_surface_area_angstrom2,
            references,
            "reference surfaces",
        )
        _validate_mapping_identity(
            target_surface_areas_angstrom2
            if target_surface_areas_angstrom2 is not None
            else target_surface_area_angstrom2,
            targets,
            "target surfaces",
        )
        _validate_bound_volume(
            reference_volumes if reference_volumes is not None else reference_volume,
            references,
            "reference volumes",
        )
        _validate_bound_volume(
            target_volumes if target_volumes is not None else target_volume,
            targets,
            "target volumes",
        )
        _validate_bound_surface(
            reference_surface_areas_angstrom2
            if reference_surface_areas_angstrom2 is not None
            else reference_surface_area_angstrom2,
            references,
            "reference surfaces",
        )
        _validate_bound_surface(
            target_surface_areas_angstrom2
            if target_surface_areas_angstrom2 is not None
            else target_surface_area_angstrom2,
            targets,
            "target surfaces",
        )
        inherited_state = next(
            (
                state
                for state in (reference_state, target_state)
                if state in {Availability.INVALID_INPUT, Availability.NUMERICAL_FAILURE}
            ),
            None,
        )
        if inherited_state is not None:
            matching = PocketMatchingResult(
                inherited_state,
                diagnostics=reference_diagnostics + target_diagnostics,
                settings=run_settings,
                transform=transform or StructuralTransform(),
            )
        else:
            matching = match_pocket_candidates(
                references,
                targets,
                correspondences,
                transform=transform,
                settings=run_settings,
            )
        comparisons: list[PocketComparison] = []
        diagnostics = list(reference_diagnostics + target_diagnostics + matching.diagnostics)
        reference_volume_input: Mapping[str, VolumeInput] | VolumeInput | None = (
            reference_volumes if reference_volumes is not None else reference_volume
        )
        target_volume_input: Mapping[str, VolumeInput] | VolumeInput | None = (
            target_volumes if target_volumes is not None else target_volume
        )
        for match in matching.matches:
            reference = match.reference_candidate
            target = match.target_candidate
            comparison = compare_pocket_match(
                match,
                correspondences=correspondences,
                reference_volume=_volume_for(reference_volume_input, reference),
                target_volume=_volume_for(target_volume_input, target),
                reference_surface_area_angstrom2=(
                    _surface_for(reference_surface_area_angstrom2, reference, "reference_surface_area_angstrom2")
                    if reference_surface_area_angstrom2 is not None
                    else _surface_for(reference_surface_areas_angstrom2, reference, "reference_surface_areas_angstrom2")
                ),
                target_surface_area_angstrom2=(
                    _surface_for(target_surface_area_angstrom2, target, "target_surface_area_angstrom2")
                    if target_surface_area_angstrom2 is not None
                    else _surface_for(target_surface_areas_angstrom2, target, "target_surface_areas_angstrom2")
                ),
                mutations=mutations,
                interaction_differences=interaction_differences,
                local_displacements=local_displacements,  # type: ignore[arg-type]
                qc=qc,
                reference_qc=reference_qc,
                target_qc=target_qc,
                qc_diagnostics=qc_diagnostics,
            )
            comparisons.append(comparison)
            diagnostics.extend(comparison.diagnostics)
        states = [item.availability for item in comparisons]
        if matching.availability is Availability.NUMERICAL_FAILURE:
            availability = Availability.NUMERICAL_FAILURE
        elif states and Availability.NUMERICAL_FAILURE in states:
            availability = Availability.NUMERICAL_FAILURE
        elif states and Availability.INVALID_INPUT in states:
            availability = Availability.INVALID_INPUT
        elif states and all(item is Availability.NOT_APPLICABLE for item in states):
            availability = Availability.NOT_APPLICABLE
        elif not comparisons:
            availability = matching.availability
        else:
            availability = Availability.AVAILABLE
        concordance = _build_concordance(matching, tuple(comparisons))
        return PocketComparisonReport(
            availability,
            matching,
            tuple(comparisons),
            tuple(_unique_diagnostics(diagnostics)),
            concordance,
        )

    analyze = compare
    compare_candidates = compare


def _unique_diagnostics(values: Sequence[Diagnostic]) -> tuple[Diagnostic, ...]:
    output: list[Diagnostic] = []
    seen: set[tuple[object, ...]] = set()
    for item in values:
        key = (item.code, item.severity, item.message, item.source_id, item.atom_id, item.residue_id, item.remediation)
        if key not in seen:
            output.append(item)
            seen.add(key)
    return tuple(output)


def _channel_state(
    values: Sequence[Availability], *, empty: Availability = Availability.NOT_APPLICABLE
) -> Availability:
    if not values:
        return empty
    for state in (Availability.NUMERICAL_FAILURE, Availability.INVALID_INPUT, Availability.DEPENDENCY_UNAVAILABLE):
        if state in values:
            return state
    if Availability.AVAILABLE in values:
        return Availability.AVAILABLE
    if Availability.NOT_DETECTED in values:
        return Availability.NOT_DETECTED
    return Availability.NOT_APPLICABLE


def _build_concordance(
    matching: PocketMatchingResult,
    comparisons: Sequence[PocketComparison],
) -> PocketEvidenceSections:
    surface_measurements = tuple(
        {
            "reference": {
                "area_angstrom2": item.reference_surface_area_angstrom2,
                "method": item.reference_surface_method,
                "provenance": item.reference_surface_provenance,
                "units": item.reference_surface_units,
            },
            "target": {
                "area_angstrom2": item.target_surface_area_angstrom2,
                "method": item.target_surface_method,
                "provenance": item.target_surface_provenance,
                "units": item.target_surface_units,
            },
        }
        for item in comparisons
        if item.reference_surface_area_angstrom2 is not None or item.target_surface_area_angstrom2 is not None
    )
    volumes = tuple(item.volume for item in comparisons if item.volume is not None)
    volume_states = tuple(item.availability for item in volumes)
    qc_states = tuple(item.qc_availability for item in comparisons if item.qc_availability is not None)
    interactions = tuple(item.interaction_changes for item in comparisons)
    interaction_state = (
        Availability.AVAILABLE
        if any(interactions)
        else Availability.NOT_DETECTED
        if comparisons
        else Availability.NOT_APPLICABLE
    )
    return PocketEvidenceSections(
        geometry=PocketEvidenceChannel(
            matching.availability,
            matching,
            units={"surface_area": "angstrom^2"} if surface_measurements else {},
            diagnostics=matching.diagnostics,
            provenance={"surface_measurements": surface_measurements} if surface_measurements else None,
        ),
        ligand_support=PocketEvidenceChannel(Availability.NOT_APPLICABLE),
        parameter_persistence=PocketEvidenceChannel(Availability.NOT_APPLICABLE),
        volume_sensitivity=PocketEvidenceChannel(
            _channel_state(volume_states),
            volumes or None,
            units={
                "volume": "angstrom^3",
                "sensitivity": "angstrom^3",
                "relative_sensitivity": "fraction",
            }
            if volumes
            else {},
            diagnostics=tuple(diagnostic for item in volumes for diagnostic in item.diagnostics),
            provenance=tuple(provenance for item in volumes for provenance in item.provenance) or None,
        ),
        match_ambiguity=PocketEvidenceChannel(
            Availability.AVAILABLE
            if any(item.state.value == "ambiguous" for item in matching.matches)
            else Availability.NOT_DETECTED
            if matching.matches
            else Availability.NOT_APPLICABLE,
            matching.matches or None,
            diagnostics=tuple(item for match in matching.matches for item in match.diagnostics),
        ),
        qc=PocketEvidenceChannel(
            _channel_state(qc_states),
            tuple(item.qc_diagnostics for item in comparisons) or None,
            diagnostics=tuple(item for comparison in comparisons for item in comparison.qc_diagnostics),
        ),
        interactions=PocketEvidenceChannel(
            interaction_state,
            interactions or None,
            diagnostics=tuple(
                item
                for comparison in comparisons
                for item in comparison.diagnostics
                if item.code.startswith("interaction")
            ),
        ),
    )


__all__ = ["PocketComparisonReport", "PocketComparisonService"]
