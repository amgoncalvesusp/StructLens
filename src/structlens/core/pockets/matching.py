"""Deterministic pocket-candidate matching in the reference coordinate frame.

Candidate identity is never inferred from source numbering alone.  Lining
residues on the target are projected through the authoritative correspondence
map, and a bounded bipartite assignment combines that overlap with the
transformed centroid distance.  Alternatives are retained as ``ambiguous``
matches instead of being hidden by a nearest-centroid tie break.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from enum import Enum
from types import MappingProxyType
from typing import Any, TypeAlias, cast

from scipy.optimize import linear_sum_assignment  # type: ignore[import-untyped]

from structlens.core.evidence import Availability, Diagnostic, DiagnosticSeverity
from structlens.core.models import ResidueCorrespondence, ResidueId, StructuralTransform
from structlens.core.provenance import freeze_bounded_string_map

from .models import PocketCandidate


class PocketMatchState(str, Enum):
    """Outcome for one candidate pair or one unmatched candidate."""

    MATCHED = "matched"
    AMBIGUOUS = "ambiguous"
    UNMATCHED_REFERENCE = "unmatched_reference"
    UNMATCHED_TARGET = "unmatched_target"


# These aliases make the contract convenient for callers that use "status" or
# "kind" terminology while retaining one canonical enum.
PocketMatchStatus = PocketMatchState
PocketMatchKind = PocketMatchState

MAX_MATCHING_ITEMS = 100_000


def _finite(value: object, name: str, *, non_negative: bool = False) -> float:
    try:
        numeric = float(cast(Any, value))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(numeric) or (non_negative and numeric < 0.0):
        suffix = " and non-negative" if non_negative else ""
        raise ValueError(f"{name} must be finite{suffix}")
    return numeric


@dataclass(frozen=True, slots=True, init=False)
class PocketMatchingSettings:
    """Bounded, explicit scoring policy for candidate assignment."""

    minimum_lining_jaccard: float = 0.1
    maximum_centroid_distance_angstrom: float = 8.0
    ambiguity_margin: float = 0.05
    lining_weight: float = 0.7
    centroid_weight: float = 0.3
    maximum_candidate_count: int = 200

    def __init__(
        self,
        minimum_lining_jaccard: float = 0.1,
        maximum_centroid_distance_angstrom: float = 8.0,
        ambiguity_margin: float = 0.05,
        lining_weight: float = 0.7,
        centroid_weight: float = 0.3,
        maximum_candidate_count: int = 200,
        *,
        min_lining_jaccard: float | None = None,
        max_centroid_distance_angstrom: float | None = None,
        lining_overlap_weight: float | None = None,
        centroid_distance_weight: float | None = None,
    ) -> None:
        if min_lining_jaccard is not None:
            minimum_lining_jaccard = min_lining_jaccard
        if max_centroid_distance_angstrom is not None:
            maximum_centroid_distance_angstrom = max_centroid_distance_angstrom
        if lining_overlap_weight is not None:
            lining_weight = lining_overlap_weight
        if centroid_distance_weight is not None:
            centroid_weight = centroid_distance_weight
        object.__setattr__(self, "minimum_lining_jaccard", minimum_lining_jaccard)
        object.__setattr__(self, "maximum_centroid_distance_angstrom", maximum_centroid_distance_angstrom)
        object.__setattr__(self, "ambiguity_margin", ambiguity_margin)
        object.__setattr__(self, "lining_weight", lining_weight)
        object.__setattr__(self, "centroid_weight", centroid_weight)
        object.__setattr__(self, "maximum_candidate_count", maximum_candidate_count)
        self.__post_init__()

    def __post_init__(self) -> None:
        minimum = _finite(self.minimum_lining_jaccard, "minimum_lining_jaccard", non_negative=True)
        if minimum > 1.0:
            raise ValueError("minimum_lining_jaccard must be between 0 and 1")
        object.__setattr__(self, "minimum_lining_jaccard", minimum)
        object.__setattr__(
            self,
            "maximum_centroid_distance_angstrom",
            _finite(self.maximum_centroid_distance_angstrom, "maximum_centroid_distance_angstrom", non_negative=True),
        )
        object.__setattr__(
            self, "ambiguity_margin", _finite(self.ambiguity_margin, "ambiguity_margin", non_negative=True)
        )
        lining = _finite(self.lining_weight, "lining_weight", non_negative=True)
        centroid = _finite(self.centroid_weight, "centroid_weight", non_negative=True)
        if lining + centroid <= 0.0:
            raise ValueError("lining_weight and centroid_weight cannot both be zero")
        object.__setattr__(self, "lining_weight", lining)
        object.__setattr__(self, "centroid_weight", centroid)
        if isinstance(self.maximum_candidate_count, bool) or not isinstance(self.maximum_candidate_count, int):
            raise ValueError("maximum_candidate_count must be a positive integer")
        if self.maximum_candidate_count <= 0:
            raise ValueError("maximum_candidate_count must be a positive integer")
        if self.maximum_candidate_count > 200:
            raise ValueError("maximum_candidate_count must not exceed 200")

    @property
    def min_lining_jaccard(self) -> float:
        return self.minimum_lining_jaccard

    @property
    def max_centroid_distance_angstrom(self) -> float:
        return self.maximum_centroid_distance_angstrom

    @property
    def centroid_distance_weight(self) -> float:
        return self.centroid_weight

    @property
    def lining_overlap_weight(self) -> float:
        return self.lining_weight

    def to_json(self) -> dict[str, object]:
        return {
            "minimum_lining_jaccard": self.minimum_lining_jaccard,
            "maximum_centroid_distance_angstrom": self.maximum_centroid_distance_angstrom,
            "ambiguity_margin": self.ambiguity_margin,
            "lining_weight": self.lining_weight,
            "centroid_weight": self.centroid_weight,
            "maximum_candidate_count": self.maximum_candidate_count,
        }


def _residue_key(residue: ResidueId) -> tuple[str, str, str, str, str, str, str]:
    return (
        residue.structure_id,
        residue.model_id,
        residue.chain_id,
        residue.auth_seq_id,
        residue.insertion_code or "",
        residue.residue_name,
        repr(residue),
    )


@dataclass(frozen=True, slots=True)
class PocketMatch:
    """One transparent candidate assignment and its component measures."""

    reference_candidate: PocketCandidate | None
    target_candidate: PocketCandidate | None
    state: PocketMatchState
    lining_jaccard: float | None = None
    centroid_distance_angstrom: float | None = None
    score: float | None = None
    alternative_score: float | None = None
    ambiguity_margin: float | None = None
    diagnostics: tuple[Diagnostic, ...] = field(default_factory=tuple)
    units: Mapping[str, str] = field(
        default_factory=lambda: MappingProxyType(
            {
                "lining_jaccard": "fraction",
                "centroid_distance_angstrom": "angstrom",
                "score": "fraction",
            }
        )
    )
    assignment_score: float | None = None
    alternative_assignment_score: float | None = None

    def __post_init__(self) -> None:
        if self.reference_candidate is not None and not isinstance(self.reference_candidate, PocketCandidate):
            raise TypeError("reference_candidate must be a PocketCandidate or None")
        if self.target_candidate is not None and not isinstance(self.target_candidate, PocketCandidate):
            raise TypeError("target_candidate must be a PocketCandidate or None")
        state = self.state
        if not isinstance(state, PocketMatchState):
            try:
                state = PocketMatchState(state)
            except (TypeError, ValueError) as exc:
                raise ValueError(f"unknown pocket match state: {self.state!r}") from exc
            object.__setattr__(self, "state", state)
        if state is PocketMatchState.UNMATCHED_REFERENCE and self.reference_candidate is None:
            raise ValueError("unmatched_reference requires reference_candidate")
        if state is PocketMatchState.UNMATCHED_TARGET and self.target_candidate is None:
            raise ValueError("unmatched_target requires target_candidate")
        if state in {PocketMatchState.MATCHED, PocketMatchState.AMBIGUOUS} and (
            self.reference_candidate is None or self.target_candidate is None
        ):
            raise ValueError("matched states require both candidates")
        for name in ("lining_jaccard", "score"):
            value = getattr(self, name)
            if value is not None:
                value = _finite(value, name, non_negative=True)
                if name == "lining_jaccard" and value > 1.0:
                    raise ValueError("lining_jaccard must be between 0 and 1")
                object.__setattr__(self, name, value)
        for name in (
            "centroid_distance_angstrom",
            "alternative_score",
            "ambiguity_margin",
            "assignment_score",
            "alternative_assignment_score",
        ):
            value = getattr(self, name)
            if value is not None:
                object.__setattr__(self, name, _finite(value, name, non_negative=True))
        if not isinstance(self.diagnostics, Sequence) or isinstance(self.diagnostics, (str, bytes)):
            raise TypeError("diagnostics must be a bounded Sequence")
        if len(self.diagnostics) > MAX_MATCHING_ITEMS:
            raise ValueError("diagnostics exceed the bounded matching limit")
        diagnostics = tuple(self.diagnostics)
        if any(not isinstance(item, Diagnostic) for item in diagnostics):
            raise TypeError("diagnostics must contain Diagnostic values")
        object.__setattr__(self, "diagnostics", diagnostics)
        object.__setattr__(self, "units", freeze_bounded_string_map(self.units, field_name="units"))

    @property
    def reference(self) -> PocketCandidate | None:
        return self.reference_candidate

    @property
    def target(self) -> PocketCandidate | None:
        return self.target_candidate

    @property
    def status(self) -> PocketMatchState:
        return self.state

    @property
    def kind(self) -> PocketMatchState:
        return self.state

    @property
    def lining_overlap_fraction(self) -> float | None:
        return self.lining_jaccard

    @property
    def lining_residue_jaccard(self) -> float | None:
        return self.lining_jaccard

    @property
    def transformed_centroid_distance_angstrom(self) -> float | None:
        return self.centroid_distance_angstrom

    @property
    def global_assignment_score(self) -> float | None:
        return self.assignment_score

    @property
    def alternative_global_assignment_score(self) -> float | None:
        return self.alternative_assignment_score

    def to_json(self) -> dict[str, object]:
        return {
            "state": self.state.value,
            "status": self.state.value,
            "reference_candidate": self.reference_candidate.to_json() if self.reference_candidate else None,
            "target_candidate": self.target_candidate.to_json() if self.target_candidate else None,
            "lining_jaccard": self.lining_jaccard,
            "lining_overlap_fraction": self.lining_overlap_fraction,
            "centroid_distance_angstrom": self.centroid_distance_angstrom,
            "score": self.score,
            "alternative_score": self.alternative_score,
            "ambiguity_margin": self.ambiguity_margin,
            "assignment_score": self.assignment_score,
            "alternative_assignment_score": self.alternative_assignment_score,
            "units": dict(self.units),
            "diagnostics": [item.to_json() for item in self.diagnostics],
        }


@dataclass(frozen=True, slots=True)
class PocketMatchingResult:
    """Immutable collection of matched and unmatched candidate states."""

    availability: Availability
    matches: tuple[PocketMatch, ...] = field(default_factory=tuple)
    diagnostics: tuple[Diagnostic, ...] = field(default_factory=tuple)
    settings: PocketMatchingSettings = field(default_factory=PocketMatchingSettings)
    transform: StructuralTransform = field(default_factory=StructuralTransform)
    optimal_score: float | None = None
    second_best_score: float | None = None

    def __post_init__(self) -> None:
        availability = self.availability
        if not isinstance(availability, Availability):
            availability = Availability(availability)
            object.__setattr__(self, "availability", availability)
        if not isinstance(self.matches, Sequence) or isinstance(self.matches, (str, bytes)):
            raise TypeError("matches must be a bounded Sequence")
        if len(self.matches) > MAX_MATCHING_ITEMS:
            raise ValueError("matches exceed the bounded matching limit")
        matches = tuple(self.matches)
        if any(not isinstance(item, PocketMatch) for item in matches):
            raise TypeError("matches must contain PocketMatch values")
        if not isinstance(self.diagnostics, Sequence) or isinstance(self.diagnostics, (str, bytes)):
            raise TypeError("diagnostics must be a bounded Sequence")
        if len(self.diagnostics) > MAX_MATCHING_ITEMS:
            raise ValueError("diagnostics exceed the bounded matching limit")
        diagnostics = tuple(self.diagnostics)
        if any(not isinstance(item, Diagnostic) for item in diagnostics):
            raise TypeError("diagnostics must contain Diagnostic values")
        if not isinstance(self.settings, PocketMatchingSettings):
            raise TypeError("settings must be PocketMatchingSettings")
        if not isinstance(self.transform, StructuralTransform):
            raise TypeError("transform must be StructuralTransform")
        object.__setattr__(self, "matches", matches)
        object.__setattr__(self, "diagnostics", diagnostics)
        for name in ("optimal_score", "second_best_score"):
            value = getattr(self, name)
            if value is not None:
                object.__setattr__(self, name, _finite(value, name, non_negative=True))

    @property
    def status(self) -> Availability:
        return self.availability

    @property
    def candidate_matches(self) -> tuple[PocketMatch, ...]:
        return self.matches

    @property
    def global_assignment_score(self) -> float | None:
        return self.optimal_score

    @property
    def alternative_global_assignment_score(self) -> float | None:
        return self.second_best_score

    def __iter__(self):  # type: ignore[no-untyped-def]
        return iter(self.matches)

    def __len__(self) -> int:
        return len(self.matches)

    def __getitem__(self, index: int) -> PocketMatch:
        return self.matches[index]

    def to_json(self) -> dict[str, object]:
        return {
            "availability": self.availability.value,
            "matches": [item.to_json() for item in self.matches],
            "diagnostics": [item.to_json() for item in self.diagnostics],
            "settings": self.settings.to_json(),
            "transform": {
                "rotation": [list(row) for row in self.transform.rotation],
                "translation": list(self.transform.translation),
            },
            "optimal_score": self.optimal_score,
            "second_best_score": self.second_best_score,
        }


MatchInput: TypeAlias = Sequence[ResidueCorrespondence] | Mapping[ResidueId, ResidueId]


def _correspondence_map(correspondences: MatchInput) -> dict[ResidueId, ResidueId]:
    if isinstance(correspondences, Mapping):
        if len(correspondences) > MAX_MATCHING_ITEMS:
            raise ValueError("correspondences exceed the bounded matching limit")
        pairs = tuple(correspondences.items())
        if any(not isinstance(left, ResidueId) or not isinstance(right, ResidueId) for left, right in pairs):
            raise TypeError("correspondence mappings must contain ResidueId pairs")
        references = tuple(left for left, _ in pairs)
        targets = tuple(right for _, right in pairs)
        if len(set(references)) != len(references):
            raise ValueError("authoritative correspondence contains duplicate reference residues")
        if len(set(targets)) != len(targets):
            raise ValueError("authoritative correspondence contains duplicate target residues")
        return dict(pairs)
    if not isinstance(correspondences, Sequence) or isinstance(correspondences, (str, bytes)):
        raise TypeError("correspondences must be a Sequence or Mapping")
    if len(correspondences) > MAX_MATCHING_ITEMS:
        raise ValueError("correspondences exceed the bounded matching limit")
    values = tuple(correspondences)
    if any(not isinstance(item, ResidueCorrespondence) for item in values):
        raise TypeError("correspondences must contain ResidueCorrespondence values")
    output: dict[ResidueId, ResidueId] = {}
    reverse: dict[ResidueId, ResidueId] = {}
    references_seen: set[ResidueId] = set()
    targets_seen: set[ResidueId] = set()
    for item in values:
        if item.reference is not None:
            if item.reference in references_seen:
                raise ValueError("authoritative correspondence contains duplicate reference residues")
            references_seen.add(item.reference)
        if item.target is not None:
            if item.target in targets_seen:
                raise ValueError("authoritative correspondence contains duplicate target residues")
            targets_seen.add(item.target)
        if item.reference is None or item.target is None:
            continue
        if item.reference in output:
            raise ValueError("authoritative correspondence contains duplicate reference residues")
        if item.target in reverse:
            raise ValueError("authoritative correspondence contains duplicate target residues")
        output[item.reference] = item.target
        reverse[item.target] = item.reference
    return output


def _transform(point: tuple[float, float, float], transform: StructuralTransform) -> tuple[float, float, float]:
    # StructuralTransform uses row-vector semantics: target @ rotation + t.
    transformed = tuple(
        sum(point[row] * transform.rotation[row][column] for row in range(3)) + transform.translation[column]
        for column in range(3)
    )
    return cast(tuple[float, float, float], transformed)


def _distance(first: tuple[float, float, float], second: tuple[float, float, float]) -> float:
    return math.sqrt(sum((left - right) ** 2 for left, right in zip(first, second, strict=True)))


def _edge(
    reference: PocketCandidate,
    target: PocketCandidate,
    correspondence: Mapping[ResidueId, ResidueId],
    transform: StructuralTransform,
    settings: PocketMatchingSettings,
) -> tuple[float, float, float] | None:
    reference_lining = set(reference.lining_residues)
    inverse = {target_id: reference_id for reference_id, target_id in correspondence.items()}
    target_lining = {inverse[item] for item in target.lining_residues if item in inverse}
    union = reference_lining | target_lining
    jaccard = len(reference_lining & target_lining) / len(union) if union else 0.0
    transformed_target_centroid = _transform(target.centroid_xyz, transform)
    distance = _distance(reference.centroid_xyz, transformed_target_centroid)
    if jaccard < settings.minimum_lining_jaccard or distance > settings.maximum_centroid_distance_angstrom:
        return None
    distance_similarity = (
        max(0.0, 1.0 - distance / settings.maximum_centroid_distance_angstrom)
        if settings.maximum_centroid_distance_angstrom > 0.0
        else 1.0
        if distance == 0.0
        else 0.0
    )
    total_weight = settings.lining_weight + settings.centroid_weight
    score = (settings.lining_weight * jaccard + settings.centroid_weight * distance_similarity) / total_weight
    return (jaccard, distance, score)


def _optimal_assignment(
    edges: Sequence[Mapping[int, tuple[float, float, float]]],
    targets: Sequence[PocketCandidate],
    *,
    forbidden: tuple[int, int] | None = None,
) -> tuple[tuple[int | None, ...], float]:
    """Solve a bounded assignment with explicit unmatched dummy rows/columns."""

    row_count = len(edges)
    target_count = len(targets)
    size = row_count + target_count
    if size == 0:
        return (), 0.0
    # A large finite penalty keeps ineligible edges out while allowing every
    # real row/column to use an explicit unmatched dummy counterpart.
    invalid_cost = 10_000.0
    costs = [[0.0 for _ in range(size)] for _ in range(size)]
    for row in range(row_count):
        for column in range(target_count):
            edge = edges[row].get(column)
            if edge is None or forbidden == (row, column):
                costs[row][column] = invalid_cost
            else:
                costs[row][column] = -edge[2]
        for column in range(target_count, size):
            costs[row][column] = 0.0
    for row in range(row_count, size):
        for column in range(size):
            costs[row][column] = 0.0
    row_indices, column_indices = linear_sum_assignment(costs)
    assignment: list[int | None] = [None] * row_count
    total = 0.0
    for row, column in zip(row_indices, column_indices, strict=True):
        if row >= row_count or column >= target_count:
            continue
        edge = edges[row].get(int(column))
        if edge is not None and forbidden != (row, int(column)):
            assignment[int(row)] = int(column)
            total += edge[2]
    return tuple(assignment), total


def match_pocket_candidates(
    reference_candidates: Sequence[PocketCandidate],
    target_candidates: Sequence[PocketCandidate],
    correspondences: MatchInput = (),
    *,
    transform: StructuralTransform | None = None,
    settings: PocketMatchingSettings | None = None,
) -> PocketMatchingResult:
    """Match candidate collections through one authoritative residue map."""

    if not isinstance(reference_candidates, Sequence) or isinstance(reference_candidates, (str, bytes)):
        raise TypeError("reference_candidates must be a bounded Sequence")
    if not isinstance(target_candidates, Sequence) or isinstance(target_candidates, (str, bytes)):
        raise TypeError("target_candidates must be a bounded Sequence")
    run_settings = settings if settings is not None else PocketMatchingSettings()
    if not isinstance(run_settings, PocketMatchingSettings):
        raise TypeError("settings must be PocketMatchingSettings or None")
    correspondence = _correspondence_map(correspondences)
    run_transform = transform if transform is not None else StructuralTransform()
    if not isinstance(run_transform, StructuralTransform):
        raise TypeError("transform must be StructuralTransform or None")
    if (
        len(reference_candidates) > run_settings.maximum_candidate_count
        or len(target_candidates) > run_settings.maximum_candidate_count
    ):
        diagnostic = Diagnostic(
            "pocket.match.resource_limit",
            DiagnosticSeverity.ERROR,
            "Pocket candidate matching exceeded the configured candidate limit.",
            remediation="Reduce the displayed candidate set or raise maximum_candidate_count deliberately.",
        )
        return PocketMatchingResult(
            Availability.NUMERICAL_FAILURE,
            diagnostics=(diagnostic,),
            settings=run_settings,
            transform=run_transform,
        )
    reference = tuple(reference_candidates)
    target = tuple(target_candidates)
    if any(not isinstance(item, PocketCandidate) for item in reference + target):
        raise TypeError("candidate collections must contain PocketCandidate values")
    if len({item.candidate_id for item in reference}) != len(reference) or len(
        {item.candidate_id for item in target}
    ) != len(target):
        raise ValueError("candidate_id collision makes pocket ordering ambiguous")
    ordered_reference = tuple(sorted(reference, key=lambda item: item.candidate_id))
    ordered_target = tuple(sorted(target, key=lambda item: item.candidate_id))
    edges = (
        tuple(
            {
                target_index: measures
                for target_index, target_candidate in enumerate(ordered_target)
                if (
                    measures := _edge(
                        reference_candidate, target_candidate, correspondence, run_transform, run_settings
                    )
                )
                is not None
            }
            for reference_candidate in ordered_reference
        )
        if correspondence
        else tuple({} for _ in ordered_reference)
    )
    chosen, optimal_score = _optimal_assignment(edges, ordered_target)
    alternatives_by_row: dict[int, tuple[tuple[int | None, ...], float]] = {}
    for row, target_index in enumerate(chosen):
        if target_index is None:
            continue
        alternative_assignment, alternative_total = _optimal_assignment(
            edges, ordered_target, forbidden=(row, target_index)
        )
        # The explicit unmatched dummy assignment is a valid global
        # alternative and must remain in the same total-score scale.
        alternatives_by_row[row] = (alternative_assignment, alternative_total)
    second_best_score = max((item[1] for item in alternatives_by_row.values()), default=None)
    matches: list[PocketMatch] = []
    used_targets: set[int] = set()
    diagnostics: list[Diagnostic] = []
    for row, reference_candidate in enumerate(ordered_reference):
        target_index = chosen[row] if row < len(chosen) else None
        measures = edges[row].get(target_index) if target_index is not None else None
        if target_index is None or measures is None:
            unmatched_diagnostic = Diagnostic(
                "pocket.match.unmatched_reference",
                DiagnosticSeverity.WARNING,
                "No target pocket candidate passed the declared matching thresholds for this reference candidate.",
                source_id=reference_candidate.candidate_id,
            )
            diagnostics.append(unmatched_diagnostic)
            matches.append(
                PocketMatch(
                    reference_candidate, None, PocketMatchState.UNMATCHED_REFERENCE, diagnostics=(unmatched_diagnostic,)
                )
            )
            continue
        used_targets.add(target_index)
        jaccard, distance, score = measures
        alternative = alternatives_by_row.get(row)
        alternative_assignment_score = alternative[1] if alternative is not None else None
        alternative_index = alternative[0][row] if alternative is not None else None
        alternative_score = (
            edges[row][alternative_index][2]
            if alternative_index is not None and alternative_index in edges[row]
            else None
        )
        is_ambiguous = (
            alternative_assignment_score is not None
            and optimal_score - alternative_assignment_score <= run_settings.ambiguity_margin + 1e-12
        )
        state = PocketMatchState.AMBIGUOUS if is_ambiguous else PocketMatchState.MATCHED
        if is_ambiguous:
            diagnostics.append(
                Diagnostic(
                    "pocket.match.ambiguous",
                    DiagnosticSeverity.WARNING,
                    "More than one candidate assignment remains within the declared score margin.",
                    source_id=reference_candidate.candidate_id,
                )
            )
        match_diagnostics: tuple[Diagnostic, ...] = ()
        if is_ambiguous:
            match_diagnostics = (diagnostics[-1],)
        matches.append(
            PocketMatch(
                reference_candidate,
                ordered_target[target_index],
                state,
                jaccard,
                distance,
                score,
                alternative_score,
                run_settings.ambiguity_margin if is_ambiguous else None,
                diagnostics=match_diagnostics,
                assignment_score=optimal_score,
                alternative_assignment_score=alternative_assignment_score,
            )
        )
    for target_index, target_candidate in enumerate(ordered_target):
        if target_index not in used_targets:
            unmatched_diagnostic = Diagnostic(
                "pocket.match.unmatched_target",
                DiagnosticSeverity.WARNING,
                "No reference pocket candidate was assigned to this target candidate.",
                source_id=target_candidate.candidate_id,
            )
            diagnostics.append(unmatched_diagnostic)
            matches.append(
                PocketMatch(
                    None, target_candidate, PocketMatchState.UNMATCHED_TARGET, diagnostics=(unmatched_diagnostic,)
                )
            )
    state_order = {
        PocketMatchState.MATCHED: 0,
        PocketMatchState.AMBIGUOUS: 0,
        PocketMatchState.UNMATCHED_REFERENCE: 1,
        PocketMatchState.UNMATCHED_TARGET: 2,
    }
    matches.sort(
        key=lambda item: (
            state_order[item.state],
            item.reference_candidate.candidate_id if item.reference_candidate else "",
            item.target_candidate.candidate_id if item.target_candidate else "",
        )
    )
    if not matches:
        availability = Availability.NOT_APPLICABLE
    elif any(item.state is PocketMatchState.AMBIGUOUS for item in matches):
        availability = Availability.AVAILABLE
    elif any(
        item.state
        in {PocketMatchState.MATCHED, PocketMatchState.UNMATCHED_REFERENCE, PocketMatchState.UNMATCHED_TARGET}
        for item in matches
    ):
        availability = Availability.AVAILABLE
    else:
        availability = Availability.NOT_DETECTED
    return PocketMatchingResult(
        availability,
        tuple(matches),
        tuple(diagnostics),
        run_settings,
        run_transform,
        optimal_score,
        second_best_score,
    )


# Friendly aliases used by application integrations and earlier design notes.
match_pockets = match_pocket_candidates
match_candidates = match_pocket_candidates


__all__ = [
    "PocketMatch",
    "PocketMatchKind",
    "PocketMatchState",
    "PocketMatchStatus",
    "PocketMatchingResult",
    "PocketMatchingSettings",
    "match_candidates",
    "match_pocket_candidates",
    "match_pockets",
]
