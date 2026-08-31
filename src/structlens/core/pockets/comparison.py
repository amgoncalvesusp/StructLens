"""Descriptive comparison of two matched pocket candidates.

This module performs bounded orchestration only.  Immutable result contracts
and canonical serializers live in :mod:`comparison_models`.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any, cast

from structlens.core.difference_maps import ResidueDisplacementVector
from structlens.core.evidence import Availability, Diagnostic, DiagnosticSeverity
from structlens.core.interactions import InteractionDifference
from structlens.core.models import MutationEvent, ResidueCorrespondence, ResidueId
from structlens.core.provenance import MethodProvenance

from .comparison_models import PocketComparison, PocketSurfaceMeasurement, residue_key
from .matching import PocketMatch, PocketMatchState, match_pocket_candidates
from .models import PocketCandidate
from .volume import compare_pocket_volumes
from .volume_models import PocketVolumeComparison, PocketVolumeResult, PocketVolumeSensitivity

MAX_COMPARISON_ITEMS = 100_000


def _finite(value: object, name: str, *, non_negative: bool = False) -> float:
    try:
        numeric = float(cast(Any, value))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(numeric) or (non_negative and numeric < 0.0):
        suffix = " and non-negative" if non_negative else ""
        raise ValueError(f"{name} must be finite{suffix}")
    return numeric


def _correspondence_maps(
    correspondences: Sequence[ResidueCorrespondence] | Mapping[ResidueId, ResidueId],
) -> tuple[dict[ResidueId, ResidueId], dict[ResidueId, ResidueId]]:
    if isinstance(correspondences, Mapping):
        if len(correspondences) > MAX_COMPARISON_ITEMS:
            raise ValueError("correspondences exceed the bounded comparison limit")
        pairs = tuple(correspondences.items())
        if any(not isinstance(left, ResidueId) or not isinstance(right, ResidueId) for left, right in pairs):
            raise TypeError("correspondence mappings must contain ResidueId pairs")
    else:
        if not isinstance(correspondences, Sequence) or isinstance(correspondences, (str, bytes)):
            raise TypeError("correspondences must be a Sequence or Mapping")
        if len(correspondences) > MAX_COMPARISON_ITEMS:
            raise ValueError("correspondences exceed the bounded comparison limit")
        rows = tuple(correspondences)
        if any(not isinstance(item, ResidueCorrespondence) for item in rows):
            raise TypeError("correspondences must contain ResidueCorrespondence values")
        pairs = tuple(
            (item.reference, item.target) for item in rows if item.reference is not None and item.target is not None
        )
    references = tuple(left for left, _ in pairs)
    targets = tuple(right for _, right in pairs)
    if len(set(references)) != len(references):
        raise ValueError("authoritative correspondence contains duplicate reference residues")
    if len(set(targets)) != len(targets):
        raise ValueError("authoritative correspondence contains duplicate target residues")
    forward = dict(pairs)
    reverse = {target: reference for reference, target in forward.items()}
    return forward, reverse


def _volume_from_input(value: PocketVolumeResult | PocketVolumeComparison | float | None, side: str) -> float | None:
    if value is None:
        return None
    if isinstance(value, PocketVolumeResult):
        return value.fine_volume_angstrom3
    if isinstance(value, PocketVolumeComparison):
        return value.reference_volume_angstrom3 if side == "reference" else value.target_volume_angstrom3
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return _finite(value, f"{side}_volume_angstrom3", non_negative=True)
    raise TypeError("volume inputs must be PocketVolumeResult, PocketVolumeComparison, or None")


def _volume_delta(
    reference: float | None,
    target: float | None,
    signature: tuple[object, ...] = (),
    *,
    reference_sensitivity: PocketVolumeSensitivity | None = None,
    target_sensitivity: PocketVolumeSensitivity | None = None,
    provenance: tuple[MethodProvenance, ...] = (),
) -> PocketVolumeComparison:
    if reference is None or target is None:
        return PocketVolumeComparison(
            Availability.NOT_APPLICABLE,
            reference_volume_angstrom3=reference,
            target_volume_angstrom3=target,
            compatibility_signature=signature,
            reference_sensitivity=reference_sensitivity,
            target_sensitivity=target_sensitivity,
            provenance=provenance,
        )
    delta = target - reference
    return PocketVolumeComparison(
        Availability.AVAILABLE,
        delta_angstrom3=delta,
        relative_delta_fraction=delta / reference if reference > 0.0 else None,
        reference_volume_angstrom3=reference,
        target_volume_angstrom3=target,
        compatibility_signature=signature,
        reference_sensitivity=reference_sensitivity,
        target_sensitivity=target_sensitivity,
        provenance=provenance,
    )


def _surface_delta(
    reference: PocketSurfaceMeasurement | float | None, target: PocketSurfaceMeasurement | float | None
) -> tuple[object, ...]:
    ref_measurement = reference if isinstance(reference, PocketSurfaceMeasurement) else None
    target_measurement = target if isinstance(target, PocketSurfaceMeasurement) else None
    ref_value = ref_measurement.surface_area_angstrom2 if ref_measurement is not None else reference
    target_value = target_measurement.surface_area_angstrom2 if target_measurement is not None else target
    ref = _finite(ref_value, "reference_surface_area_angstrom2", non_negative=True) if ref_value is not None else None
    tar = (
        _finite(target_value, "target_surface_area_angstrom2", non_negative=True) if target_value is not None else None
    )
    metadata: tuple[object, ...] = (
        ref_measurement.method if ref_measurement is not None else None,
        target_measurement.method if target_measurement is not None else None,
        ref_measurement.provenance if ref_measurement is not None else None,
        target_measurement.provenance if target_measurement is not None else None,
        ref_measurement.units if ref_measurement is not None else None,
        target_measurement.units if target_measurement is not None else None,
    )
    if ref is None or tar is None:
        diagnostics = (
            (
                Diagnostic(
                    "pocket.surface.incomplete",
                    DiagnosticSeverity.WARNING,
                    "Surface delta is unavailable because one side has no validated surface measurement.",
                ),
            )
            if (ref is not None or tar is not None)
            else ()
        )
        return None, None, ref, tar, diagnostics, *metadata
    if (ref_measurement is None) != (target_measurement is None):
        return (
            None,
            None,
            ref,
            tar,
            (
                Diagnostic(
                    "pocket.surface.incompatible_method",
                    DiagnosticSeverity.WARNING,
                    "Surface delta requires compatible typed measurements on both sides.",
                ),
            ),
            *metadata,
        )
    if (
        ref_measurement is not None
        and target_measurement is not None
        and (
            ref_measurement.method != target_measurement.method
            or dict(ref_measurement.units) != dict(target_measurement.units)
        )
    ):
        return (
            None,
            None,
            ref,
            tar,
            (
                Diagnostic(
                    "pocket.surface.incompatible_method",
                    DiagnosticSeverity.WARNING,
                    "Surface delta is unavailable because surface methods or units differ.",
                ),
            ),
            *metadata,
        )
    delta = tar - ref
    return delta, delta / ref if ref > 0.0 else None, ref, tar, (), *metadata


def _normalise_displacements(value: object) -> tuple[tuple[ResidueId, float], ...]:
    if value is None:
        return ()
    if isinstance(value, Mapping):
        if len(value) > MAX_COMPARISON_ITEMS:
            raise ValueError("local displacements exceed the bounded comparison limit")
        pairs = tuple(value.items())
    else:
        if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
            raise TypeError("local displacements must be a Mapping or Sequence")
        if len(value) > MAX_COMPARISON_ITEMS:
            raise ValueError("local displacements exceed the bounded comparison limit")
        pairs_list: list[tuple[object, object]] = []
        for item in value:
            if isinstance(item, ResidueDisplacementVector):
                pairs_list.append((item.reference_residue, item.magnitude_angstrom))
            elif isinstance(item, tuple) and len(item) == 2:
                pairs_list.append((item[0], item[1]))
            else:
                residue = getattr(item, "reference_residue", None)
                magnitude = getattr(item, "magnitude_angstrom", None)
                if residue is None or magnitude is None:
                    raise TypeError("local displacements must contain typed displacement vectors")
                pairs_list.append((residue, magnitude))
        pairs = tuple(pairs_list)
    output: list[tuple[ResidueId, float]] = []
    for residue, magnitude in pairs:
        if not isinstance(residue, ResidueId):
            raise TypeError("local displacement mappings must be keyed by ResidueId")
        output.append((residue, _finite(magnitude, "magnitude_angstrom", non_negative=True)))
    if len({item[0] for item in output}) != len(output):
        raise ValueError("local displacements must not contain duplicate residues")
    return tuple(sorted(output, key=lambda item: residue_key(item[0])))


def _validate_displacement_correspondence(value: object, correspondence: Mapping[ResidueId, ResidueId]) -> None:
    """Reject typed vectors whose target is not the authoritative mapped target."""

    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        return
    for item in value:
        if isinstance(item, ResidueDisplacementVector) and item.reference_residue in correspondence:
            if item.target_residue != correspondence[item.reference_residue]:
                raise ValueError("local displacement target does not match authoritative correspondence")


def _qc_values(
    qc: Availability | Sequence[Diagnostic] | object | None,
) -> tuple[Availability | None, tuple[Diagnostic, ...]]:
    if qc is None:
        return None, ()
    if isinstance(qc, Availability):
        return qc, ()
    if isinstance(qc, Sequence) and not isinstance(qc, (str, bytes)):
        if len(qc) > MAX_COMPARISON_ITEMS:
            raise ValueError("qc diagnostics exceed the bounded comparison limit")
        diagnostics = tuple(qc)
        if any(not isinstance(item, Diagnostic) for item in diagnostics):
            raise TypeError("qc diagnostics must contain Diagnostic values")
        state = (
            Availability.INVALID_INPUT
            if any(item.severity is DiagnosticSeverity.ERROR for item in diagnostics)
            else Availability.AVAILABLE
        )
        return state, diagnostics
    if not hasattr(qc, "availability"):
        raise TypeError("qc must be Availability, diagnostics Sequence, or typed QC result")
    availability = cast(Any, qc).availability
    if not isinstance(availability, Availability):
        try:
            availability = Availability(availability)
        except (TypeError, ValueError) as exc:
            raise ValueError("typed QC result must expose Availability") from exc
    diagnostics = getattr(qc, "diagnostics", ())
    if not isinstance(diagnostics, Sequence) or isinstance(diagnostics, (str, bytes)):
        raise TypeError("qc diagnostics must be a bounded Sequence")
    if len(diagnostics) > MAX_COMPARISON_ITEMS:
        raise ValueError("qc diagnostics exceed the bounded comparison limit")
    values = tuple(diagnostics)
    if any(not isinstance(item, Diagnostic) for item in values):
        raise TypeError("qc diagnostics must contain Diagnostic values")
    return availability, values


def _interaction_touches_reference_lining(
    item: InteractionDifference,
    reference_lining: set[ResidueId],
    target_lining: set[ResidueId],
    reverse: Mapping[ResidueId, ResidueId],
) -> bool:
    for record in (item.reference_record, item.target_record):
        if record is None:
            continue
        for residue in (record.residue_a, record.residue_b):
            if (
                residue in reference_lining
                or residue in target_lining
                or (residue is not None and reverse.get(residue) in reference_lining)
            ):
                return True
    return False


def _precomputed_volume_pair(
    reference: PocketVolumeComparison, target: PocketVolumeComparison
) -> PocketVolumeComparison:
    if reference.compatibility_signature != target.compatibility_signature:
        return PocketVolumeComparison(
            Availability.NOT_APPLICABLE,
            diagnostics=(
                Diagnostic(
                    "pocket.volume.incompatible_settings",
                    DiagnosticSeverity.WARNING,
                    "Pocket volume comparisons with different method settings cannot be combined.",
                ),
            ),
        )
    if reference.availability is not Availability.AVAILABLE or target.availability is not Availability.AVAILABLE:
        state = (
            Availability.NUMERICAL_FAILURE
            if Availability.NUMERICAL_FAILURE in {reference.availability, target.availability}
            else Availability.INVALID_INPUT
            if Availability.INVALID_INPUT in {reference.availability, target.availability}
            else Availability.NOT_APPLICABLE
        )
        return PocketVolumeComparison(
            state,
            reference_volume_angstrom3=reference.reference_volume_angstrom3,
            target_volume_angstrom3=target.target_volume_angstrom3,
            diagnostics=reference.diagnostics + target.diagnostics,
            compatibility_signature=reference.compatibility_signature,
            reference_sensitivity=reference.reference_sensitivity,
            target_sensitivity=target.target_sensitivity,
            provenance=reference.provenance + target.provenance,
        )
    return _volume_delta(
        reference.reference_volume_angstrom3,
        target.target_volume_angstrom3,
        reference.compatibility_signature,
        reference_sensitivity=reference.reference_sensitivity,
        target_sensitivity=target.target_sensitivity,
        provenance=reference.provenance + target.provenance,
    )


def _validate_volume_input(value: object, name: str) -> None:
    if value is None or isinstance(value, (PocketVolumeResult, PocketVolumeComparison)):
        return
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        _finite(value, name, non_negative=True)
        return
    raise TypeError("volume inputs must be typed volume results/comparisons or finite numbers")


def _validate_surface_input(value: object, name: str) -> None:
    if value is None or isinstance(value, PocketSurfaceMeasurement):
        return
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        _finite(value, name, non_negative=True)
        return
    raise TypeError("surface inputs must be PocketSurfaceMeasurement or finite numbers")


def _validate_pre_matching_inputs(correspondences: object, values: Mapping[str, object]) -> None:
    """Validate all comparison evidence before candidate matching or iteration."""

    forward, _ = _correspondence_maps(cast(Any, correspondences))
    for name, item_type in (
        ("mutations", MutationEvent),
        ("interaction_differences", InteractionDifference),
    ):
        value = values.get(name, ())
        if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
            raise TypeError(f"{name} must be a bounded Sequence")
        if len(value) > MAX_COMPARISON_ITEMS:
            raise ValueError(f"{name} exceed the bounded comparison limit")
        if any(not isinstance(item, item_type) for item in value):
            raise TypeError(f"{name} must contain {item_type.__name__} values")
    displacement = values.get("local_displacements")
    _validate_displacement_correspondence(displacement, forward)
    _normalise_displacements(displacement)
    for name in ("qc", "reference_qc", "target_qc", "qc_diagnostics"):
        if name in values:
            _qc_values(values[name])
    for name in ("reference_volume", "target_volume"):
        _validate_volume_input(values.get(name), name)
    for name in ("reference_surface_area_angstrom2", "target_surface_area_angstrom2"):
        _validate_surface_input(values.get(name), name)


def compare_pocket_match(
    match: PocketMatch,
    *,
    correspondences: Sequence[ResidueCorrespondence] | Mapping[ResidueId, ResidueId] = (),
    reference_volume: PocketVolumeResult | PocketVolumeComparison | float | None = None,
    target_volume: PocketVolumeResult | PocketVolumeComparison | float | None = None,
    reference_surface_area_angstrom2: PocketSurfaceMeasurement | float | None = None,
    target_surface_area_angstrom2: PocketSurfaceMeasurement | float | None = None,
    mutations: Sequence[MutationEvent] = (),
    interaction_differences: Sequence[InteractionDifference] = (),
    local_displacements: Mapping[ResidueId, float]
    | Sequence[ResidueDisplacementVector]
    | Sequence[tuple[ResidueId, float]]
    | None = None,
    qc: Availability | Sequence[Diagnostic] | object | None = None,
    reference_qc: Availability | Sequence[Diagnostic] | object | None = None,
    target_qc: Availability | Sequence[Diagnostic] | object | None = None,
    qc_diagnostics: Sequence[Diagnostic] | None = None,
) -> PocketComparison:
    if not isinstance(match, PocketMatch):
        raise TypeError("match must be a PocketMatch")
    _validate_pre_matching_inputs(
        correspondences,
        {
            "mutations": mutations,
            "interaction_differences": interaction_differences,
            "local_displacements": local_displacements,
            "qc": qc,
            "reference_qc": reference_qc,
            "target_qc": target_qc,
            "qc_diagnostics": qc_diagnostics,
            "reference_volume": reference_volume,
            "target_volume": target_volume,
            "reference_surface_area_angstrom2": reference_surface_area_angstrom2,
            "target_surface_area_angstrom2": target_surface_area_angstrom2,
        },
    )
    forward, reverse = _correspondence_maps(correspondences)
    reference_candidate = match.reference_candidate
    target_candidate = match.target_candidate
    reference_lining = set(reference_candidate.lining_residues) if reference_candidate is not None else set()
    target_lining = set(target_candidate.lining_residues) if target_candidate is not None else set()
    authoritative = bool(forward)
    comparable_match = match.state in {PocketMatchState.MATCHED, PocketMatchState.AMBIGUOUS}
    projected_target = {reverse[item]: item for item in target_lining if item in reverse} if authoritative else {}
    conserved = reference_lining & set(projected_target) if authoritative and comparable_match else set()
    losses = reference_lining - set(projected_target) if authoritative and comparable_match else set()
    gains = (
        {
            target_id
            for target_id in target_lining
            if target_id not in reverse or reverse[target_id] not in reference_lining
        }
        if authoritative and comparable_match
        else set()
    )
    mutation_values = tuple(mutations)
    associated = tuple(
        item
        for item in mutation_values
        if authoritative
        and comparable_match
        and (
            (item.reference is not None and item.reference in reference_lining)
            or (item.target is not None and item.target in target_lining)
        )
    )
    interaction_values = tuple(interaction_differences)
    interaction_values = tuple(
        item
        for item in interaction_values
        if authoritative
        and comparable_match
        and _interaction_touches_reference_lining(item, reference_lining, target_lining, reverse)
    )
    displacement_values = _normalise_displacements(local_displacements)
    sidechain_values: tuple[tuple[ResidueId, float], ...] = ()
    if authoritative and comparable_match and not isinstance(correspondences, Mapping):
        rows = tuple(correspondences)
        derived_ca = tuple(
            (item.reference, item.ca_displacement_angstrom)
            for item in rows
            if item.reference is not None and item.target is not None and item.ca_displacement_angstrom is not None
        )
        if not displacement_values and derived_ca:
            displacement_values = _normalise_displacements(cast(Sequence[tuple[ResidueId, float]], derived_ca))
        derived_sidechain = tuple(
            (item.reference, item.sidechain_rmsd_angstrom)
            for item in rows
            if item.reference is not None and item.target is not None and item.sidechain_rmsd_angstrom is not None
        )
        if derived_sidechain:
            sidechain_values = _normalise_displacements(cast(Sequence[tuple[ResidueId, float]], derived_sidechain))
    valid_reference_residues = set(forward)
    relevant_displacements = tuple(
        item
        for item in displacement_values
        if authoritative and comparable_match and item[0] in reference_lining and item[0] in valid_reference_residues
    )
    relevant_sidechain = tuple(
        item
        for item in sidechain_values
        if authoritative and comparable_match and item[0] in reference_lining and item[0] in valid_reference_residues
    )
    local_displacement = max((item[1] for item in relevant_displacements), default=None)
    sidechain_displacement = max((item[1] for item in relevant_sidechain), default=None)
    if (
        authoritative
        and comparable_match
        and isinstance(reference_volume, PocketVolumeResult)
        and isinstance(target_volume, PocketVolumeResult)
    ):
        volume = compare_pocket_volumes(reference_volume, target_volume)
    elif (
        authoritative
        and comparable_match
        and isinstance(reference_volume, PocketVolumeComparison)
        and isinstance(target_volume, PocketVolumeComparison)
    ):
        volume = _precomputed_volume_pair(reference_volume, target_volume)
    elif authoritative and comparable_match and (reference_volume is not None or target_volume is not None):
        result_inputs = tuple(
            item for item in (reference_volume, target_volume) if isinstance(item, PocketVolumeResult)
        )
        unavailable_results = tuple(item for item in result_inputs if item.availability is not Availability.AVAILABLE)
        if unavailable_results:
            state = (
                Availability.NUMERICAL_FAILURE
                if any(item.availability is Availability.NUMERICAL_FAILURE for item in unavailable_results)
                else Availability.INVALID_INPUT
                if any(item.availability is Availability.INVALID_INPUT for item in unavailable_results)
                else Availability.NOT_APPLICABLE
            )
            volume = PocketVolumeComparison(
                state,
                diagnostics=tuple(diagnostic for item in unavailable_results for diagnostic in item.diagnostics),
            )
        elif isinstance(reference_volume, PocketVolumeResult) or isinstance(target_volume, PocketVolumeResult):
            volume = PocketVolumeComparison(
                Availability.NOT_APPLICABLE,
                diagnostics=(
                    Diagnostic(
                        "pocket.volume.incomplete",
                        DiagnosticSeverity.WARNING,
                        "Pocket volume delta requires a result on both sides.",
                    ),
                ),
            )
        elif isinstance(reference_volume, PocketVolumeComparison) != isinstance(target_volume, PocketVolumeComparison):
            volume = PocketVolumeComparison(
                Availability.NOT_APPLICABLE,
                diagnostics=(
                    Diagnostic(
                        "pocket.volume.incompatible_settings",
                        DiagnosticSeverity.WARNING,
                        "Pocket volume comparisons cannot be combined with single measurements.",
                    ),
                ),
            )
        elif isinstance(reference_volume, (int, float)) and isinstance(target_volume, (int, float)):
            volume = _volume_delta(
                _volume_from_input(reference_volume, "reference"), _volume_from_input(target_volume, "target")
            )
        else:
            volume = PocketVolumeComparison(Availability.NOT_APPLICABLE)
    else:
        volume = None
    if not authoritative or not comparable_match:
        volume = None
    surface_values = _surface_delta(reference_surface_area_angstrom2, target_surface_area_angstrom2)
    surface_delta, relative_surface, reference_surface, target_surface, surface_diagnostics = surface_values[:5]
    (
        reference_surface_method,
        target_surface_method,
        reference_surface_provenance,
        target_surface_provenance,
        reference_surface_units,
        target_surface_units,
    ) = surface_values[5:]
    if not authoritative or not comparable_match:
        surface_delta = None
        relative_surface = None
    qc_states: list[Availability] = []
    qc_items: list[Diagnostic] = []
    for item in (qc, reference_qc, target_qc):
        qc_state, item_diagnostics = _qc_values(item)
        if qc_state is not None:
            qc_states.append(qc_state)
        qc_items.extend(item_diagnostics)
    if qc_diagnostics is not None:
        extra_state, extra_diagnostics = _qc_values(qc_diagnostics)
        if extra_state is not None:
            qc_states.append(extra_state)
        qc_items.extend(extra_diagnostics)
    qc_availability = None
    if qc_states:
        qc_availability = (
            Availability.NUMERICAL_FAILURE
            if Availability.NUMERICAL_FAILURE in qc_states
            else Availability.INVALID_INPUT
            if Availability.INVALID_INPUT in qc_states
            else Availability.AVAILABLE
            if all(item is Availability.AVAILABLE for item in qc_states)
            else Availability.DEPENDENCY_UNAVAILABLE
            if Availability.DEPENDENCY_UNAVAILABLE in qc_states
            else Availability.NOT_APPLICABLE
            if Availability.NOT_APPLICABLE in qc_states
            else Availability.NOT_DETECTED
        )
    diagnostics = list(match.diagnostics) + list(cast(tuple[Diagnostic, ...], surface_diagnostics)) + qc_items
    if not authoritative and comparable_match:
        diagnostics.append(
            Diagnostic(
                "pocket.compare.no_correspondence",
                DiagnosticSeverity.WARNING,
                "Pocket lining and dependent deltas are unavailable without an authoritative correspondence.",
            )
        )
    if qc_availability is not None and qc_availability is not Availability.AVAILABLE:
        availability = cast(Availability, qc_availability)
    elif match.state in {PocketMatchState.UNMATCHED_REFERENCE, PocketMatchState.UNMATCHED_TARGET} or not authoritative:
        availability = Availability.NOT_APPLICABLE
    elif volume is not None and volume.availability is not Availability.AVAILABLE:
        availability = volume.availability
        diagnostics.extend(volume.diagnostics)
    else:
        availability = Availability.AVAILABLE
    return PocketComparison(
        availability=availability,
        match=match,
        volume=volume,
        surface_delta_angstrom2=cast(float | None, surface_delta),
        relative_surface_delta_fraction=cast(float | None, relative_surface),
        reference_surface_area_angstrom2=cast(float | None, reference_surface),
        target_surface_area_angstrom2=cast(float | None, target_surface),
        reference_surface_method=cast(str | None, reference_surface_method),
        target_surface_method=cast(str | None, target_surface_method),
        reference_surface_provenance=cast(Any, reference_surface_provenance),
        target_surface_provenance=cast(Any, target_surface_provenance),
        reference_surface_units=cast(Any, reference_surface_units),
        target_surface_units=cast(Any, target_surface_units),
        lining_residue_conserved=tuple(conserved),
        lining_residue_gains=tuple(gains),
        lining_residue_losses=tuple(losses),
        associated_mutations=associated,
        interaction_changes=interaction_values,
        local_displacement_angstrom=local_displacement,
        ca_displacement_angstrom=local_displacement,
        sidechain_displacement_angstrom=sidechain_displacement,
        local_displacements=relevant_displacements,
        qc_availability=qc_availability,
        qc_diagnostics=tuple(qc_items),
        diagnostics=tuple(diagnostics),
    )


def compare_pocket_candidates(
    reference_candidate: PocketCandidate,
    target_candidate: PocketCandidate,
    correspondences: Sequence[ResidueCorrespondence] | Mapping[ResidueId, ResidueId] = (),
    **kwargs: object,
) -> PocketComparison:
    if not isinstance(reference_candidate, PocketCandidate) or not isinstance(target_candidate, PocketCandidate):
        raise TypeError("candidate inputs must be PocketCandidate values")
    _validate_pre_matching_inputs(correspondences, kwargs)
    transform = kwargs.pop("transform", None)
    settings = kwargs.pop("settings", None)
    result = match_pocket_candidates(
        (reference_candidate,),
        (target_candidate,),
        correspondences,
        transform=cast(Any, transform),
        settings=cast(Any, settings),
    )
    match = (
        result.matches[0]
        if result.matches
        else PocketMatch(reference_candidate, target_candidate, PocketMatchState.MATCHED)
    )
    return compare_pocket_match(match, correspondences=correspondences, **kwargs)  # type: ignore[arg-type]


compare_pockets = compare_pocket_candidates

__all__ = [
    "PocketComparison",
    "PocketSurfaceMeasurement",
    "compare_pocket_candidates",
    "compare_pocket_match",
    "compare_pockets",
]
