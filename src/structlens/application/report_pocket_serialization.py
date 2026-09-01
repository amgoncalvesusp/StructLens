"""Typed wire reconstruction for pocket evidence in canonical reports."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, cast

from structlens.core.evidence import Availability
from structlens.core.pockets.comparison_models import PocketComparison
from structlens.core.pockets.matching import PocketMatch, PocketMatchingResult, PocketMatchingSettings, PocketMatchState
from structlens.core.pockets.models import (
    AlphaSphere,
    AlphaSphereDetectionResult,
    PocketCandidate,
    PocketDetectionSettings,
    PocketGeometrySettings,
)
from structlens.core.pockets.volume_models import (
    PocketVolumeComparison,
    PocketVolumeResult,
    PocketVolumeSensitivity,
    PocketVolumeSettings,
)
from structlens.core.provenance import FrozenJSON, MethodProvenance, freeze_json
from structlens.core.reports import (
    PocketDetectionSnapshot,
    PocketLiningSnapshot,
    PocketReportSnapshot,
    PocketVolumeSnapshot,
)


def pocket_report_from_dict(payload: Mapping[str, Any]) -> PocketReportSnapshot:
    """Decode every supported pocket subtype through existing bounded wire helpers."""

    from structlens.application import report_serialization as wire

    return PocketReportSnapshot(
        detections=tuple(
            _detection(wire._object(item, "pocket detection"), wire) for item in wire._array(payload, "detections")
        ),
        volumes=tuple(
            _volume_snapshot(wire._object(item, "pocket volume"), wire) for item in wire._array(payload, "volumes")
        ),
        matching=None
        if payload.get("matching") is None
        else _matching(wire._object(payload["matching"], "pocket matching"), wire),
        comparisons=tuple(
            _comparison(wire._object(item, "pocket comparison"), wire) for item in wire._array(payload, "comparisons")
        ),
        lining_residues=tuple(
            _lining(wire._object(item, "pocket lining"), wire) for item in wire._array(payload, "lining_residues")
        ),
        concordance=wire._pocket_sections_from_dict(payload.get("concordance")),
        diagnostics=tuple(
            wire._diagnostic_from_dict(wire._object(item, "pocket diagnostic"))
            for item in wire._array(payload, "diagnostics")
        ),
        provenance=None
        if payload.get("provenance") is None
        else MethodProvenance.from_json(wire._object(payload["provenance"], "pocket provenance")),
        units=dict(wire._object(payload.get("units", {}), "pocket units")),
        availability=Availability(str(payload.get("availability", "not_applicable"))),
    )


def _detection(payload: Mapping[str, Any], wire: Any) -> PocketDetectionSnapshot:
    detection_value = payload.get("detection")
    detection = None
    if detection_value is not None:
        item = wire._object(detection_value, "pocket detection result")
        detection = AlphaSphereDetectionResult(
            tuple(_sphere(wire._object(value, "alpha sphere"), wire) for value in wire._array(item, "spheres")),
            tuple(
                wire._diagnostic_from_dict(wire._object(value, "detection diagnostic"))
                for value in wire._array(item, "diagnostics")
            ),
        )
    settings_value = payload.get("settings")
    settings = (
        None if settings_value is None else _detection_settings(wire._object(settings_value, "detection settings"))
    )
    return PocketDetectionSnapshot(
        str(payload["role"]),
        Availability(str(payload["availability"])),
        tuple(_candidate(wire._object(item, "pocket candidate"), wire) for item in wire._array(payload, "candidates")),
        detection,
        settings,
        tuple(
            wire._diagnostic_from_dict(wire._object(item, "detection diagnostic"))
            for item in wire._array(payload, "diagnostics")
        ),
        {str(key): int(value) for key, value in wire._object(payload.get("counts", {}), "detection counts").items()},
        None
        if payload.get("provenance") is None
        else MethodProvenance.from_json(wire._object(payload["provenance"], "detection provenance")),
    )


def _detection_settings(payload: Mapping[str, Any]) -> PocketDetectionSettings:
    geometry = payload.get("geometry", {})
    return PocketDetectionSettings(
        geometry=PocketGeometrySettings(**dict(geometry)),
        **{str(key): value for key, value in payload.items() if key != "geometry"},
    )


def _sphere(payload: Mapping[str, Any], wire: Any) -> AlphaSphere:
    return AlphaSphere(
        wire._coordinate(payload, "center_xyz"),
        float(payload["radius_angstrom"]),
        tuple(str(item) for item in wire._array(payload, "touching_atom_ids")),
        tuple(wire._required_residue(item, "lining residue") for item in wire._array(payload, "lining_residues")),
        tuple(str(item) for item in wire._array(payload, "source_simplex_atom_ids")),
    )


def _candidate(payload: Mapping[str, Any], wire: Any) -> PocketCandidate:
    lineage = wire._object(payload.get("lineage", {}), "candidate lineage")
    candidate = PocketCandidate(
        tuple(
            _sphere(wire._object(item, "candidate alpha sphere"), wire)
            for item in wire._array(payload, "alpha_spheres")
        ),
        lineage.get("source_content_id"),
        lineage.get("selection_id"),
    )
    if payload.get("candidate_id") != candidate.candidate_id:
        raise ValueError("candidate_id does not match canonical candidate content")
    return candidate


def _volume_snapshot(payload: Mapping[str, Any], wire: Any) -> PocketVolumeSnapshot:
    result_value = payload.get("result")
    result = None if result_value is None else _volume_result(wire._object(result_value, "volume result"), wire)
    return PocketVolumeSnapshot(str(payload["role"]), result, Availability(str(payload["availability"])))


def _volume_result(payload: Mapping[str, Any], wire: Any) -> PocketVolumeResult:
    grid = wire._object(payload.get("grid", {}), "volume grid")
    coarse = wire._object(payload.get("coarse", {}), "coarse volume")
    fine = wire._object(payload.get("fine", {}), "fine volume")
    settings = payload.get("settings")
    sensitivity = payload.get("sensitivity")
    return PocketVolumeResult(
        Availability(str(payload["availability"])),
        coarse_voxel_count=wire._optional_int(payload.get("coarse_voxel_count", coarse.get("voxel_count"))),
        fine_voxel_count=wire._optional_int(payload.get("fine_voxel_count", fine.get("voxel_count"))),
        coarse_volume_angstrom3=wire._optional_float(
            payload.get("coarse_volume_angstrom3", coarse.get("volume_angstrom3"))
        ),
        fine_volume_angstrom3=wire._optional_float(payload.get("fine_volume_angstrom3", fine.get("volume_angstrom3"))),
        grid_shape=None if grid.get("shape") is None else _shape(grid["shape"], wire, "grid shape"),
        grid_origin_xyz=None if grid.get("origin_xyz") is None else wire._coordinate(grid, "origin_xyz"),
        grid_phase=str(grid.get("phase", "cell_center")),
        coarse_grid_shape=None
        if grid.get("coarse_shape") is None
        else _shape(grid["coarse_shape"], wire, "coarse grid shape"),
        fine_grid_shape=None if grid.get("fine_shape") is None else _shape(grid["fine_shape"], wire, "fine grid shape"),
        rotation_error_bound_angstrom3=wire._optional_float(payload.get("rotation_error_bound_angstrom3")),
        sensitivity=None
        if sensitivity is None
        else PocketVolumeSensitivity(**wire._object(sensitivity, "volume sensitivity")),
        units=dict(wire._object(payload.get("units", {}), "volume units")),
        provenance_parameters=dict(wire._object(payload.get("provenance_parameters", {}), "volume parameters")),
        diagnostics=tuple(
            wire._diagnostic_from_dict(wire._object(item, "volume diagnostic"))
            for item in wire._array(payload, "diagnostics")
        ),
        candidate_id=payload.get("candidate_id"),
        compatibility_signature=tuple(
            tuple(item) if isinstance(item, (list, tuple)) else item
            for item in wire._array(payload, "compatibility_signature")
        ),
        settings=None
        if settings is None
        else PocketVolumeSettings(
            **{
                str(key): value
                for key, value in wire._object(settings, "volume settings").items()
                if key not in {"grid_phase", "grid_sampling"}
            }
        ),
        provenance=None
        if payload.get("provenance") is None
        else MethodProvenance.from_json(wire._object(payload["provenance"], "volume provenance")),
    )


def _shape(value: object, wire: Any, name: str) -> tuple[int, int, int]:
    values = tuple(int(item) for item in wire._array_value(value, name))
    if len(values) != 3:
        raise ValueError(f"{name} must contain exactly three values")
    return values


def _matching(payload: Mapping[str, Any], wire: Any) -> PocketMatchingResult:
    return PocketMatchingResult(
        Availability(str(payload["availability"])),
        tuple(_match(wire._object(item, "pocket match"), wire) for item in wire._array(payload, "matches")),
        tuple(
            wire._diagnostic_from_dict(wire._object(item, "matching diagnostic"))
            for item in wire._array(payload, "diagnostics")
        ),
        PocketMatchingSettings(**wire._object(payload["settings"], "matching settings")),
        wire._transform_from_dict(wire._object(payload["transform"], "matching transform")),
        wire._optional_float(payload.get("optimal_score")),
        wire._optional_float(payload.get("second_best_score")),
    )


def _match(payload: Mapping[str, Any], wire: Any) -> PocketMatch:
    return PocketMatch(
        None
        if payload.get("reference_candidate") is None
        else _candidate(wire._object(payload["reference_candidate"], "reference candidate"), wire),
        None
        if payload.get("target_candidate") is None
        else _candidate(wire._object(payload["target_candidate"], "target candidate"), wire),
        PocketMatchState(str(payload.get("state", payload.get("status")))),
        wire._optional_float(payload.get("lining_jaccard")),
        wire._optional_float(payload.get("centroid_distance_angstrom")),
        wire._optional_float(payload.get("score")),
        wire._optional_float(payload.get("alternative_score")),
        wire._optional_float(payload.get("ambiguity_margin")),
        tuple(
            wire._diagnostic_from_dict(wire._object(item, "match diagnostic"))
            for item in wire._array(payload, "diagnostics")
        ),
        dict(wire._object(payload.get("units", {}), "match units")),
        wire._optional_float(payload.get("assignment_score")),
        wire._optional_float(payload.get("alternative_assignment_score")),
    )


def _lining(payload: Mapping[str, Any], wire: Any) -> PocketLiningSnapshot:
    return PocketLiningSnapshot(
        str(payload["role"]),
        str(payload["candidate_id"]),
        tuple(wire._required_residue(item, "lining residue") for item in wire._array(payload, "residues")),
        Availability(str(payload["availability"])),
        tuple(
            wire._diagnostic_from_dict(wire._object(item, "lining diagnostic"))
            for item in wire._array(payload, "diagnostics")
        ),
        dict(wire._object(payload.get("units", {}), "lining units")),
    )


def _comparison(payload: Mapping[str, Any], wire: Any) -> PocketComparison:
    associated_mutations = tuple(
        wire._mutation_from_dict(wire._object(item, "associated mutation"))
        for item in wire._array(payload, "associated_mutations")
    )
    match = _match(wire._object(payload["match"], "comparison match"), wire)
    volume_value = payload.get("volume")
    volume = None
    if volume_value is not None:
        value = wire._object(volume_value, "volume comparison")
        volume = PocketVolumeComparison(
            Availability(str(value["availability"])),
            wire._optional_float(value.get("delta_angstrom3")),
            wire._optional_float(value.get("relative_delta_fraction", value.get("relative_delta"))),
            wire._optional_float(value.get("reference_volume_angstrom3")),
            wire._optional_float(value.get("target_volume_angstrom3")),
            tuple(
                wire._diagnostic_from_dict(wire._object(item, "volume comparison diagnostic"))
                for item in wire._array(value, "diagnostics")
            ),
            dict(wire._object(value.get("units", {}), "volume comparison units")),
            tuple(
                tuple(item) if isinstance(item, (list, tuple)) else item
                for item in wire._array(value, "compatibility_signature")
            ),
            None
            if value.get("reference_sensitivity") is None
            else PocketVolumeSensitivity(**wire._object(value["reference_sensitivity"], "reference sensitivity")),
            None
            if value.get("target_sensitivity") is None
            else PocketVolumeSensitivity(**wire._object(value["target_sensitivity"], "target sensitivity")),
            tuple(
                MethodProvenance.from_json(wire._object(item, "comparison provenance"))
                for item in wire._array(value, "provenance")
            ),
        )

    def surface_provenance(value: object) -> FrozenJSON | MethodProvenance | None:
        if value is None:
            return None
        if isinstance(value, MethodProvenance):
            return value
        if not isinstance(value, Mapping):
            raise ValueError("surface provenance must be an object or null")
        if {"method_id", "method_version", "artifact_id"}.issubset(value):
            return MethodProvenance.from_json(cast(Mapping[str, object], value))
        return freeze_json(cast(Mapping[str, object], value), path="surface provenance")

    return PocketComparison(
        Availability(str(payload["availability"])),
        match,
        volume,
        wire._optional_float(payload.get("surface_delta_angstrom2")),
        wire._optional_float(payload.get("relative_surface_delta_fraction")),
        wire._optional_float(payload.get("reference_surface_area_angstrom2")),
        wire._optional_float(payload.get("target_surface_area_angstrom2")),
        payload.get("reference_surface_method"),
        payload.get("target_surface_method"),
        surface_provenance(payload.get("reference_surface_provenance")),
        surface_provenance(payload.get("target_surface_provenance")),
        None
        if payload.get("reference_surface_units") is None
        else dict(wire._object(payload["reference_surface_units"], "reference surface units")),
        None
        if payload.get("target_surface_units") is None
        else dict(wire._object(payload["target_surface_units"], "target surface units")),
        tuple(
            wire._required_residue(item, "conserved lining residue")
            for item in wire._array(payload, "lining_residue_conserved")
        ),
        tuple(
            wire._required_residue(item, "gained lining residue")
            for item in wire._array(payload, "lining_residue_gains")
        ),
        tuple(
            wire._required_residue(item, "lost lining residue")
            for item in wire._array(payload, "lining_residue_losses")
        ),
        associated_mutations,
        tuple(
            wire._interaction_difference_from_dict(wire._object(item, "interaction change"))
            for item in wire._array(payload, "interaction_changes")
        ),
        wire._optional_float(payload.get("local_displacement_angstrom")),
        wire._optional_float(payload.get("ca_displacement_angstrom")),
        wire._optional_float(payload.get("sidechain_displacement_angstrom")),
        tuple(
            (
                wire._required_residue(
                    wire._object(item, "local displacement").get("residue"), "local displacement residue"
                ),
                float(wire._object(item, "local displacement")["magnitude_angstrom"]),
            )
            for item in wire._array(payload, "local_displacements")
        ),
        None if payload.get("qc_availability") is None else Availability(str(payload["qc_availability"])),
        tuple(
            wire._diagnostic_from_dict(wire._object(item, "QC diagnostic"))
            for item in wire._array(payload, "qc_diagnostics")
        ),
        tuple(
            wire._diagnostic_from_dict(wire._object(item, "comparison diagnostic"))
            for item in wire._array(payload, "diagnostics")
        ),
        dict(wire._object(payload.get("units", {}), "comparison units")),
    )


__all__ = ["pocket_report_from_dict"]
