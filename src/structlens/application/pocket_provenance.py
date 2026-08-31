"""Deterministic provenance snapshots for application pocket workflows."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from importlib.metadata import version as package_version
from typing import cast

from structlens.core.models import ResidueId
from structlens.core.parsing import ParsedStructure
from structlens.core.pockets import (
    POCKET_LIGAND_RULES_VERSION,
    POCKET_RADII_VERSION,
    FocusedPocketSelection,
    PocketCandidate,
    PocketDetectionSettings,
    PocketVolumeSettings,
)
from structlens.core.provenance import FrozenJSON, MethodProvenance
from structlens.core.sites import SiteDefinition

VOLUME_HYDROGEN_POLICY = "deposited_heavy_atoms_only"

_DETECTION_METHOD_ID = "structlens.pocket.detect"
_DETECTION_METHOD_VERSION = "0.4.0"
_VOLUME_METHOD_ID = "structlens.pocket.volume"
_VOLUME_METHOD_VERSION = "0.4.0"
_FOCUS_METHOD_ID = "structlens.pocket.focus"
_FOCUS_METHOD_VERSION = "0.4.0"


def build_detection_provenance(
    parsed: ParsedStructure,
    settings: PocketDetectionSettings,
    *,
    input_atoms: int,
) -> MethodProvenance:
    """Bind blind detection to its selected source and declared settings."""

    settings_json = settings.to_json()
    parameters: dict[str, object] = {
        **settings_json,
        "input_atoms": input_atoms,
        "selection_id": parsed.selection.selection_id,
        "selection": {
            "model_id": parsed.selection.model_id,
            "author_chain_ids": parsed.selection.author_chain_ids,
            "label_chain_ids": parsed.selection.label_chain_ids,
            "chain_locators": tuple(
                {
                    "author_chain_id": locator.author_chain_id,
                    "label_chain_id": locator.label_chain_id,
                    "entity_id": locator.entity_id,
                }
                for locator in parsed.selection.chain_locators
            ),
            "altloc_policy": parsed.selection.altloc_policy.value,
            "assembly_scope": parsed.selection.assembly_scope.value,
        },
        "atom_scope": "selected_polymer_heavy_atoms",
        "component_scope": "selected_model_and_chain_polymer_residue_records",
        "solvent_exposure_policy": "grid_probe_connectivity_with_boundary_rejection",
    }
    return MethodProvenance(
        method_id=_DETECTION_METHOD_ID,
        method_version=_DETECTION_METHOD_VERSION,
        parameters=cast(Mapping[str, FrozenJSON], parameters),
        units=_settings_units(settings_json),
        backend_versions={
            "scipy": package_version("scipy"),
            "pocket_radii": POCKET_RADII_VERSION,
        },
        input_hashes={
            "raw_source": parsed.raw_source_hash,
            "logical_content": parsed.selection.content_id,
        },
        analyzed_representation=parsed.selection.assembly_scope.value,
    )


def build_volume_provenance(
    parsed: ParsedStructure,
    candidate: PocketCandidate | None,
    settings: PocketVolumeSettings,
    *,
    polymer_atom_count: int,
    component_ids: tuple[str, ...],
) -> MethodProvenance:
    """Build deterministic provenance for a bounded volume measurement."""

    settings_json = settings.to_json()
    selection = parsed.selection
    parameters: dict[str, object] = {
        **settings_json,
        "source_content_id": selection.content_id,
        "selection_id": selection.selection_id,
        "selection": _selection_snapshot(parsed),
        "model_id": selection.model_id,
        "author_chain_ids": selection.author_chain_ids,
        "label_chain_ids": selection.label_chain_ids,
        "altloc_policy": selection.altloc_policy.value,
        "assembly_scope": selection.assembly_scope.value,
        "atom_scope": "selected_primary_polymer_heavy_atoms",
        "component_scope": "selected_retained_ligand_ion_other_components",
        "component_rules_version": POCKET_LIGAND_RULES_VERSION,
        "hydrogen_policy": VOLUME_HYDROGEN_POLICY,
        "polymer_atom_count": polymer_atom_count,
        "component_ids": component_ids,
        **_candidate_snapshot(candidate),
    }
    candidate_hashes: dict[str, str] = {}
    if candidate is not None:
        if candidate.source_content_id is None:
            raise ValueError("candidate requires source content lineage")
        candidate_hashes = {
            "candidate": candidate.candidate_id,
            "source_content": candidate.source_content_id,
        }
    return MethodProvenance(
        method_id=_VOLUME_METHOD_ID,
        method_version=_VOLUME_METHOD_VERSION,
        parameters=cast(Mapping[str, FrozenJSON], parameters),
        units={
            "volume": "angstrom^3",
            "length": "angstrom",
            "coarse_grid_spacing_angstrom": "angstrom",
            "fine_grid_spacing_angstrom": "angstrom",
            "boundary_margin_angstrom": "angstrom",
            "sphere_radii_angstrom": "angstrom",
        },
        backend_versions={
            "numpy": package_version("numpy"),
            "scipy": package_version("scipy"),
            "pocket_radii": POCKET_RADII_VERSION,
        },
        input_hashes={
            "raw_source": parsed.raw_source_hash,
            "logical_content": selection.content_id,
            **candidate_hashes,
        },
        analyzed_representation=selection.assembly_scope.value,
    )


def build_volume_result_parameters(
    parsed: ParsedStructure,
    candidate: PocketCandidate | None,
    measured_parameters: Mapping[str, object],
) -> dict[str, object]:
    """Snapshot producer lineage beside the immutable method provenance."""

    return {
        **dict(measured_parameters),
        "raw_source_hash": parsed.raw_source_hash,
        "logical_content": parsed.selection.content_id,
        "altloc_policy": parsed.selection.altloc_policy.value,
        "hydrogen_policy": VOLUME_HYDROGEN_POLICY,
        **_sphere_snapshot(candidate),
    }


def build_focus_provenance(
    parsed: ParsedStructure,
    candidates: Sequence[PocketCandidate],
    definition: SiteDefinition,
    result: FocusedPocketSelection,
) -> MethodProvenance:
    """Bind a focused-pocket choice to its inputs, rules, and support evidence."""

    selection = parsed.selection
    parameters: dict[str, object] = {
        "selection_id": selection.selection_id,
        "model_id": selection.model_id,
        "author_chain_ids": selection.author_chain_ids,
        "label_chain_ids": selection.label_chain_ids,
        "altloc_policy": selection.altloc_policy.value,
        "assembly_scope": selection.assembly_scope.value,
        "hydrogen_policy": VOLUME_HYDROGEN_POLICY,
        "focused_site_atom_scope": "selected_primary_heavy_atoms",
        "ligand_rules_version": POCKET_LIGAND_RULES_VERSION,
        "site_definition": {
            "site_id": definition.site_id,
            "name": definition.name,
            "mode": definition.mode.value,
            "reference_residues": tuple(_residue_snapshot(item) for item in definition.reference_residues),
            "center_residue": (
                _residue_snapshot(definition.center_residue) if definition.center_residue is not None else None
            ),
            "ligand_id": definition.ligand_id,
            "radius_angstrom": definition.radius_angstrom,
        },
        "candidate_ids": tuple(sorted(candidate.candidate_id for candidate in candidates)),
        "candidate_lineage": tuple(
            {
                "candidate_id": candidate.candidate_id,
                "source_content_id": candidate.source_content_id,
                "selection_id": candidate.selection_id,
            }
            for candidate in sorted(candidates, key=lambda item: item.candidate_id)
        ),
        "selected_candidate_id": result.candidate.candidate_id if result.candidate is not None else None,
        "availability": result.availability.value,
        "ligand_support": result.ligand_support.to_json() if result.ligand_support is not None else None,
    }
    return MethodProvenance(
        method_id=_FOCUS_METHOD_ID,
        method_version=_FOCUS_METHOD_VERSION,
        parameters=cast(Mapping[str, FrozenJSON], parameters),
        units={
            "site_definition.radius_angstrom": "angstrom",
            "ligand_support.ligand_center_distance_angstrom": "angstrom",
            "ligand_support.atom_coverage_fraction": "fraction",
            "ligand_support.lining_residue_overlap_fraction": "fraction",
        },
        backend_versions={"pocket_ligand_rules": POCKET_LIGAND_RULES_VERSION},
        input_hashes={
            "raw_source": parsed.raw_source_hash,
            "logical_content": selection.content_id,
        },
        analyzed_representation=selection.assembly_scope.value,
    )


def _selection_snapshot(parsed: ParsedStructure) -> dict[str, object]:
    selection = parsed.selection
    return {
        "selection_id": selection.selection_id,
        "format": selection.format.value,
        "model_id": selection.model_id,
        "author_chain_ids": selection.author_chain_ids,
        "label_chain_ids": selection.label_chain_ids,
        "chain_locators": tuple(
            {
                "author_chain_id": locator.author_chain_id,
                "label_chain_id": locator.label_chain_id,
                "entity_id": locator.entity_id,
            }
            for locator in selection.chain_locators
        ),
        "altloc_policy": selection.altloc_policy.value,
        "assembly_scope": selection.assembly_scope.value,
    }


def _candidate_snapshot(candidate: PocketCandidate | None) -> dict[str, object]:
    return {
        "candidate_id": candidate.candidate_id if candidate is not None else None,
        "candidate_lineage": (
            {
                "candidate_id": candidate.candidate_id,
                "source_content_id": candidate.source_content_id,
                "selection_id": candidate.selection_id,
            }
            if candidate is not None
            else None
        ),
        **_sphere_snapshot(candidate),
    }


def _sphere_snapshot(candidate: PocketCandidate | None) -> dict[str, object]:
    spheres = candidate.alpha_spheres if candidate is not None else ()
    return {
        "sphere_count": len(spheres),
        "sphere_ids": tuple(sphere.sphere_id for sphere in spheres),
        "sphere_radii_angstrom": tuple(sphere.radius_angstrom for sphere in spheres),
    }


def _residue_snapshot(residue: ResidueId) -> dict[str, str | None]:
    return {
        "structure_id": residue.structure_id,
        "model_id": residue.model_id,
        "chain_id": residue.chain_id,
        "auth_seq_id": residue.auth_seq_id,
        "insertion_code": residue.insertion_code,
        "residue_name": residue.residue_name,
    }


def _settings_units(payload: Mapping[str, object], prefix: str = "") -> dict[str, str]:
    units: dict[str, str] = {}
    for key, value in payload.items():
        path = f"{prefix}.{key}" if prefix else key
        if isinstance(value, Mapping):
            units.update(_settings_units(value, path))
            continue
        leaf = key.casefold()
        if "angstrom" in leaf or leaf.endswith("_radius") or "spacing" in leaf or "padding" in leaf or "margin" in leaf:
            units[path] = "angstrom"
        elif (
            "count" in leaf
            or "candidate" in leaf
            or "atom" in leaf
            or "simplex" in leaf
            or "cluster_size" in leaf
            or "cell" in leaf
            or "check" in leaf
        ):
            units[path] = "count"
    return units


__all__ = [
    "VOLUME_HYDROGEN_POLICY",
    "build_detection_provenance",
    "build_focus_provenance",
    "build_volume_provenance",
    "build_volume_result_parameters",
]
