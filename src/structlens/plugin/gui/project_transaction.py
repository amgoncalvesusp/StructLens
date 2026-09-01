"""Transactional canonical project persistence for the desktop adapter."""

from __future__ import annotations

import os
import shutil
import tempfile
from collections.abc import Mapping
from dataclasses import asdict, dataclass
from pathlib import Path

from structlens.application.dto import AnalysisReportRequest
from structlens.application.project_state import ProjectState
from structlens.application.report_serialization import load_report, save_report
from structlens.application.report_snapshot_io import load_snapshot, store_snapshot
from structlens.core.errors import ProjectSchemaError
from structlens.core.interactions import InteractionThresholds
from structlens.core.models import ResidueId
from structlens.core.msa import MSASettings
from structlens.core.reports.safe_io import validate_safe_ancestors
from structlens.core.sites import SiteDefinition, SiteDefinitionMode
from structlens.plugin.visualization.renderer import VisualizationState, visualization_state_from_mapping

from .loaded_source import LoadedSource
from .model import CanonicalReportBinding


@dataclass(frozen=True, slots=True)
class CanonicalProjectCandidate:
    project: ProjectState
    binding: CanonicalReportBinding
    reference_source: LoadedSource
    target_source: LoadedSource
    visualization_state: VisualizationState


@dataclass(frozen=True, slots=True)
class LegacyProjectCandidate:
    project: ProjectState


ProjectCandidate = CanonicalProjectCandidate | LegacyProjectCandidate


def save_canonical_project(
    binding: CanonicalReportBinding,
    project_path: str | Path,
    *,
    visualization_state: Mapping[str, object],
) -> ProjectState:
    """Persist one verified generation, committing project JSON last."""

    if not isinstance(binding, CanonicalReportBinding):
        raise TypeError("a verified canonical report binding is required")
    visual_state = parse_visualization_state(visualization_state)
    output = Path(project_path)
    artifact_parent = Path(f"{output}.artifacts")
    generation = artifact_parent / binding.report.report_id
    validate_safe_ancestors(artifact_parent, label="project artifacts")
    artifact_parent.mkdir(parents=True, exist_ok=True)
    validate_safe_ancestors(artifact_parent, label="project artifacts")

    if generation.exists():
        _verify_generation(binding, generation)
    else:
        staging = Path(tempfile.mkdtemp(prefix=".structlens-stage-", dir=artifact_parent))
        try:
            _write_generation(binding, staging)
            _verify_generation(binding, staging)
            try:
                os.replace(staging, generation)
            except OSError:
                if not generation.exists():
                    raise
                _verify_generation(binding, generation)
        finally:
            if staging.exists():
                shutil.rmtree(staging)

    references = _generation_references(output, generation)
    report_relative = _relative_to_project(output, generation / "report.json")
    request = binding.request
    project = ProjectState(
        settings=request.analysis_settings,
        analysis_results=(),
        source_hashes={
            "reference_content_id": request.reference_snapshot.content_id,
            "reference_raw_sha256": request.reference_snapshot.raw_sha256,
            "target_content_id": request.target_snapshot.content_id,
            "target_raw_sha256": request.target_snapshot.raw_sha256,
        },
        visualization_state=visualization_payload(visual_state),
        comparison_mode=binding.report.comparison_mode,
        msa_settings=asdict(request.msa_settings),
        interaction_thresholds=asdict(request.interaction_thresholds),
        site_definitions=tuple(_site_to_payload(item) for item in request.site_definitions),
        evidence_sources={"canonical_request": _request_metadata(request)},
        reference_selection=request.reference_selection,
        target_selection=request.target_selection,
        input_selections={
            "reference": request.reference_selection,
            "target": request.target_selection,
        },
        report_references=references,
        report_path=report_relative,
        report_hash=binding.report.report_id,
        report_verification="verified",
    )
    project.verify_report(
        output.parent / report_relative,
        snapshots=binding.snapshots,
    )
    project.save(output)
    return project


def stage_project(project_path: str | Path) -> ProjectCandidate:
    """Load and fully validate a candidate without mutating GUI state."""

    source = Path(project_path)
    project = ProjectState.load(source)
    if not _has_canonical_metadata(project):
        return LegacyProjectCandidate(project)
    if (
        project.report_verification != "verified"
        or project.report_hash is None
        or project.report_path is None
        or project.reference_selection is None
        or project.target_selection is None
    ):
        raise ProjectSchemaError("canonical project metadata is incomplete or unverified")

    generation = Path(f"{source}.artifacts") / project.report_hash
    _validate_generation_root(generation)
    report_path = _validated_reference(source, generation, project.report_path)
    if tuple(project.report_references) != tuple(sorted(set(project.report_references))):
        raise ProjectSchemaError("canonical artifact references must be unique and sorted")
    if project.report_path not in project.report_references:
        raise ProjectSchemaError("canonical report is missing from artifact references")
    for reference in project.report_references:
        artifact = _validated_reference(source, generation, reference)
        _validate_regular_artifact(artifact)

    reference_snapshot = load_snapshot(
        project.reference_selection.content_id,
        _snapshot_directory(generation, "reference"),
        _source_hash(project, "reference_raw_sha256"),
    )
    target_snapshot = load_snapshot(
        project.target_selection.content_id,
        _snapshot_directory(generation, "target"),
        _source_hash(project, "target_raw_sha256"),
    )
    snapshots = (reference_snapshot, target_snapshot)
    report = load_report(report_path, snapshots=snapshots)
    if report.report_id != project.report_hash:
        raise ProjectSchemaError("canonical report ID does not match project metadata")
    if report.reference_selection.selection_id != project.reference_selection.selection_id:
        raise ProjectSchemaError("reference selection does not match the canonical report")
    if report.target_selection.selection_id != project.target_selection.selection_id:
        raise ProjectSchemaError("target selection does not match the canonical report")

    request = _restore_request(project, snapshots)
    binding = CanonicalReportBinding.create(report, request)
    visual_state = parse_visualization_state(project.visualization_state)
    reference_source = LoadedSource.from_snapshot(
        project.reference_selection.path or reference_snapshot.display_name,
        reference_snapshot,
    )
    target_source = LoadedSource.from_snapshot(
        project.target_selection.path or target_snapshot.display_name,
        target_snapshot,
    )
    return CanonicalProjectCandidate(project, binding, reference_source, target_source, visual_state)


def _write_generation(binding: CanonicalReportBinding, generation: Path) -> None:
    store_snapshot(binding.snapshots[0], _snapshot_directory(generation, "reference"))
    store_snapshot(binding.snapshots[1], _snapshot_directory(generation, "target"))
    save_report(binding.report, generation / "report.json", snapshots=binding.snapshots)


def _verify_generation(binding: CanonicalReportBinding, generation: Path) -> None:
    _validate_generation_root(generation)
    snapshots = tuple(
        load_snapshot(
            snapshot.content_id,
            _snapshot_directory(generation, role),
            snapshot.raw_sha256,
        )
        for role, snapshot in zip(("reference", "target"), binding.snapshots, strict=True)
    )
    loaded = load_report(generation / "report.json", snapshots=snapshots)
    if loaded.report_id != binding.report.report_id:
        raise ValueError("existing canonical artifact generation conflicts with this report")


def _generation_references(project: Path, generation: Path) -> tuple[str, ...]:
    references = tuple(
        sorted(_relative_to_project(project, artifact) for artifact in generation.rglob("*") if artifact.is_file())
    )
    if not references:
        raise ValueError("canonical artifact generation is empty")
    return references


def _relative_to_project(project: Path, artifact: Path) -> str:
    return artifact.absolute().relative_to(project.parent.absolute()).as_posix()


def _validated_reference(project: Path, generation: Path, value: str) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise ProjectSchemaError("canonical artifact reference must be a non-empty relative path")
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        raise ProjectSchemaError("canonical artifact reference escapes its generation")
    candidate = (project.parent / relative).absolute()
    try:
        candidate.relative_to(generation.absolute())
    except ValueError as exc:
        raise ProjectSchemaError("canonical artifact reference escapes its generation") from exc
    return candidate


def _validate_generation_root(generation: Path) -> None:
    validate_safe_ancestors(generation, label="canonical artifact generation")
    if not generation.is_dir() or _is_link_or_reparse(generation):
        raise ProjectSchemaError("canonical artifact generation is missing or unsafe")


def _validate_regular_artifact(path: Path) -> None:
    validate_safe_ancestors(path.parent, label="canonical artifact")
    if not path.is_file() or _is_link_or_reparse(path):
        raise ProjectSchemaError(f"canonical artifact is missing or unsafe: {path.name}")


def _is_link_or_reparse(path: Path) -> bool:
    try:
        stat = path.lstat()
    except OSError:
        return True
    reparse_flag = getattr(stat, "st_file_attributes", 0) & getattr(stat, "FILE_ATTRIBUTE_REPARSE_POINT", 0x400)
    return path.is_symlink() or bool(reparse_flag)


def _has_canonical_metadata(project: ProjectState) -> bool:
    return bool(
        project.report_path
        or project.report_hash
        or project.report_references
        or project.reference_selection
        or project.target_selection
    )


def _source_hash(project: ProjectState, name: str) -> str:
    value = project.source_hashes.get(name)
    if not isinstance(value, str) or not value:
        raise ProjectSchemaError(f"canonical project is missing {name}")
    return value


def _snapshot_directory(generation: Path, role: str) -> Path:
    root = generation / "snapshots"
    return root if role == "reference" else root / role


def _request_metadata(request: AnalysisReportRequest) -> dict[str, object]:
    return {
        "minimum_vector_magnitude_angstrom": request.minimum_vector_magnitude_angstrom,
        "maximum_vectors": request.maximum_vectors,
        "manual_pairs": [
            {"reference": asdict(reference), "target": asdict(target)} for reference, target in request.manual_pairs
        ],
    }


def _restore_request(
    project: ProjectState,
    snapshots: tuple[object, object],
) -> AnalysisReportRequest:
    reference_snapshot, target_snapshot = snapshots
    if project.reference_selection is None or project.target_selection is None:
        raise ProjectSchemaError("canonical project selections are missing")
    metadata = project.evidence_sources.get("canonical_request", {})
    if not isinstance(metadata, Mapping):
        raise ProjectSchemaError("canonical request metadata is invalid")
    manual_values = metadata.get("manual_pairs", ())
    if not isinstance(manual_values, (list, tuple)):
        raise ProjectSchemaError("canonical manual pairs metadata is invalid")
    manual_pairs = tuple(_manual_pair(value) for value in manual_values)
    sites = tuple(_site_from_payload(value) for value in project.site_definitions)
    try:
        return AnalysisReportRequest(
            reference_snapshot,  # type: ignore[arg-type]
            target_snapshot,  # type: ignore[arg-type]
            project.reference_selection,
            project.target_selection,
            analysis_settings=project.settings,
            msa_settings=MSASettings(**dict(project.msa_settings)),
            interaction_thresholds=InteractionThresholds(**dict(project.interaction_thresholds)),
            site_definitions=sites,
            manual_pairs=manual_pairs,
            minimum_vector_magnitude_angstrom=float(metadata.get("minimum_vector_magnitude_angstrom", 0.5)),
            maximum_vectors=int(metadata.get("maximum_vectors", 100)),
        )
    except (TypeError, ValueError, OverflowError) as exc:
        raise ProjectSchemaError(f"canonical request metadata is invalid: {exc}") from exc


def _manual_pair(value: object) -> tuple[ResidueId, ResidueId]:
    if not isinstance(value, Mapping):
        raise ProjectSchemaError("canonical manual pair must be an object")
    try:
        return ResidueId(**dict(value["reference"])), ResidueId(**dict(value["target"]))
    except (KeyError, TypeError, ValueError) as exc:
        raise ProjectSchemaError("canonical manual pair is invalid") from exc


def _site_to_payload(site: SiteDefinition) -> dict[str, object]:
    return {
        "site_id": site.site_id,
        "name": site.name,
        "mode": site.mode.value,
        "reference_residues": [asdict(item) for item in site.reference_residues],
        "center_residue": None if site.center_residue is None else asdict(site.center_residue),
        "ligand_id": site.ligand_id,
        "radius_angstrom": site.radius_angstrom,
    }


def _site_from_payload(value: object) -> SiteDefinition:
    if not isinstance(value, Mapping):
        raise ProjectSchemaError("canonical site definition must be an object")
    try:
        residues = tuple(ResidueId(**dict(item)) for item in value.get("reference_residues", ()))
        center = value.get("center_residue")
        return SiteDefinition(
            str(value["site_id"]),
            name=str(value["name"]),
            mode=SiteDefinitionMode(str(value["mode"])),
            reference_residues=residues,
            center_residue=None if center is None else ResidueId(**dict(center)),
            ligand_id=None if value.get("ligand_id") is None else str(value["ligand_id"]),
            radius_angstrom=None if value.get("radius_angstrom") is None else float(value["radius_angstrom"]),
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise ProjectSchemaError("canonical site definition is invalid") from exc


def parse_visualization_state(payload: object) -> VisualizationState:
    try:
        return visualization_state_from_mapping(payload)
    except (TypeError, ValueError) as exc:
        raise ProjectSchemaError(f"visualization state is invalid: {exc}") from exc


def visualization_payload(state: VisualizationState) -> dict[str, object]:
    return {
        "highlight_filter": state.highlight_filter.value,
        "color_mode": state.color_mode.value,
        "representation": state.representation.value,
        "show_labels": state.show_labels,
        "show_reference": state.show_reference,
        "show_target": state.show_target,
        "local_radius_angstrom": state.local_radius_angstrom,
        "preset": state.preset,
    }


__all__ = [
    "CanonicalProjectCandidate",
    "LegacyProjectCandidate",
    "ProjectCandidate",
    "parse_visualization_state",
    "save_canonical_project",
    "stage_project",
    "visualization_payload",
]
