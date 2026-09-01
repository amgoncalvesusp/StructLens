"""Versioned, JSON-serializable project state."""

from __future__ import annotations

import hashlib
import json
import math
import re
from collections.abc import Mapping
from dataclasses import asdict, dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Any

import numpy as np

from structlens.core.errors import ProjectSchemaError
from structlens.core.metrics.sequence_metrics import SequenceAlignmentMetrics
from structlens.core.metrics.structural_metrics import StructuralMetrics
from structlens.core.models import (
    AlignmentMode,
    AnalysisResult,
    AnalysisSettings,
    ComparisonMode,
    CorrespondenceStatus,
    MutationEvent,
    MutationKind,
    ReferenceVsManyAnalysis,
    ResidueCorrespondence,
    ResidueId,
    StructuralTransform,
    TargetAnalysis,
)
from structlens.core.parsing import AltlocPolicy, AssemblyScope, ChainLocator, InputSelection, StructureFormat
from structlens.core.provenance import MethodProvenance
from structlens.core.reports.safe_io import atomic_write_text, read_bounded_bytes

_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_MAX_PROJECT_ITEMS = 100_000
_MAX_PROJECT_BYTES = 100 * 1024 * 1024


@dataclass(frozen=True, slots=True)
class ProjectState:
    reference_source: str | None = None
    target_sources: tuple[str, ...] = ()
    settings: AnalysisSettings = field(default_factory=AnalysisSettings)
    key_residues: tuple[ResidueId, ...] = ()
    analysis_results: tuple[AnalysisResult, ...] = ()
    source_hashes: dict[str, str] = field(default_factory=dict)
    visualization_state: dict[str, Any] = field(default_factory=dict)
    schema_version: str = "3.0"
    source_objects: dict[str, str] = field(default_factory=dict)
    comparison_mode: ComparisonMode = ComparisonMode.PAIRWISE
    reference_vs_many: ReferenceVsManyAnalysis | None = None
    msa_settings: dict[str, Any] = field(default_factory=dict)
    interaction_thresholds: dict[str, float] = field(default_factory=dict)
    site_definitions: tuple[dict[str, Any], ...] = ()
    distance_difference_metadata: dict[str, Any] = field(default_factory=dict)
    evidence_sources: dict[str, Any] = field(default_factory=dict)
    # v0.4 persistence references.  These are metadata only; the canonical
    # scientific report remains the content-addressed artifact on disk.
    reference_selection: InputSelection | None = None
    target_selection: InputSelection | None = None
    input_selections: dict[str, InputSelection] = field(default_factory=dict)
    report_references: tuple[str, ...] = ()
    report_path: str | None = None
    report_hash: str | None = None
    pocket_settings: dict[str, Any] = field(default_factory=dict)
    report_verification: str = "unverified"

    def __post_init__(self) -> None:
        object.__setattr__(self, "target_sources", tuple(self.target_sources))
        object.__setattr__(self, "key_residues", tuple(self.key_residues))
        object.__setattr__(self, "analysis_results", tuple(self.analysis_results))
        object.__setattr__(self, "source_hashes", MappingProxyType(dict(self.source_hashes)))
        object.__setattr__(self, "source_objects", MappingProxyType(dict(self.source_objects)))
        if not isinstance(self.comparison_mode, ComparisonMode):
            object.__setattr__(self, "comparison_mode", ComparisonMode(self.comparison_mode))
        object.__setattr__(self, "visualization_state", MappingProxyType(dict(self.visualization_state)))
        object.__setattr__(self, "msa_settings", MappingProxyType(dict(self.msa_settings)))
        object.__setattr__(self, "interaction_thresholds", MappingProxyType(dict(self.interaction_thresholds)))
        object.__setattr__(self, "site_definitions", tuple(dict(item) for item in self.site_definitions))
        object.__setattr__(
            self, "distance_difference_metadata", MappingProxyType(dict(self.distance_difference_metadata))
        )
        object.__setattr__(self, "evidence_sources", MappingProxyType(dict(self.evidence_sources)))
        if self.reference_selection is not None and not isinstance(self.reference_selection, InputSelection):
            raise TypeError("reference_selection must be an InputSelection or None")
        if self.target_selection is not None and not isinstance(self.target_selection, InputSelection):
            raise TypeError("target_selection must be an InputSelection or None")
        selections = dict(self.input_selections)
        if any(not isinstance(key, str) or not key.strip() for key in selections):
            raise ValueError("input selection names must be non-empty strings")
        if any(not isinstance(value, InputSelection) for value in selections.values()):
            raise TypeError("input_selections must contain InputSelection values")
        object.__setattr__(self, "input_selections", MappingProxyType(selections))
        references = tuple(self.report_references)
        if any(not isinstance(item, str) or not item.strip() for item in references):
            raise ValueError("report_references must contain non-empty strings")
        object.__setattr__(self, "report_references", references)
        if self.report_hash is not None and (
            not isinstance(self.report_hash, str) or not _SHA256.fullmatch(self.report_hash.lower())
        ):
            raise ValueError("report_hash must be a lowercase SHA-256 hexadecimal digest")
        if isinstance(self.report_hash, str):
            object.__setattr__(self, "report_hash", self.report_hash.lower())
        pocket_settings = _json_value(self.pocket_settings)
        if not isinstance(pocket_settings, Mapping):
            raise TypeError("pocket_settings must be a mapping or typed settings object")
        object.__setattr__(self, "pocket_settings", MappingProxyType(dict(pocket_settings)))
        if self.report_verification not in {"unverified", "verified", "legacy_unverified"}:
            raise ValueError("report_verification must be unverified, verified, or legacy_unverified")

    def to_dict(self) -> dict[str, Any]:
        has_v04_state = bool(
            self.reference_selection
            or self.target_selection
            or self.input_selections
            or self.report_references
            or self.report_path
            or self.report_hash
            or self.pocket_settings
        )
        return {
            "schema_version": "4.0" if has_v04_state else self.schema_version,
            "reference_source": self.reference_source,
            "target_sources": list(self.target_sources),
            "settings": {
                "alignment_mode": self.settings.alignment_mode.value,
                "minimum_sequence_identity": self.settings.minimum_sequence_identity,
                "minimum_sequence_coverage": self.settings.minimum_sequence_coverage,
                "substitution_matrix": self.settings.substitution_matrix,
                "gap_open": self.settings.gap_open,
                "gap_extend": self.settings.gap_extend,
                "refined_rmsd": self.settings.refined_rmsd,
                "refinement_cutoff_angstrom": self.settings.refinement_cutoff_angstrom,
                "refinement_max_iterations": self.settings.refinement_max_iterations,
                "usalign_executable": self.settings.usalign_executable,
            },
            "key_residues": [asdict(residue) for residue in self.key_residues],
            "analysis_results": [_analysis_to_dict(result) for result in self.analysis_results],
            "source_hashes": dict(self.source_hashes),
            "source_objects": dict(self.source_objects),
            "visualization_state": dict(self.visualization_state),
            "comparison_mode": self.comparison_mode.value,
            "reference_vs_many": (
                None if self.reference_vs_many is None else _reference_vs_many_to_dict(self.reference_vs_many)
            ),
            "v03": {
                "msa_settings": dict(self.msa_settings),
                "interaction_thresholds": dict(self.interaction_thresholds),
                "site_definitions": [dict(item) for item in self.site_definitions],
                "distance_difference_metadata": dict(self.distance_difference_metadata),
                "evidence_sources": dict(self.evidence_sources),
            },
            "v04": {
                "reference_selection": _selection_to_dict(self.reference_selection),
                "target_selection": _selection_to_dict(self.target_selection),
                "input_selections": {key: _selection_to_dict(value) for key, value in self.input_selections.items()},
                "report_references": list(self.report_references),
                "report_path": self.report_path,
                "report_hash": self.report_hash,
                "pocket_settings": _json_value(self.pocket_settings),
                "report_verification": self.report_verification,
            },
        }

    def to_json(self, *, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent, sort_keys=True, allow_nan=False)

    def save(self, path: str | Path) -> None:
        atomic_write_text(
            path,
            self.to_json() + "\n",
            max_bytes=_MAX_PROJECT_BYTES,
            label="project file",
            idempotent=False,
            replace_existing=True,
        )

    def verify_report(
        self,
        report_path: str | Path | None = None,
        *,
        snapshots: tuple[Any, ...] = (),
        snapshot_dir: str | Path | None = None,
    ) -> bool:
        """Load the canonical report and verify its persisted ID and evidence."""

        if self.report_hash is None:
            raise ProjectSchemaError("project report reference has no persisted report hash; state is unverifiable")
        source = report_path if report_path is not None else self.report_path
        if source is None:
            raise ProjectSchemaError("project has no canonical report path")
        from structlens.application.report_serialization import load_report

        try:
            report = load_report(source, snapshots=snapshots, snapshot_dir=snapshot_dir)
        except (OSError, ValueError, TypeError) as exc:
            raise ProjectSchemaError(f"canonical report verification failed: {exc}") from exc
        if report.report_id != self.report_hash:
            raise ProjectSchemaError("canonical report ID does not match the project report hash")
        return True

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> ProjectState:
        if not isinstance(payload, dict):
            raise ProjectSchemaError("StructLens project payload must be an object")
        try:
            normalized_payload = _json_value(payload)
        except (TypeError, ValueError) as exc:
            raise ProjectSchemaError(f"Invalid StructLens project payload: {exc}") from exc
        if not isinstance(normalized_payload, dict):
            raise ProjectSchemaError("StructLens project payload must be an object")
        payload = normalized_payload
        schema_version = str(payload.get("schema_version", "3.0"))
        if schema_version not in {"0.1", "0.2", "0.3", "3.0", "4.0"}:
            raise ProjectSchemaError("Unsupported StructLens project schema version")
        settings_payload = payload.get("settings", {})
        if not isinstance(settings_payload, dict):
            raise ProjectSchemaError("settings project payload must be an object")
        try:
            settings = AnalysisSettings(
                alignment_mode=AlignmentMode(settings_payload.get("alignment_mode", "auto")),
                minimum_sequence_identity=float(settings_payload.get("minimum_sequence_identity", 0.30)),
                minimum_sequence_coverage=float(settings_payload.get("minimum_sequence_coverage", 0.70)),
                substitution_matrix=str(settings_payload.get("substitution_matrix", "BLOSUM62")),
                gap_open=float(settings_payload.get("gap_open", -10.0)),
                gap_extend=float(settings_payload.get("gap_extend", -0.5)),
                refined_rmsd=bool(settings_payload.get("refined_rmsd", False)),
                refinement_cutoff_angstrom=float(settings_payload.get("refinement_cutoff_angstrom", 2.0)),
                refinement_max_iterations=int(settings_payload.get("refinement_max_iterations", 10)),
                usalign_executable=settings_payload.get("usalign_executable"),
            )
        except (TypeError, ValueError, OverflowError) as exc:
            raise ProjectSchemaError(f"Invalid analysis settings: {exc}") from exc
        residues = tuple(ResidueId(**item) for item in _bounded_items(payload.get("key_residues", []), "key_residues"))
        analyses = tuple(
            _analysis_from_dict(item)
            for item in _bounded_items(payload.get("analysis_results", []), "analysis_results")
        )
        multi_payload = payload.get("reference_vs_many")
        v03 = payload.get("v03", {})
        if not isinstance(v03, dict):
            raise ProjectSchemaError("v03 project payload must be an object")
        v04 = payload.get("v04", {})
        if not isinstance(v04, dict):
            raise ProjectSchemaError("v04 project payload must be an object")
        selection_values = v04.get("input_selections", payload.get("input_selections", {}))
        if not isinstance(selection_values, dict) or len(selection_values) > _MAX_PROJECT_ITEMS:
            raise ProjectSchemaError("input_selections project payload must be an object")
        reference_selection = _selection_from_dict(v04.get("reference_selection", payload.get("reference_selection")))
        target_selection = _selection_from_dict(v04.get("target_selection", payload.get("target_selection")))
        pocket_payload = v04.get("pocket_settings", payload.get("pocket_settings", {}))
        if not isinstance(pocket_payload, dict):
            raise ProjectSchemaError("pocket_settings project payload must be an object")
        stored_version = (
            "4.0"
            if schema_version in {"0.1", "0.2", "0.3"} or (schema_version == "3.0" and "v04" not in payload)
            else schema_version
        )
        report_verification = v04.get("report_verification")
        if report_verification is None:
            report_verification = "legacy_unverified" if schema_version != "4.0" else "unverified"
        return cls(
            reference_source=payload.get("reference_source"),
            target_sources=tuple(
                str(item) for item in _bounded_items(payload.get("target_sources", []), "target_sources")
            ),
            settings=settings,
            key_residues=residues,
            analysis_results=analyses,
            source_hashes=dict(payload.get("source_hashes", {})),
            visualization_state=dict(payload.get("visualization_state", {})),
            schema_version=stored_version,
            source_objects=dict(payload.get("source_objects", {})),
            comparison_mode=ComparisonMode(payload.get("comparison_mode", "pairwise")),
            reference_vs_many=(None if multi_payload is None else _reference_vs_many_from_dict(multi_payload)),
            msa_settings=dict(v03.get("msa_settings", {})),
            interaction_thresholds={
                str(key): float(value) for key, value in dict(v03.get("interaction_thresholds", {})).items()
            },
            site_definitions=tuple(
                dict(item) for item in _bounded_items(v03.get("site_definitions", []), "site_definitions")
            ),
            distance_difference_metadata=dict(v03.get("distance_difference_metadata", {})),
            evidence_sources=dict(v03.get("evidence_sources", {})),
            reference_selection=reference_selection,
            target_selection=target_selection,
            input_selections={
                str(key): _required_selection_from_dict(value) for key, value in selection_values.items()
            },
            report_references=tuple(
                str(item)
                for item in _bounded_items(
                    v04.get("report_references", payload.get("report_references", [])), "report_references"
                )
            ),
            report_path=v04.get("report_path", payload.get("report_path")),
            report_hash=v04.get("report_hash", payload.get("report_hash")),
            pocket_settings=dict(pocket_payload),
            report_verification=str(report_verification),
        )

    @classmethod
    def from_json(cls, payload: str) -> ProjectState:
        try:
            decoded = json.loads(payload)
        except json.JSONDecodeError as exc:
            raise ProjectSchemaError(f"Invalid StructLens project JSON: {exc}") from exc
        if not isinstance(decoded, dict):
            raise ProjectSchemaError("StructLens project JSON must contain an object")
        return cls.from_dict(decoded)

    @classmethod
    def load(cls, path: str | Path) -> ProjectState:
        try:
            data = read_bounded_bytes(path, max_bytes=_MAX_PROJECT_BYTES, label="project file")
            return cls.from_json(data.decode("utf-8"))
        except (OSError, UnicodeDecodeError, ValueError) as exc:
            raise ProjectSchemaError(f"unable to load project file: {exc}") from exc

    def with_analysis(self, analysis: AnalysisResult) -> ProjectState:
        retained = tuple(item for item in self.analysis_results if item.target_id != analysis.target_id)
        return ProjectState(
            reference_source=self.reference_source,
            target_sources=self.target_sources,
            settings=self.settings,
            key_residues=self.key_residues,
            analysis_results=retained + (analysis,),
            source_hashes=dict(self.source_hashes),
            visualization_state=dict(self.visualization_state),
            schema_version=self.schema_version,
            source_objects=dict(self.source_objects),
            comparison_mode=self.comparison_mode,
            reference_vs_many=self.reference_vs_many,
            msa_settings=dict(self.msa_settings),
            interaction_thresholds=dict(self.interaction_thresholds),
            site_definitions=self.site_definitions,
            distance_difference_metadata=dict(self.distance_difference_metadata),
            evidence_sources=dict(self.evidence_sources),
            reference_selection=self.reference_selection,
            target_selection=self.target_selection,
            input_selections=dict(self.input_selections),
            report_references=self.report_references,
            report_path=self.report_path,
            report_hash=self.report_hash,
            pocket_settings=dict(self.pocket_settings),
            report_verification=self.report_verification,
        )

    def with_key_residue(self, residue: ResidueId) -> ProjectState:
        if residue in self.key_residues:
            return self
        return ProjectState(
            reference_source=self.reference_source,
            target_sources=self.target_sources,
            settings=self.settings,
            key_residues=self.key_residues + (residue,),
            analysis_results=self.analysis_results,
            source_hashes=dict(self.source_hashes),
            visualization_state=dict(self.visualization_state),
            schema_version=self.schema_version,
            source_objects=dict(self.source_objects),
            comparison_mode=self.comparison_mode,
            reference_vs_many=self.reference_vs_many,
            msa_settings=dict(self.msa_settings),
            interaction_thresholds=dict(self.interaction_thresholds),
            site_definitions=self.site_definitions,
            distance_difference_metadata=dict(self.distance_difference_metadata),
            evidence_sources=dict(self.evidence_sources),
            reference_selection=self.reference_selection,
            target_selection=self.target_selection,
            input_selections=dict(self.input_selections),
            report_references=self.report_references,
            report_path=self.report_path,
            report_hash=self.report_hash,
            pocket_settings=dict(self.pocket_settings),
            report_verification=self.report_verification,
        )

    def with_source_hashes(self) -> ProjectState:
        paths = tuple(path for path in ((self.reference_source,) + self.target_sources) if path)
        hashes = {path: _sha256_file(path) for path in paths if Path(path).is_file()}
        return ProjectState(
            reference_source=self.reference_source,
            target_sources=self.target_sources,
            settings=self.settings,
            key_residues=self.key_residues,
            analysis_results=self.analysis_results,
            source_hashes=hashes,
            visualization_state=dict(self.visualization_state),
            schema_version=self.schema_version,
            source_objects=dict(self.source_objects),
            comparison_mode=self.comparison_mode,
            reference_vs_many=self.reference_vs_many,
            msa_settings=dict(self.msa_settings),
            interaction_thresholds=dict(self.interaction_thresholds),
            site_definitions=self.site_definitions,
            distance_difference_metadata=dict(self.distance_difference_metadata),
            evidence_sources=dict(self.evidence_sources),
            reference_selection=self.reference_selection,
            target_selection=self.target_selection,
            input_selections=dict(self.input_selections),
            report_references=self.report_references,
            report_path=self.report_path,
            report_hash=self.report_hash,
            pocket_settings=dict(self.pocket_settings),
            report_verification=self.report_verification,
        )

    def migrate(self) -> ProjectState:
        """Return a v0.4 project while retaining all supported legacy data."""

        return ProjectState(
            reference_source=self.reference_source,
            target_sources=self.target_sources,
            settings=self.settings,
            key_residues=self.key_residues,
            analysis_results=self.analysis_results,
            source_hashes=dict(self.source_hashes),
            visualization_state=dict(self.visualization_state),
            schema_version="4.0",
            source_objects=dict(self.source_objects),
            comparison_mode=self.comparison_mode,
            reference_vs_many=self.reference_vs_many,
            msa_settings=dict(self.msa_settings),
            interaction_thresholds=dict(self.interaction_thresholds),
            site_definitions=self.site_definitions,
            distance_difference_metadata=dict(self.distance_difference_metadata),
            evidence_sources=dict(self.evidence_sources),
            reference_selection=self.reference_selection,
            target_selection=self.target_selection,
            input_selections=dict(self.input_selections),
            report_references=self.report_references,
            report_path=self.report_path,
            report_hash=self.report_hash,
            pocket_settings=dict(self.pocket_settings),
            report_verification=self.report_verification,
        )


def _residue_to_dict(residue: ResidueId | None) -> dict[str, Any] | None:
    return None if residue is None else asdict(residue)


def _residue_from_dict(payload: dict[str, Any] | None) -> ResidueId | None:
    return None if payload is None else ResidueId(**payload)


def _analysis_to_dict(result: AnalysisResult) -> dict[str, Any]:
    return {
        "reference_id": result.reference_id,
        "target_id": result.target_id,
        "sequence_identity": result.sequence_identity,
        "sequence_similarity": result.sequence_similarity,
        "sequence_coverage": result.sequence_coverage,
        "alignment_decision": result.alignment_decision,
        "strict_rmsd_angstrom": result.strict_rmsd_angstrom,
        "refined_rmsd_angstrom": result.refined_rmsd_angstrom,
        "mapped_residue_count": result.mapped_residue_count,
        "refined_residue_count": result.refined_residue_count,
        "excluded_alignment_indices": list(result.excluded_alignment_indices),
        "tm_score": result.tm_score,
        "provenance": dict(result.provenance),
        "method_provenance": (None if result.method_provenance is None else result.method_provenance.to_json()),
        "transform": (
            None
            if result.transform is None
            else {
                "rotation": [list(row) for row in result.transform.rotation],
                "translation": list(result.transform.translation),
            }
        ),
        "correspondences": [
            {
                **asdict(item),
                "reference": _residue_to_dict(item.reference),
                "target": _residue_to_dict(item.target),
                "status": item.status.value,
            }
            for item in result.correspondences
        ],
        "mutations": [
            {
                **asdict(item),
                "reference": _residue_to_dict(item.reference),
                "target": _residue_to_dict(item.target),
                "kind": item.kind.value,
            }
            for item in result.mutations
        ],
    }


def _analysis_from_dict(payload: dict[str, Any]) -> AnalysisResult:
    correspondences = tuple(
        ResidueCorrespondence(
            **{
                **item,
                "reference": _residue_from_dict(item.get("reference")),
                "target": _residue_from_dict(item.get("target")),
                "status": CorrespondenceStatus(item["status"]),
            }
        )
        for item in payload.get("correspondences", [])
    )
    mutations = tuple(
        MutationEvent(
            **{
                **item,
                "reference": _residue_from_dict(item.get("reference")),
                "target": _residue_from_dict(item.get("target")),
                "kind": MutationKind(item["kind"]),
            }
        )
        for item in payload.get("mutations", [])
    )
    return AnalysisResult(
        reference_id=payload["reference_id"],
        target_id=payload["target_id"],
        correspondences=correspondences,
        mutations=mutations,
        sequence_identity=float(payload["sequence_identity"]),
        sequence_similarity=(
            None if payload.get("sequence_similarity") is None else float(payload["sequence_similarity"])
        ),
        sequence_coverage=float(payload["sequence_coverage"]),
        alignment_decision=payload["alignment_decision"],
        strict_rmsd_angstrom=payload.get("strict_rmsd_angstrom"),
        refined_rmsd_angstrom=payload.get("refined_rmsd_angstrom"),
        mapped_residue_count=int(payload.get("mapped_residue_count", 0)),
        refined_residue_count=payload.get("refined_residue_count"),
        excluded_alignment_indices=tuple(payload.get("excluded_alignment_indices", [])),
        tm_score=payload.get("tm_score"),
        provenance=payload.get("provenance", {}),
        method_provenance=(
            None
            if payload.get("method_provenance") is None
            else MethodProvenance.from_json(payload["method_provenance"])
        ),
        transform=(
            None
            if payload.get("transform") is None
            else StructuralTransform(
                tuple(tuple(row) for row in payload["transform"]["rotation"]),
                tuple(payload["transform"]["translation"]),
            )
        ),
    )


def _reference_vs_many_to_dict(analysis: ReferenceVsManyAnalysis) -> dict[str, Any]:
    return {
        "reference_id": analysis.reference_id,
        "comparison_mode": analysis.comparison_mode.value,
        "targets": {target_id: _target_analysis_to_dict(target) for target_id, target in analysis.targets.items()},
    }


def _reference_vs_many_from_dict(payload: dict[str, Any]) -> ReferenceVsManyAnalysis:
    return ReferenceVsManyAnalysis(
        reference_id=str(payload["reference_id"]),
        targets={
            str(target_id): _target_analysis_from_dict(target_payload)
            for target_id, target_payload in payload.get("targets", {}).items()
        },
        comparison_mode=ComparisonMode(payload.get("comparison_mode", "reference_vs_many")),
    )


def _target_analysis_to_dict(target: TargetAnalysis) -> dict[str, Any]:
    structural = target.structural_metrics
    return {
        "target_id": target.target_id,
        "correspondence": [_correspondence_to_dict(item) for item in target.correspondence],
        "mutations": [_mutation_to_dict(item) for item in target.mutations],
        "sequence_metrics": {
            "identity": target.sequence_metrics.identity,
            "similarity": target.sequence_metrics.similarity,
            "coverage": target.sequence_metrics.coverage,
            "reference_coverage": target.sequence_metrics.reference_coverage,
            "target_coverage": target.sequence_metrics.target_coverage,
            "aligned_canonical_residue_count": target.sequence_metrics.aligned_canonical_residue_count,
        },
        "structural_metrics": (
            None
            if structural is None
            else {
                "strict_ca_rmsd_angstrom": structural.strict_ca_rmsd_angstrom,
                "mapped_residue_count": structural.mapped_residue_count,
                "rotation": structural.rotation.tolist(),
                "translation": structural.translation.tolist(),
            }
        ),
        "transform": {
            "rotation": [list(row) for row in target.transform.rotation],
            "translation": list(target.transform.translation),
        },
        "provenance": dict(target.provenance),
        "method_provenance": (None if target.method_provenance is None else target.method_provenance.to_json()),
    }


def _target_analysis_from_dict(payload: dict[str, Any]) -> TargetAnalysis:
    sequence_payload = payload["sequence_metrics"]
    structural_payload = payload.get("structural_metrics")
    structural = (
        None
        if structural_payload is None
        else StructuralMetrics(
            float(structural_payload["strict_ca_rmsd_angstrom"]),
            int(structural_payload["mapped_residue_count"]),
            np.asarray(structural_payload["rotation"], dtype=float),
            np.asarray(structural_payload["translation"], dtype=float),
        )
    )
    transform_payload = payload.get("transform", {})
    return TargetAnalysis(
        target_id=str(payload["target_id"]),
        correspondence=tuple(_correspondence_from_dict(item) for item in payload.get("correspondence", [])),
        mutations=tuple(_mutation_from_dict(item) for item in payload.get("mutations", [])),
        sequence_metrics=SequenceAlignmentMetrics(
            float(sequence_payload["identity"]),
            float(sequence_payload["similarity"]),
            float(sequence_payload["coverage"]),
            float(sequence_payload["reference_coverage"]),
            float(sequence_payload["target_coverage"]),
            int(sequence_payload["aligned_canonical_residue_count"]),
        ),
        structural_metrics=structural,
        transform=StructuralTransform(
            tuple(tuple(row) for row in transform_payload.get("rotation", StructuralTransform().rotation)),
            tuple(transform_payload.get("translation", (0.0, 0.0, 0.0))),
        ),
        provenance=payload.get("provenance", {}),
        method_provenance=(
            None
            if payload.get("method_provenance") is None
            else MethodProvenance.from_json(payload["method_provenance"])
        ),
    )


def _correspondence_to_dict(item: ResidueCorrespondence) -> dict[str, Any]:
    payload = asdict(item)
    payload["reference"] = _residue_to_dict(item.reference)
    payload["target"] = _residue_to_dict(item.target)
    payload["status"] = item.status.value
    return payload


def _correspondence_from_dict(payload: dict[str, Any]) -> ResidueCorrespondence:
    return ResidueCorrespondence(
        **{
            **payload,
            "reference": _residue_from_dict(payload.get("reference")),
            "target": _residue_from_dict(payload.get("target")),
            "status": CorrespondenceStatus(payload["status"]),
        }
    )


def _mutation_to_dict(item: MutationEvent) -> dict[str, Any]:
    payload = asdict(item)
    payload["reference"] = _residue_to_dict(item.reference)
    payload["target"] = _residue_to_dict(item.target)
    payload["kind"] = item.kind.value
    return payload


def _mutation_from_dict(payload: dict[str, Any]) -> MutationEvent:
    return MutationEvent(
        **{
            **payload,
            "reference": _residue_from_dict(payload.get("reference")),
            "target": _residue_from_dict(payload.get("target")),
            "kind": MutationKind(payload["kind"]),
        }
    )


def _sha256_file(path: str) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _selection_to_dict(selection: InputSelection | None) -> dict[str, Any] | None:
    if selection is None:
        return None
    return {
        "content_id": selection.content_id,
        "display_name": selection.display_name,
        "format": selection.format.value,
        "model_id": selection.model_id,
        "author_chain_ids": list(selection.author_chain_ids),
        "label_chain_ids": list(selection.label_chain_ids),
        "altloc_policy": selection.altloc_policy.value,
        "assembly_scope": selection.assembly_scope.value,
        "path": selection.path,
        "chain_locators": [
            {
                "author_chain_id": locator.author_chain_id,
                "label_chain_id": locator.label_chain_id,
                "entity_id": locator.entity_id,
            }
            for locator in selection.chain_locators
        ],
        "selection_id": selection.selection_id,
    }


def _bounded_items(value: object, name: str) -> list[Any]:
    if not isinstance(value, (list, tuple)):
        raise ProjectSchemaError(f"{name} project payload must be an array")
    if len(value) > _MAX_PROJECT_ITEMS:
        raise ProjectSchemaError(f"{name} exceeds the maximum item limit")
    return list(value)


def _selection_from_dict(payload: object) -> InputSelection | None:
    if payload is None:
        return None
    if not isinstance(payload, dict):
        raise ProjectSchemaError("InputSelection project payload must be an object")
    try:
        locators = tuple(ChainLocator(**dict(item)) for item in payload.get("chain_locators", []))
        selection = InputSelection(
            content_id=str(payload["content_id"]),
            display_name=str(payload.get("display_name", payload["content_id"])),
            format=StructureFormat(str(payload["format"])),
            model_id=str(payload["model_id"]),
            author_chain_ids=tuple(payload.get("author_chain_ids", [])),
            label_chain_ids=tuple(payload.get("label_chain_ids", [])),
            altloc_policy=AltlocPolicy(str(payload.get("altloc_policy", AltlocPolicy.HIGHEST_OCCUPANCY.value))),
            assembly_scope=AssemblyScope(str(payload.get("assembly_scope", AssemblyScope.ASYMMETRIC_UNIT.value))),
            path=payload.get("path"),
            chain_locators=locators,
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise ProjectSchemaError(f"Invalid InputSelection project payload: {exc}") from exc
    stored_id = payload.get("selection_id")
    if stored_id is not None and stored_id != selection.selection_id:
        raise ProjectSchemaError("InputSelection selection_id does not match its content and settings")
    return selection


def _required_selection_from_dict(payload: object) -> InputSelection:
    selection = _selection_from_dict(payload)
    if selection is None:
        raise ProjectSchemaError("input selection must not be null")
    return selection


def _json_value(value: object, *, depth: int = 0, nodes: list[int] | None = None) -> Any:
    """Copy project metadata while allowing typed settings with ``to_json``."""

    counter = nodes if nodes is not None else [0]
    counter[0] += 1
    if depth > 64 or counter[0] > 1_000_000:
        raise ValueError("project metadata exceeds its nesting or size limit")
    if hasattr(value, "to_json"):
        return _json_value(value.to_json(), depth=depth + 1, nodes=counter)
    if isinstance(value, Mapping):
        if len(value) > _MAX_PROJECT_ITEMS:
            raise ValueError("project metadata contains too many fields")
        result: dict[str, Any] = {}
        for key, item in value.items():
            if not isinstance(key, str) or not key or len(key) > 1_024:
                raise ValueError("project metadata contains an invalid field name")
            result[key] = _json_value(item, depth=depth + 1, nodes=counter)
        return result
    if isinstance(value, (list, tuple)):
        if len(value) > _MAX_PROJECT_ITEMS:
            raise ValueError("project metadata contains too many items")
        return [_json_value(item, depth=depth + 1, nodes=counter) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError("project metadata must contain only finite numbers")
    if isinstance(value, str) and len(value) > 1_000_000:
        raise ValueError("project metadata contains an oversized string")
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    raise TypeError(f"project metadata contains unsupported value {type(value).__name__}")


__all__ = ["ProjectState"]
