"""Versioned, JSON-serializable project state."""

from __future__ import annotations

import hashlib
import json
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


@dataclass(frozen=True, slots=True)
class ProjectState:
    reference_source: str | None = None
    target_sources: tuple[str, ...] = ()
    settings: AnalysisSettings = field(default_factory=AnalysisSettings)
    key_residues: tuple[ResidueId, ...] = ()
    analysis_results: tuple[AnalysisResult, ...] = ()
    source_hashes: dict[str, str] = field(default_factory=dict)
    visualization_state: dict[str, Any] = field(default_factory=dict)
    schema_version: str = "0.2"
    source_objects: dict[str, str] = field(default_factory=dict)
    comparison_mode: ComparisonMode = ComparisonMode.PAIRWISE
    reference_vs_many: ReferenceVsManyAnalysis | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "target_sources", tuple(self.target_sources))
        object.__setattr__(self, "key_residues", tuple(self.key_residues))
        object.__setattr__(self, "analysis_results", tuple(self.analysis_results))
        object.__setattr__(self, "source_hashes", MappingProxyType(dict(self.source_hashes)))
        object.__setattr__(self, "source_objects", MappingProxyType(dict(self.source_objects)))
        if not isinstance(self.comparison_mode, ComparisonMode):
            object.__setattr__(self, "comparison_mode", ComparisonMode(self.comparison_mode))
        object.__setattr__(
            self, "visualization_state", MappingProxyType(dict(self.visualization_state))
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
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
            "analysis_results": [
                _analysis_to_dict(result) for result in self.analysis_results
            ],
            "source_hashes": dict(self.source_hashes),
            "source_objects": dict(self.source_objects),
            "visualization_state": dict(self.visualization_state),
            "comparison_mode": self.comparison_mode.value,
            "reference_vs_many": (
                None
                if self.reference_vs_many is None
                else _reference_vs_many_to_dict(self.reference_vs_many)
            ),
        }

    def to_json(self, *, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent, sort_keys=True)

    def save(self, path: str | Path) -> None:
        Path(path).write_text(self.to_json() + "\n", encoding="utf-8")

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> ProjectState:
        if payload.get("schema_version") not in {"0.1", "0.2"}:
            raise ProjectSchemaError("Unsupported StructLens project schema version")
        settings_payload = payload.get("settings", {})
        settings = AnalysisSettings(
            alignment_mode=AlignmentMode(
                settings_payload.get("alignment_mode", "auto")
            ),
            minimum_sequence_identity=float(
                settings_payload.get("minimum_sequence_identity", 0.30)
            ),
            minimum_sequence_coverage=float(
                settings_payload.get("minimum_sequence_coverage", 0.70)
            ),
            substitution_matrix=str(
                settings_payload.get("substitution_matrix", "BLOSUM62")
            ),
            gap_open=float(settings_payload.get("gap_open", -10.0)),
            gap_extend=float(settings_payload.get("gap_extend", -0.5)),
            refined_rmsd=bool(settings_payload.get("refined_rmsd", False)),
            refinement_cutoff_angstrom=float(
                settings_payload.get("refinement_cutoff_angstrom", 2.0)
            ),
            refinement_max_iterations=int(
                settings_payload.get("refinement_max_iterations", 10)
            ),
            usalign_executable=settings_payload.get("usalign_executable"),
        )
        residues = tuple(ResidueId(**item) for item in payload.get("key_residues", []))
        analyses = tuple(
            _analysis_from_dict(item) for item in payload.get("analysis_results", [])
        )
        multi_payload = payload.get("reference_vs_many")
        return cls(
            reference_source=payload.get("reference_source"),
            target_sources=tuple(payload.get("target_sources", [])),
            settings=settings,
            key_residues=residues,
            analysis_results=analyses,
            source_hashes=dict(payload.get("source_hashes", {})),
            visualization_state=dict(payload.get("visualization_state", {})),
            schema_version=str(payload.get("schema_version", "0.2")),
            source_objects=dict(payload.get("source_objects", {})),
            comparison_mode=ComparisonMode(payload.get("comparison_mode", "pairwise")),
            reference_vs_many=(
                None if multi_payload is None else _reference_vs_many_from_dict(multi_payload)
            ),
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
        return cls.from_json(Path(path).read_text(encoding="utf-8"))

    def with_analysis(self, analysis: AnalysisResult) -> ProjectState:
        retained = tuple(
            item
            for item in self.analysis_results
            if item.target_id != analysis.target_id
        )
        return ProjectState(
            self.reference_source,
            self.target_sources,
            self.settings,
            self.key_residues,
            retained + (analysis,),
            self.source_hashes,
            self.visualization_state,
            self.schema_version,
            self.source_objects,
            self.comparison_mode,
            self.reference_vs_many,
        )

    def with_key_residue(self, residue: ResidueId) -> ProjectState:
        if residue in self.key_residues:
            return self
        return ProjectState(
            self.reference_source,
            self.target_sources,
            self.settings,
            self.key_residues + (residue,),
            self.analysis_results,
            self.source_hashes,
            self.visualization_state,
            self.schema_version,
            self.source_objects,
            self.comparison_mode,
            self.reference_vs_many,
        )

    def with_source_hashes(self) -> ProjectState:
        paths = tuple(
            path
            for path in ((self.reference_source,) + self.target_sources)
            if path
        )
        hashes = {
            path: _sha256_file(path)
            for path in paths
            if Path(path).is_file()
        }
        return ProjectState(
            self.reference_source,
            self.target_sources,
            self.settings,
            self.key_residues,
            self.analysis_results,
            hashes,
            self.visualization_state,
            self.schema_version,
            self.source_objects,
            self.comparison_mode,
            self.reference_vs_many,
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
            None
            if payload.get("sequence_similarity") is None
            else float(payload["sequence_similarity"])
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
        "targets": {
            target_id: _target_analysis_to_dict(target)
            for target_id, target in analysis.targets.items()
        },
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


__all__ = ["ProjectState"]
