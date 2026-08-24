"""Orchestration of mapping, superposition, mutation and residue metrics."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import replace
from pathlib import Path
from threading import Event
from typing import Any, cast

import numpy as np

from structlens.core.alignment.sequence import GlobalSequenceAlignmentEngine
from structlens.core.alignment.superposition import SuperpositionResult, superpose
from structlens.core.errors import AnalysisCancelledError, ChainNotFoundError, MappingError, StructLensError
from structlens.core.geometry.displacement import ca_displacement
from structlens.core.geometry.kabsch import apply_transform
from structlens.core.geometry.rmsd import residue_rmsds
from structlens.core.mapping.manual_mapper import ManualResidueMapper
from structlens.core.mapping.sequence_mapper import SequenceResidueMapper
from structlens.core.metrics.sequence_metrics import calculate_sequence_metrics
from structlens.core.metrics.structural_metrics import StructuralMetrics
from structlens.core.models import (
    AllVsAllAnalysis,
    AnalysisResult,
    AnalysisSettings,
    ComparisonMode,
    MultipleStructureAnalysis,
    MultiStructurePosition,
    PairwiseMatrix,
    ProteinChain,
    ProteinStructure,
    ReferenceVsManyAnalysis,
    ResidueCorrespondence,
    ResidueId,
    ResidueRecord,
    SequenceAlignmentSettings,
    StructuralAlignmentSettings,
    StructuralTransform,
    TargetAnalysis,
)
from structlens.core.mutations.detector import detect_mutations


class AnalysisService:
    """Run a complete pairwise analysis without importing PyMOL."""

    def __init__(self, structural_adapter: object | None = None) -> None:
        self._sequence_engine = GlobalSequenceAlignmentEngine()
        self._mapper = SequenceResidueMapper()
        self._structural_adapter = structural_adapter

    def analyze(
        self,
        reference: ProteinStructure | ProteinChain,
        target: ProteinStructure | ProteinChain,
        settings: AnalysisSettings | None = None,
        *,
        reference_chain_id: str | None = None,
        target_chain_id: str | None = None,
        manual_pairs: list[tuple[object, object]] | None = None,
        cancel_event: Event | None = None,
    ) -> AnalysisResult:
        _check_cancel(cancel_event)
        settings = settings or AnalysisSettings()
        reference_chain = _select_chain(reference, reference_chain_id)
        target_chain = _select_chain(target, target_chain_id)
        if settings.alignment_mode.value == "manual":
            if not manual_pairs:
                raise MappingError("Manual mode requires explicit residue pairs")
            initial = ManualResidueMapper().build_correspondence(
                reference_chain,
                target_chain,
                manual_pairs,  # type: ignore[arg-type]
            )
        else:
            alignment = self._sequence_engine.align(
                reference_chain,
                target_chain,
                settings=_sequence_settings(settings),
            )
            initial = self._mapper.build_correspondence(
                reference_chain, target_chain, alignment
            )
        _check_cancel(cancel_event)
        sequence_metrics = calculate_sequence_metrics(
            initial, reference_chain.sequence, target_chain.sequence
        )
        decision = _alignment_decision(
            settings, sequence_metrics.identity, sequence_metrics.coverage
        )
        tm_score: float | None = None
        structural_result: Any | None = None
        provenance = {
            "mapping_source": "sequence",
            "engine": "Biopython PairwiseAligner",
        }
        if settings.alignment_mode.value == "structure" or (
            settings.alignment_mode.value == "auto"
            and decision.startswith("structure-guided")
        ):
            adapter = self._structural_adapter
            if adapter is None:
                if not reference_chain.source_path or not target_chain.source_path:
                    raise MappingError(
                        "Structure-guided mapping requires source file paths and US-align"
                    )
                from structlens.integrations.usalign.adapter import USAlignAdapter

                adapter = USAlignAdapter(
                    executable=settings.usalign_executable,
                    structure_paths={
                        reference_chain.structure_id: reference_chain.source_path,
                        target_chain.structure_id: target_chain.source_path,
                    },
                )
            structural_result = cast(Any, adapter).align(
                reference_chain,
                target_chain,
                StructuralAlignmentSettings(
                    executable=settings.usalign_executable or "USalign"
                ),
            )
            initial = list(structural_result.correspondences)
            tm_score = structural_result.tm_score
            provenance = {
                "mapping_source": "US-align",
                "engine": "US-align",
                "executable_version": structural_result.executable_version or "unknown",
                **dict(getattr(structural_result, "metadata", {})),
            }
            sequence_metrics = calculate_sequence_metrics(
                initial, reference_chain.sequence, target_chain.sequence
            )

        _check_cancel(cancel_event)
        geometrized, strict, refined, excluded = _calculate_geometry(
            initial, reference_chain, target_chain, settings
        )
        mutations = tuple(detect_mutations(geometrized))
        return AnalysisResult(
            reference_id=reference_chain.structure_id,
            target_id=target_chain.structure_id,
            correspondences=tuple(geometrized),
            mutations=mutations,
            sequence_identity=sequence_metrics.identity,
            sequence_similarity=sequence_metrics.similarity,
            sequence_coverage=sequence_metrics.coverage,
            alignment_decision=decision,
            strict_rmsd_angstrom=strict.strict_rmsd_angstrom if strict else None,
            refined_rmsd_angstrom=refined.strict_rmsd_angstrom if refined else None,
            mapped_residue_count=strict.residue_count if strict else 0,
            refined_residue_count=refined.residue_count if refined else None,
            excluded_alignment_indices=tuple(excluded),
            tm_score=tm_score,
            provenance=provenance,
            transform=_structural_transform(structural_result, provenance, settings),
        )

    def analyze_reference_vs_many(
        self,
        reference: ProteinStructure | ProteinChain,
        targets: Mapping[str, ProteinStructure | ProteinChain]
        | Sequence[ProteinStructure | ProteinChain],
        settings: AnalysisSettings | None = None,
        *,
        reference_chain_id: str | None = None,
        target_chain_ids: Mapping[str, str] | None = None,
        cancel_event: Event | None = None,
    ) -> ReferenceVsManyAnalysis:
        """Analyze every target independently while preserving reference numbering."""

        settings = settings or AnalysisSettings()
        target_items = _target_items(targets)
        results: dict[str, TargetAnalysis] = {}
        reference_id = _select_chain(reference, reference_chain_id).structure_id
        for requested_id, target in target_items:
            _check_cancel(cancel_event)
            target_chain_id = target_chain_ids.get(requested_id) if target_chain_ids else None
            result = self.analyze(
                reference,
                target,
                settings,
                reference_chain_id=reference_chain_id,
                target_chain_id=target_chain_id,
                cancel_event=cancel_event,
            )
            results[requested_id] = _target_analysis_from_result(result, target_id=requested_id)
        return ReferenceVsManyAnalysis(
            reference_id=reference_id,
            targets=results,
            comparison_mode=ComparisonMode.REFERENCE_VS_MANY,
        )

    def analyze_all_vs_all(
        self,
        structures: Mapping[str, ProteinStructure | ProteinChain]
        | Sequence[ProteinStructure | ProteinChain],
        settings: AnalysisSettings | None = None,
        *,
        chain_ids: Mapping[str, str] | None = None,
        cancel_event: Event | None = None,
    ) -> AllVsAllAnalysis:
        """Compute each unordered pair once and expose mirrored matrices."""

        settings = settings or AnalysisSettings()
        items = _target_items(structures)
        structure_ids = tuple(identifier for identifier, _ in items)
        pair_results: dict[tuple[str, str], AnalysisResult | None] = {}
        metric_values: dict[str, dict[tuple[str, str], float | None]] = {
            "sequence_identity": {},
            "sequence_similarity": {},
            "sequence_coverage": {},
            "tm_score": {},
            "strict_ca_rmsd": {},
            "refined_ca_rmsd": {},
        }
        for index, (left_id, left) in enumerate(items):
            for right_id, right in items[index + 1 :]:
                _check_cancel(cancel_event)
                key = (left_id, right_id)
                try:
                    result = self.analyze(
                        left,
                        right,
                        settings,
                        reference_chain_id=chain_ids.get(left_id) if chain_ids else None,
                        target_chain_id=chain_ids.get(right_id) if chain_ids else None,
                        cancel_event=cancel_event,
                    )
                except AnalysisCancelledError:
                    raise
                except StructLensError:
                    pair_results[key] = None
                    continue
                pair_results[key] = result
                metric_values["sequence_identity"][key] = result.sequence_identity
                metric_values["sequence_similarity"][key] = result.sequence_similarity
                metric_values["sequence_coverage"][key] = result.sequence_coverage
                metric_values["tm_score"][key] = result.tm_score
                metric_values["strict_ca_rmsd"][key] = result.strict_rmsd_angstrom
                metric_values["refined_ca_rmsd"][key] = result.refined_rmsd_angstrom
        matrices = {
            name: PairwiseMatrix.from_pairs(
                name,
                structure_ids,
                values,
                unit=("fraction" if name.startswith("sequence_") else "Å")
                if name != "tm_score"
                else "score",
            )
            for name, values in metric_values.items()
        }
        return AllVsAllAnalysis(structure_ids, pair_results, matrices)

    def analyze_multiple_structure_alignment(
        self,
        reference: ProteinStructure | ProteinChain,
        targets: Mapping[str, ProteinStructure | ProteinChain]
        | Sequence[ProteinStructure | ProteinChain],
        settings: AnalysisSettings | None = None,
        *,
        reference_chain_id: str | None = None,
        target_chain_ids: Mapping[str, str] | None = None,
        cancel_event: Event | None = None,
    ) -> MultipleStructureAnalysis:
        """Build a reference-indexed multi-structure summary.

        A dedicated MSTA parser can provide richer transforms later.  The
        deterministic fallback uses the authoritative reference-vs-many rows
        and records that provenance explicitly rather than renumbering input.
        """

        reference_chain = _select_chain(reference, reference_chain_id)
        many = self.analyze_reference_vs_many(
            reference,
            targets,
            settings,
            reference_chain_id=reference_chain_id,
            target_chain_ids=target_chain_ids,
            cancel_event=cancel_event,
        )
        by_index: dict[int, dict[str, ResidueId | None]] = {}
        reference_by_index: dict[int, ResidueId | None] = {}
        deviations: dict[int, dict[str, float | None]] = {}
        for target_id, target in many.targets.items():
            for item in target.correspondence:
                reference_by_index.setdefault(item.alignment_index, item.reference)
                by_index.setdefault(item.alignment_index, {})[target_id] = item.target
                deviations.setdefault(item.alignment_index, {})[target_id] = item.ca_displacement_angstrom
        positions: list[MultiStructurePosition] = []
        for index in sorted(by_index):
            residues = by_index[index]
            mapped = sum(residue is not None for residue in residues.values())
            total = len(many.target_ids) + 1
            variation_values = [value for value in deviations[index].values() if value is not None]
            variability = (
                float(np.std(np.asarray([0.0, *variation_values], dtype=float)))
                if variation_values
                else None
            )
            reference_residue = reference_by_index.get(index)
            positions.append(
                MultiStructurePosition(
                    alignment_index=index,
                    reference_residue=reference_residue,
                    residues={reference_chain.structure_id: reference_residue, **residues},
                    coverage=(mapped + int(reference_residue is not None)) / total,
                    ca_positional_variability_angstrom=variability,
                    per_structure_deviation_angstrom={
                        reference_chain.structure_id: 0.0,
                        **deviations[index],
                    },
                )
            )
        transforms = {
            reference_chain.structure_id: StructuralTransform(),
            **{target_id: target.transform for target_id, target in many.targets.items()},
        }
        return MultipleStructureAnalysis(
            (reference_chain.structure_id, *many.target_ids),
            tuple(positions),
            transforms,
            provenance={"backend": "StructLens reference-indexed fallback"},
        )

    def analyze_paths(
        self,
        reference_path: str | Path,
        target_path: str | Path,
        settings: AnalysisSettings | None = None,
    ) -> AnalysisResult:
        from structlens.core.parsing import load_structure

        reference = load_structure(Path(reference_path))
        target = load_structure(Path(target_path))
        return self.analyze(reference, target, settings)


def _select_chain(
    structure_or_chain: ProteinStructure | ProteinChain,
    chain_id: str | None,
) -> ProteinChain:
    if isinstance(structure_or_chain, ProteinChain):
        return structure_or_chain
    if chain_id is not None:
        for chain in structure_or_chain.chains:
            if chain.chain_id == chain_id:
                return chain
        raise ChainNotFoundError(f"Chain {chain_id!r} was not found")
    if not structure_or_chain.chains:
        raise ChainNotFoundError("Structure has no protein chains")
    return structure_or_chain.chains[0]


def _sequence_settings(settings: AnalysisSettings) -> SequenceAlignmentSettings:
    return SequenceAlignmentSettings(
        settings.substitution_matrix, settings.gap_open, settings.gap_extend
    )


def _check_cancel(cancel_event: Event | None) -> None:
    if cancel_event is not None and cancel_event.is_set():
        raise AnalysisCancelledError("Comparison cancelled by the user")


def _alignment_decision(
    settings: AnalysisSettings, identity: float, coverage: float
) -> str:
    mode = settings.alignment_mode.value
    if mode == "sequence":
        return f"sequence-guided (explicit mode; identity={identity:.3f}, coverage={coverage:.3f})"
    if mode == "manual":
        return "manual (explicit correspondences required)"
    if mode == "structure":
        return "structure-guided (explicit mode)"
    if (
        identity >= settings.minimum_sequence_identity
        and coverage >= settings.minimum_sequence_coverage
    ):
        return f"sequence-guided (AUTO: identity={identity:.3f} >= {settings.minimum_sequence_identity:.3f}, coverage={coverage:.3f} >= {settings.minimum_sequence_coverage:.3f})"
    return f"structure-guided (AUTO: identity={identity:.3f} or coverage={coverage:.3f} below thresholds)"


def _record_map(chain: ProteinChain) -> dict[ResidueId, ResidueRecord]:
    return {record.residue_id: record for record in chain.residue_records}


def _calculate_geometry(
    correspondences: list[ResidueCorrespondence],
    reference: ProteinChain,
    target: ProteinChain,
    settings: AnalysisSettings,
) -> tuple[
    list[ResidueCorrespondence],
    SuperpositionResult | None,
    SuperpositionResult | None,
    list[int],
]:
    reference_records = _record_map(reference)
    target_records = _record_map(target)
    pairs: list[tuple[int, np.ndarray, np.ndarray]] = []
    for item in correspondences:
        if item.reference is None or item.target is None:
            continue
        ref_record = reference_records.get(item.reference)
        target_record = target_records.get(item.target)
        if ref_record is None or target_record is None:
            continue
        ref_ca = next(
            (atom.coordinate for atom in ref_record.atoms if atom.name.upper() == "CA"),
            None,
        )
        target_ca = next(
            (
                atom.coordinate
                for atom in target_record.atoms
                if atom.name.upper() == "CA"
            ),
            None,
        )
        if ref_ca is not None and target_ca is not None:
            pairs.append(
                (item.alignment_index, np.asarray(ref_ca), np.asarray(target_ca))
            )
    if not pairs:
        return correspondences, None, None, []
    ref_coords = np.array([pair[1] for pair in pairs], dtype=float)
    target_coords = np.array([pair[2] for pair in pairs], dtype=float)
    strict = superpose(ref_coords, target_coords, residue_count=len(pairs))
    fitted = apply_transform(target_coords, strict.rotation, strict.translation)
    updated = list(correspondences)
    for pair_index, (alignment_index, ref_coord, _) in enumerate(pairs):
        item = updated[alignment_index]
        if item.reference is None or item.target is None:
            continue
        transformed_target_record = target_records.get(item.target)
        transformed_target_atoms = (
            {
                atom.name: tuple(
                    float(value)
                    for value in apply_transform(
                        np.asarray([atom.coordinate]), strict.rotation, strict.translation
                    )[0]
                )
                for atom in transformed_target_record.atoms
            }
            if transformed_target_record is not None
            else {}
        )
        reference_record = reference_records.get(item.reference)
        residue_metrics = (
            residue_rmsds(
                {atom.name: atom.coordinate for atom in reference_record.atoms},
                transformed_target_atoms,
                reference_record.residue_name,
            )
            if reference_record is not None and transformed_target_record is not None
            else None
        )
        updated[alignment_index] = replace(
            item,
            ca_displacement_angstrom=ca_displacement(ref_coord, fitted[pair_index]),
            backbone_rmsd_angstrom=(
                residue_metrics.backbone_rmsd_angstrom if residue_metrics else None
            ),
            sidechain_rmsd_angstrom=(
                residue_metrics.sidechain_rmsd_angstrom if residue_metrics else None
            ),
            all_heavy_atom_rmsd_angstrom=(
                residue_metrics.all_heavy_atom_rmsd_angstrom if residue_metrics else None
            ),
        )
    excluded: list[int] = []
    refined = strict
    if settings.refined_rmsd and len(pairs) >= 3:
        keep = np.ones(len(pairs), dtype=bool)
        for _ in range(settings.refinement_max_iterations):
            candidate = superpose(
                ref_coords[keep], target_coords[keep], residue_count=int(keep.sum())
            )
            candidate_fitted = apply_transform(
                target_coords, candidate.rotation, candidate.translation
            )
            distances = np.linalg.norm(ref_coords - candidate_fitted, axis=1)
            new_keep = distances <= settings.refinement_cutoff_angstrom
            if new_keep.sum() < 1 or np.array_equal(new_keep, keep):
                refined = candidate
                keep = new_keep
                break
            keep = new_keep
            refined = candidate
        excluded = [pairs[index][0] for index, kept in enumerate(keep) if not kept]
        for alignment_index in excluded:
            updated[alignment_index] = replace(
                updated[alignment_index], is_outlier=True
            )
    return updated, strict, refined if settings.refined_rmsd else None, excluded


def _target_items(
    structures: Mapping[str, ProteinStructure | ProteinChain]
    | Sequence[ProteinStructure | ProteinChain],
) -> tuple[tuple[str, ProteinStructure | ProteinChain], ...]:
    if isinstance(structures, Mapping):
        values = tuple((str(identifier), value) for identifier, value in structures.items())
    else:
        values = tuple((value.structure_id, value) for value in structures)
    identifiers = [identifier for identifier, _ in values]
    if len(set(identifiers)) != len(identifiers):
        raise ValueError("structure identifiers must be unique")
    return values


def _target_analysis_from_result(
    result: AnalysisResult,
    *,
    target_id: str | None = None,
) -> TargetAnalysis:
    structural = None
    if result.strict_rmsd_angstrom is not None:
        structural = StructuralMetrics(
            result.strict_rmsd_angstrom,
            result.mapped_residue_count,
            np.eye(3),
            np.zeros(3),
        )
    from structlens.core.metrics.sequence_metrics import SequenceAlignmentMetrics

    sequence = SequenceAlignmentMetrics(
        result.sequence_identity,
        result.sequence_similarity if result.sequence_similarity is not None else 0.0,
        result.sequence_coverage,
        result.sequence_coverage,
        result.sequence_coverage,
        result.mapped_residue_count,
    )
    return TargetAnalysis(
        target_id=target_id or result.target_id,
        correspondence=result.correspondences,
        mutations=result.mutations,
        sequence_metrics=sequence,
        structural_metrics=structural,
        transform=result.transform or StructuralTransform(),
        provenance=result.provenance,
    )


def _structural_transform(
    structural_result: Any | None,
    provenance: Mapping[str, str],
    settings: AnalysisSettings,
) -> StructuralTransform | None:
    if settings.alignment_mode.value not in {"structure", "auto"}:
        return None
    if provenance.get("mapping_source") != "US-align" or structural_result is None:
        return None
    transform = getattr(structural_result, "transform", None)
    if transform is None:
        return None
    return StructuralTransform(
        tuple(tuple(row) for row in transform.rotation),
        tuple(transform.translation),
    )


__all__ = ["AnalysisService"]
