"""Canonical, report-driven pairwise analysis orchestration.

This service is the application boundary that combines parser-boundary QC,
pairwise analysis, MSA, interactions, site metrics, distance maps,
displacement vectors, and residue Evidence Cards.  Pure section helpers live
in the adjacent ``report_*`` modules so this class remains easy to audit.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import Protocol, TypeVar

from structlens.application.analysis_service import AnalysisService
from structlens.application.dto import AnalysisReportRequest
from structlens.application.interaction_service import InteractionAnalysisService
from structlens.application.msa_service import MuscleAlignmentEngine
from structlens.application.pocket_comparison_service import PocketComparisonService
from structlens.application.pocket_service import PocketDetectionReport, StructurePocketService
from structlens.application.quality_service import StructureQualityService
from structlens.application.report_evidence import (
    evidence_cards as _evidence_cards,
)
from structlens.application.report_evidence import (
    reference_ligand_atoms as _reference_ligand_atoms,
)
from structlens.application.report_evidence import (
    site_metrics as _site_metrics,
)
from structlens.application.report_failures import (
    analysis_failure_report as _analysis_failure_report,
)
from structlens.application.report_failures import (
    invalid_input_report as _invalid_input_report,
)
from structlens.application.report_failures import (
    merge_diagnostics as _merge_diagnostics,
)
from structlens.application.report_failures import (
    section_failure as _section_failure,
)
from structlens.application.report_geometry import (
    analysis_sequence as _analysis_sequence,
)
from structlens.application.report_geometry import (
    displacement_vectors as _displacement_vectors,
)
from structlens.application.report_geometry import (
    distance_map as _distance_map,
)
from structlens.application.report_geometry import (
    position_lookup as _position_lookup,
)
from structlens.application.report_input import (
    apply_selection_diagnostics as _apply_selection_diagnostics,
)
from structlens.application.report_input import (
    canonical_manual_pairs as _canonical_manual_pairs,
)
from structlens.application.report_input import (
    canonical_selection as _canonical_selection,
)
from structlens.application.report_input import (
    canonical_site_definition as _canonical_site_definition,
)
from structlens.application.report_input import (
    canonical_snapshot as _canonical_snapshot,
)
from structlens.application.report_input import (
    inputs_are_usable as _inputs_are_usable,
)
from structlens.application.report_input import (
    selection_diagnostics as _selection_diagnostics,
)
from structlens.application.report_input import (
    single_chain as _single_chain,
)
from structlens.application.report_provenance import report_provenance as _report_provenance
from structlens.application.site_service import define_site
from structlens.core.difference_maps import ResidueDisplacementVector
from structlens.core.evidence import Availability, Diagnostic, DiagnosticSeverity, InteractionEvidence
from structlens.core.interactions import InteractionDifference, InteractionRecord, InteractionThresholds
from structlens.core.models import (
    AnalysisResult,
    AnalysisSettings,
    ProteinChain,
    ProteinStructure,
    ResidueCorrespondence,
    ResidueId,
    ResidueRecord,
)
from structlens.core.msa import MultipleSequenceAlignmentEngine
from structlens.core.parsing import (
    CoordinateQualityError,
    InputSelection,
    ParsedStructure,
    SnapshotError,
    SourceSnapshot,
    StructureParseError,
    load_structure_evidence,
)
from structlens.core.provenance import MethodProvenance
from structlens.core.quality import StructureQualityReport
from structlens.core.reports import AnalysisReport, AnalysisSnapshot, InputQualityBundle, SectionAvailability
from structlens.core.reports.pockets import PocketDetectionSnapshot, PocketReportSnapshot, PocketVolumeSnapshot
from structlens.core.sites import SiteDefinition, SiteMetrics

_SectionValue = TypeVar("_SectionValue")


class _AnalysisRunner(Protocol):
    def analyze(
        self,
        reference: ProteinStructure | ProteinChain,
        target: ProteinStructure | ProteinChain,
        settings: AnalysisSettings | None = None,
        *,
        reference_chain_id: str | None = None,
        target_chain_id: str | None = None,
        manual_pairs: list[tuple[object, object]] | None = None,
    ) -> AnalysisResult: ...


class _QualityRunner(Protocol):
    def analyze(self, parsed: ParsedStructure) -> StructureQualityReport: ...


class _InteractionRunner(Protocol):
    def detect(
        self,
        residues: Sequence[ResidueRecord],
        thresholds: InteractionThresholds | None = None,
        *,
        structure_id: str | None = None,
    ) -> tuple[InteractionRecord, ...]: ...

    def compare(
        self,
        reference: Sequence[InteractionRecord],
        target: Sequence[InteractionRecord],
        position_for: Callable[[ResidueId], str | None],
    ) -> tuple[InteractionDifference, ...]: ...


class ReportService:
    """Build one immutable report while containing failures by section."""

    def __init__(
        self,
        *,
        analysis_service: _AnalysisRunner | None = None,
        quality_service: _QualityRunner | None = None,
        msa_engine: MultipleSequenceAlignmentEngine | None = None,
        interaction_service: _InteractionRunner | None = None,
        pocket_service: StructurePocketService | None = None,
        pocket_comparison_service: PocketComparisonService | None = None,
    ) -> None:
        self._analysis = analysis_service or AnalysisService()
        self._quality = quality_service or StructureQualityService()
        self._msa = msa_engine or MuscleAlignmentEngine()
        self._interactions = interaction_service or InteractionAnalysisService()
        self._pockets = pocket_service or StructurePocketService()
        self._pocket_comparison = pocket_comparison_service or PocketComparisonService()

    def analyze(self, request: AnalysisReportRequest) -> AnalysisReport:
        if not isinstance(request, AnalysisReportRequest):
            raise TypeError("request must be an AnalysisReportRequest")

        reference, reference_qc = self._load_input(
            request.reference_snapshot,
            request.reference_selection,
            source_role="reference",
        )
        target, target_qc = self._load_input(
            request.target_snapshot,
            request.target_selection,
            source_role="target",
        )
        provenance = _report_provenance(request)
        selection_diagnostics = _selection_diagnostics(request, reference, target)
        reference_qc, target_qc = _apply_selection_diagnostics(
            reference_qc,
            target_qc,
            selection_diagnostics,
        )
        input_quality = InputQualityBundle(reference_qc, target_qc)
        input_diagnostics = _merge_diagnostics(
            reference_qc.diagnostics,
            target_qc.diagnostics,
            selection_diagnostics,
        )
        if not _inputs_are_usable(reference, target, reference_qc, target_qc):
            return _invalid_input_report(request, input_quality, provenance, input_diagnostics)

        assert reference is not None and target is not None
        try:
            result = self._analysis.analyze(
                reference.protein_structure,
                target.protein_structure,
                request.analysis_settings,
                reference_chain_id=_single_chain(reference).chain_id,
                target_chain_id=_single_chain(target).chain_id,
                manual_pairs=_canonical_manual_pairs(request.manual_pairs) or None,
            )
        except Exception as error:
            return _analysis_failure_report(
                request,
                input_quality,
                provenance,
                _merge_diagnostics(input_diagnostics, (_section_failure("analysis", error),)),
            )

        return self._downstream_report(
            request,
            reference,
            target,
            input_quality,
            result,
            provenance,
            input_diagnostics,
        )

    def _load_input(
        self,
        snapshot: SourceSnapshot,
        selection: InputSelection,
        *,
        source_role: str,
    ) -> tuple[ParsedStructure | None, StructureQualityReport]:
        canonical_snapshot = _canonical_snapshot(snapshot, source_role)
        canonical_selection = _canonical_selection(selection, source_role)
        try:
            parsed = load_structure_evidence(canonical_snapshot, selection=canonical_selection)
        except CoordinateQualityError as error:
            return None, error.report
        except (SnapshotError, StructureParseError, ValueError):
            diagnostic = Diagnostic(
                code=f"report.input.{source_role}.invalid",
                severity=DiagnosticSeverity.ERROR,
                message=f"The selected {source_role} structure could not be normalized.",
                source_id=canonical_snapshot.display_name,
                remediation="Review the selected model, chain, coordinate format, and parser diagnostics.",
            )
            return None, StructureQualityReport(Availability.INVALID_INPUT, (diagnostic,))
        try:
            quality = self._quality.analyze(parsed)
            if not isinstance(quality, StructureQualityReport):
                raise TypeError("quality runner returned an invalid report")
        except Exception as error:
            diagnostic = Diagnostic(
                code=f"report.input.{source_role}.quality.failed",
                severity=DiagnosticSeverity.ERROR,
                message=(
                    f"The {source_role} input quality check could not be completed "
                    f"({type(error).__name__})."
                ),
                source_id=canonical_snapshot.display_name,
                remediation="Review the input and quality-control configuration before analysis.",
            )
            return parsed, StructureQualityReport(Availability.NUMERICAL_FAILURE, (diagnostic,))
        return parsed, quality

    def _downstream_report(
        self,
        request: AnalysisReportRequest,
        reference: ParsedStructure,
        target: ParsedStructure,
        input_quality: InputQualityBundle,
        result: AnalysisResult,
        provenance: MethodProvenance,
        input_diagnostics: tuple[Diagnostic, ...],
    ) -> AnalysisReport:
        reference_chain = _single_chain(reference)
        target_chain = _single_chain(target)
        site_definitions = tuple(_canonical_site_definition(item) for item in request.site_definitions)
        diagnostics = list(input_diagnostics)
        msa, msa_state = self._run_section(
            "msa",
            diagnostics,
            lambda: self._msa.align(
                (_analysis_sequence(reference_chain), _analysis_sequence(target_chain)),
                request.msa_settings,
            ),
        )
        interaction_evidence, interaction_state = self._run_section(
            "interactions",
            diagnostics,
            lambda: self._interaction_evidence(reference_chain, target_chain, result.correspondences, request),
        )
        sites, site_state = self._site_section(site_definitions, diagnostics, reference, target, result)
        distance_map, map_state = self._run_section(
            "distance_map", diagnostics, lambda: _distance_map(reference_chain, target_chain, result.correspondences)
        )
        vectors, vector_state = self._vector_section(request, diagnostics, reference_chain, target_chain, result)
        pockets, pocket_state = self._pocket_section(
            request,
            reference,
            target,
            result,
            interaction_evidence,
            diagnostics,
            provenance,
        )
        cards, card_state = self._run_section(
            "evidence_cards",
            diagnostics,
            lambda: _evidence_cards(
                result,
                reference_chain,
                msa,
                interaction_evidence,
                sites,
                site_definitions,
                input_quality,
                interaction_state,
                site_state,
            ),
        )
        return AnalysisReport(
            reference_selection=request.reference_selection,
            target_selection=request.target_selection,
            input_quality=input_quality,
            analysis=AnalysisSnapshot.from_result(result),
            msa=msa,
            interactions=interaction_evidence,
            sites=sites,
            distance_map=distance_map,
            displacement_vectors=vectors or (),
            evidence_cards=cards or (),
            site_definitions=site_definitions,
            availability=SectionAvailability(
                input_quality=Availability.AVAILABLE,
                analysis=Availability.AVAILABLE,
                msa=msa_state,
                interactions=interaction_state,
                sites=site_state,
                distance_map=map_state,
                displacement_vectors=vector_state,
                evidence_cards=card_state,
                pockets=pocket_state,
            ),
            diagnostics=_merge_diagnostics(tuple(diagnostics)),
            provenance=provenance,
            pockets=pockets,
        )

    def _pocket_section(
        self,
        request: AnalysisReportRequest,
        reference: ParsedStructure,
        target: ParsedStructure,
        result: AnalysisResult,
        interactions: InteractionEvidence | None,
        diagnostics: list[Diagnostic],
        provenance: MethodProvenance,
    ) -> tuple[PocketReportSnapshot, Availability]:
        """Run pocket detection and comparison as one contained report section."""

        try:
            reference_detection = self._pockets.analyze(reference)
            target_detection = self._pockets.analyze(target)
            reference_volumes = tuple(
                self._pockets.measure_volume(reference, candidate) for candidate in reference_detection.candidates
            )
            target_volumes = tuple(
                self._pockets.measure_volume(target, candidate) for candidate in target_detection.candidates
            )
            reference_volume_map = {item.candidate_id: item for item in reference_volumes if item.candidate_id}
            target_volume_map = {item.candidate_id: item for item in target_volumes if item.candidate_id}
            comparison = self._pocket_comparison.compare(
                reference_detection.candidates,
                target_detection.candidates,
                result.correspondences,
                transform=result.transform,
                reference_volumes=reference_volume_map or None,
                target_volumes=target_volume_map or None,
                mutations=result.mutations,
                interaction_differences=(() if interactions is None else interactions.differences),
            )
            reference_snapshot = _pocket_detection_snapshot("reference", reference_detection, self._pockets)
            target_snapshot = _pocket_detection_snapshot("target", target_detection, self._pockets)
            volume_snapshots = tuple(
                PocketVolumeSnapshot("reference", item) for item in reference_volumes
            ) + tuple(PocketVolumeSnapshot("target", item) for item in target_volumes)
            section_diagnostics = (
                tuple(reference_detection.diagnostics)
                + tuple(target_detection.diagnostics)
                + tuple(diagnostic for item in reference_volumes for diagnostic in item.diagnostics)
                + tuple(diagnostic for item in target_volumes for diagnostic in item.diagnostics)
                + tuple(comparison.diagnostics)
            )
            # Keep the report-level diagnostic stream complete without relying
            # on a GUI-side reconstruction of pocket failures.
            diagnostics.extend(item for item in section_diagnostics if item not in diagnostics)
            states = (
                reference_detection.availability,
                target_detection.availability,
                *(item.availability for item in reference_volumes),
                *(item.availability for item in target_volumes),
                comparison.availability,
            )
            availability = _aggregate_availability(states)
            snapshot = PocketReportSnapshot(
                detections=(reference_snapshot, target_snapshot),
                volumes=volume_snapshots,
                matching=comparison.matching,
                comparisons=comparison.comparisons,
                concordance=comparison.concordance,
                diagnostics=section_diagnostics,
                provenance=provenance,
                units={"length": "angstrom", "volume": "angstrom^3"},
                availability=availability,
            )
            return snapshot, availability
        except Exception as error:
            diagnostic = _section_failure("pockets", error)
            diagnostics.append(diagnostic)
            snapshot = PocketReportSnapshot(
                diagnostics=(diagnostic,),
                provenance=provenance,
                units={"length": "angstrom", "volume": "angstrom^3"},
                availability=Availability.NUMERICAL_FAILURE,
            )
            return snapshot, Availability.NUMERICAL_FAILURE

    def _interaction_evidence(
        self,
        reference: ProteinChain,
        target: ProteinChain,
        correspondences: Sequence[ResidueCorrespondence],
        request: AnalysisReportRequest,
    ) -> InteractionEvidence:
        reference_records = self._interactions.detect(
            reference.residue_records,
            request.interaction_thresholds,
            structure_id=reference.structure_id,
        )
        target_records = self._interactions.detect(
            target.residue_records,
            request.interaction_thresholds,
            structure_id=target.structure_id,
        )
        positions = _position_lookup(correspondences)
        differences = self._interactions.compare(reference_records, target_records, positions.get)
        return InteractionEvidence(differences, reference_records, target_records)

    def _site_section(
        self,
        definitions: tuple[SiteDefinition, ...],
        diagnostics: list[Diagnostic],
        reference: ParsedStructure,
        target: ParsedStructure,
        result: AnalysisResult,
    ) -> tuple[tuple[SiteMetrics, ...], Availability]:
        if not definitions:
            return (), Availability.NOT_APPLICABLE
        ligand_atoms = _reference_ligand_atoms(reference)
        reference_records = _single_chain(reference).residue_records
        resolved: list[SiteDefinition] = []
        for definition in definitions:
            if define_site(definition, reference_records, ligand_atoms=ligand_atoms):
                resolved.append(definition)
            else:
                diagnostics.append(
                    Diagnostic(
                        code=f"report.sites.{definition.site_id}.not_detected",
                        severity=DiagnosticSeverity.WARNING,
                        message=f"Site definition {definition.site_id!r} resolved to no reference residues.",
                        remediation="Review the residue identity, ligand identity, center, and search radius.",
                    )
                )
        if not resolved:
            return (), Availability.NOT_DETECTED
        values, state = self._run_section(
            "sites", diagnostics, lambda: _site_metrics(tuple(resolved), reference, target, result)
        )
        return values or (), state

    def _vector_section(
        self,
        request: AnalysisReportRequest,
        diagnostics: list[Diagnostic],
        reference: ProteinChain,
        target: ProteinChain,
        result: AnalysisResult,
    ) -> tuple[tuple[ResidueDisplacementVector, ...], Availability]:
        if result.transform is None:
            diagnostics.append(
                Diagnostic(
                    "report.vectors.transform_unavailable",
                    DiagnosticSeverity.WARNING,
                    "Displacement vectors require a target-to-reference rigid transform.",
                    remediation="Use mapped residues with C-alpha coordinates to obtain a rigid fit.",
                )
            )
            return (), Availability.NUMERICAL_FAILURE
        values, state = self._run_section(
            "vectors", diagnostics, lambda: _displacement_vectors(reference, target, result, request)
        )
        return values or (), state

    @staticmethod
    def _run_section(
        name: str,
        diagnostics: list[Diagnostic],
        operation: Callable[[], _SectionValue],
    ) -> tuple[_SectionValue | None, Availability]:
        try:
            return operation(), Availability.AVAILABLE
        except Exception as error:
            diagnostics.append(_section_failure(name, error))
            return None, Availability.NUMERICAL_FAILURE


def _pocket_detection_snapshot(
    role: str,
    report: PocketDetectionReport,
    service: StructurePocketService,
) -> PocketDetectionSnapshot:
    """Adapt the application detector result to the report's typed contract."""

    return PocketDetectionSnapshot(
        role=role,
        availability=report.availability,
        candidates=report.candidates,
        settings=service.settings,
        diagnostics=report.diagnostics,
        counts=report.counts,
        provenance=report.provenance,
    )


def _aggregate_availability(statuses: tuple[Availability, ...]) -> Availability:
    """Preserve fatal pocket evidence states while retaining no-result states."""

    values = tuple(statuses)
    for fatal in (Availability.INVALID_INPUT, Availability.DEPENDENCY_UNAVAILABLE, Availability.NUMERICAL_FAILURE):
        if fatal in values:
            return fatal
    if Availability.NOT_DETECTED in values:
        return Availability.NOT_DETECTED
    if Availability.NOT_APPLICABLE in values:
        return Availability.NOT_APPLICABLE
    return Availability.AVAILABLE


__all__ = ["ReportService"]
