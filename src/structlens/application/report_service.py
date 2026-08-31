"""Canonical, report-driven pairwise analysis orchestration.

This service is the only application path that combines parser-boundary QC,
pairwise analysis, MSA, interactions, site metrics, distance maps,
displacement vectors, and residue Evidence Cards.  Every calculation consumes
the same captured source snapshots and selected normalized structures.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import fields
from importlib.metadata import PackageNotFoundError, version
from typing import Protocol, TypeVar

import numpy as np

from structlens.application.analysis_service import AnalysisService
from structlens.application.difference_map_service import (
    build_displacement_vectors,
    calculate_distance_difference,
)
from structlens.application.dto import AnalysisReportRequest
from structlens.application.interaction_service import InteractionAnalysisService
from structlens.application.msa_service import MuscleAlignmentEngine
from structlens.application.quality_service import StructureQualityService
from structlens.application.site_service import calculate_site_metrics, define_site
from structlens.core.difference_maps import ResidueDisplacementVector
from structlens.core.evidence import (
    Availability,
    Diagnostic,
    DiagnosticSeverity,
    EvidenceCard,
    InteractionEvidence,
    SequenceEvidence,
    SiteEvidence,
    StructureEvidence,
    build_evidence_card,
    quality_for_sections,
)
from structlens.core.geometry.kabsch import apply_transform
from structlens.core.interactions import InteractionDifference, InteractionRecord, InteractionThresholds
from structlens.core.models import (
    AnalysisResult,
    AnalysisSettings,
    AtomRecord,
    ProteinChain,
    ProteinStructure,
    ResidueCorrespondence,
    ResidueId,
    ResidueRecord,
    StructuralTransform,
)
from structlens.core.msa import (
    AnalysisSequence,
    MultipleSequenceAlignment,
    MultipleSequenceAlignmentEngine,
    SequenceResidueRef,
)
from structlens.core.parsing import (
    CoordinateQualityError,
    InputSelection,
    ParsedStructure,
    SnapshotError,
    SourceSnapshot,
    StructureParseError,
    load_structure_evidence,
)
from structlens.core.provenance import FrozenJSON, MethodProvenance
from structlens.core.quality import StructureQualityReport
from structlens.core.reports import (
    AnalysisReport,
    AnalysisSnapshot,
    DistanceMapSnapshot,
    InputQualityBundle,
    SectionAvailability,
)
from structlens.core.sites import SiteDefinition, SiteMetrics

_METHOD_ID = "structlens.analysis_report"
_METHOD_VERSION = "1"
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
    ) -> None:
        self._analysis = analysis_service or AnalysisService()
        self._quality = quality_service or StructureQualityService()
        self._msa = msa_engine or MuscleAlignmentEngine()
        self._interactions = interaction_service or InteractionAnalysisService()

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
        input_quality = InputQualityBundle(reference_qc, target_qc)
        provenance = _report_provenance(request)
        input_diagnostics = _merge_diagnostics(
            reference_qc.diagnostics,
            target_qc.diagnostics,
            _selection_diagnostics(request, reference, target),
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
            diagnostic = _section_failure("analysis", error)
            return _analysis_failure_report(
                request,
                input_quality,
                provenance,
                _merge_diagnostics(input_diagnostics, (diagnostic,)),
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
        except (SnapshotError, StructureParseError, ValueError) as error:
            diagnostic = Diagnostic(
                code=f"report.input.{source_role}.invalid",
                severity=DiagnosticSeverity.ERROR,
                message=f"The selected {source_role} structure could not be normalized.",
                source_id=canonical_snapshot.display_name,
                remediation="Review the selected model, chain, coordinate format, and parser diagnostics.",
            )
            del error
            return None, StructureQualityReport(Availability.INVALID_INPUT, (diagnostic,))
        quality = self._quality.analyze(parsed)
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
            lambda: self._interaction_evidence(
                reference_chain,
                target_chain,
                result.correspondences,
                request,
            ),
        )
        sites, site_state = self._site_section(
            site_definitions,
            diagnostics,
            reference,
            target,
            result,
        )
        distance_map, map_state = self._run_section(
            "distance_map",
            diagnostics,
            lambda: _distance_map(reference_chain, target_chain, result.correspondences),
        )
        vectors, vector_state = self._vector_section(
            request,
            diagnostics,
            reference_chain,
            target_chain,
            result,
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
            ),
            diagnostics=_merge_diagnostics(tuple(diagnostics)),
            provenance=provenance,
        )

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
            "sites",
            diagnostics,
            lambda: _site_metrics(tuple(resolved), reference, target, result),
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
            "vectors",
            diagnostics,
            lambda: _displacement_vectors(reference, target, result, request),
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


def _inputs_are_usable(
    reference: ParsedStructure | None,
    target: ParsedStructure | None,
    reference_qc: StructureQualityReport,
    target_qc: StructureQualityReport,
) -> bool:
    return (
        reference is not None
        and target is not None
        and reference_qc.availability is Availability.AVAILABLE
        and target_qc.availability is Availability.AVAILABLE
        and len(reference.protein_structure.chains) == 1
        and len(target.protein_structure.chains) == 1
    )


def _selection_diagnostics(
    request: AnalysisReportRequest,
    reference: ParsedStructure | None,
    target: ParsedStructure | None,
) -> tuple[Diagnostic, ...]:
    diagnostics: list[Diagnostic] = []
    for role, parsed, selection in (
        ("reference", reference, request.reference_selection),
        ("target", target, request.target_selection),
    ):
        if parsed is not None and len(parsed.protein_structure.chains) != 1:
            diagnostics.append(
                Diagnostic(
                    code=f"report.input.{role}.chain_count",
                    severity=DiagnosticSeverity.ERROR,
                    message=(
                        f"The {role} selection resolved to {len(parsed.protein_structure.chains)} protein chains; "
                        "pairwise reporting requires exactly one."
                    ),
                    source_id=f"{role}.{selection.format.value}",
                    remediation="Select exactly one author or label chain for this comparison.",
                )
            )
    return tuple(diagnostics)


def _single_chain(parsed: ParsedStructure) -> ProteinChain:
    if len(parsed.protein_structure.chains) != 1:
        raise ValueError("the report workflow requires exactly one selected chain per structure")
    return parsed.protein_structure.chains[0]


def _canonical_snapshot(snapshot: SourceSnapshot, role: str) -> SourceSnapshot:
    suffix = "pdb" if snapshot.logical_format == "pdb" else "cif"
    return SourceSnapshot(
        decompressed_bytes=snapshot.decompressed_bytes,
        raw_sha256=snapshot.raw_sha256,
        content_id=snapshot.content_id,
        decompressed_sha256=snapshot.decompressed_sha256,
        display_name=f"{role}.{suffix}",
        logical_format=snapshot.logical_format,
        is_gzip=snapshot.is_gzip,
    )


def _canonical_selection(selection: InputSelection, role: str) -> InputSelection:
    suffix = "pdb" if selection.format.value == "pdb" else "cif"
    return InputSelection(
        selection.content_id,
        f"{role}.{suffix}",
        selection.format,
        selection.model_id,
        author_chain_ids=selection.author_chain_ids,
        label_chain_ids=selection.label_chain_ids,
        altloc_policy=selection.altloc_policy,
        assembly_scope=selection.assembly_scope,
        chain_locators=selection.chain_locators,
    )


def _canonical_residue_id(residue: ResidueId, role: str = "reference") -> ResidueId:
    return ResidueId(
        role,
        residue.model_id,
        residue.chain_id,
        residue.auth_seq_id,
        residue.insertion_code,
        residue.residue_name,
    )


def _canonical_optional_residue_id(
    residue: ResidueId | None,
    role: str = "reference",
) -> ResidueId | None:
    return None if residue is None else _canonical_residue_id(residue, role)


def _canonical_site_definition(definition: SiteDefinition) -> SiteDefinition:
    return SiteDefinition(
        definition.site_id,
        definition.name,
        definition.mode,
        tuple(_canonical_residue_id(item) for item in definition.reference_residues),
        _canonical_optional_residue_id(definition.center_residue),
        definition.ligand_id,
        definition.radius_angstrom,
    )


def _canonical_manual_pairs(
    pairs: Sequence[tuple[ResidueId, ResidueId]],
) -> list[tuple[object, object]]:
    return [
        (_canonical_residue_id(reference, "reference"), _canonical_residue_id(target, "target"))
        for reference, target in pairs
    ]


def _analysis_sequence(chain: ProteinChain) -> AnalysisSequence:
    residues = tuple(
        SequenceResidueRef(index, record.one_letter or "X", record.residue_id)
        for index, record in enumerate(chain.residue_records)
    )
    return AnalysisSequence(chain.structure_id, chain.chain_id, chain.sequence, residues, "structure")


def _reference_position(residue: ResidueId) -> str:
    insertion = residue.insertion_code or ""
    return f"{residue.chain_id}:{residue.auth_seq_id}{insertion}"


def _position_lookup(correspondences: Sequence[ResidueCorrespondence]) -> dict[ResidueId, str]:
    output: dict[ResidueId, str] = {}
    for item in correspondences:
        if item.reference is None:
            continue
        position = _reference_position(item.reference)
        output[item.reference] = position
        if item.target is not None:
            output[item.target] = position
    return output


def _ca(record: ResidueRecord) -> tuple[float, float, float] | None:
    coordinate = next((atom.coordinate for atom in record.atoms if atom.name.upper() == "CA"), None)
    if coordinate is None:
        return None
    return (float(coordinate[0]), float(coordinate[1]), float(coordinate[2]))


def _coordinate_maps(
    reference: ProteinChain,
    target: ProteinChain,
    correspondences: Sequence[ResidueCorrespondence],
) -> tuple[
    tuple[str, ...],
    dict[str, ResidueId],
    dict[str, ResidueId],
    dict[str, tuple[float, float, float]],
    dict[str, tuple[float, float, float]],
]:
    reference_by_id = {record.residue_id: record for record in reference.residue_records}
    target_by_id = {record.residue_id: record for record in target.residue_records}
    labels: list[str] = []
    reference_ids: dict[str, ResidueId] = {}
    target_ids: dict[str, ResidueId] = {}
    reference_ca: dict[str, tuple[float, float, float]] = {}
    target_ca: dict[str, tuple[float, float, float]] = {}
    for item in correspondences:
        if item.reference is None:
            continue
        label = _reference_position(item.reference)
        labels.append(label)
        reference_ids[label] = item.reference
        reference_record = reference_by_id.get(item.reference)
        coordinate = _ca(reference_record) if reference_record is not None else None
        if coordinate is not None:
            reference_ca[label] = coordinate
        if item.target is None:
            continue
        target_ids[label] = item.target
        target_record = target_by_id.get(item.target)
        coordinate = _ca(target_record) if target_record is not None else None
        if coordinate is not None:
            target_ca[label] = coordinate
    return tuple(labels), reference_ids, target_ids, reference_ca, target_ca


def _distance_map(
    reference: ProteinChain,
    target: ProteinChain,
    correspondences: Sequence[ResidueCorrespondence],
) -> DistanceMapSnapshot:
    labels, _, _, reference_ca, target_ca = _coordinate_maps(reference, target, correspondences)
    matrix = calculate_distance_difference(labels, reference_ca, target_ca)
    return DistanceMapSnapshot.from_matrix(matrix)


def _displacement_vectors(
    reference: ProteinChain,
    target: ProteinChain,
    result: AnalysisResult,
    request: AnalysisReportRequest,
) -> tuple[ResidueDisplacementVector, ...]:
    assert result.transform is not None
    labels, reference_ids, target_ids, reference_ca, target_ca = _coordinate_maps(
        reference,
        target,
        result.correspondences,
    )
    transformed = _transform_coordinates(target_ca, result.transform)
    return build_displacement_vectors(
        labels,
        reference_ids,
        target_ids,
        reference_ca,
        transformed,
        minimum_magnitude_angstrom=request.minimum_vector_magnitude_angstrom,
        maximum_vectors=request.maximum_vectors,
    )


def _transform_coordinates(
    coordinates: Mapping[str, Sequence[float]],
    transform: StructuralTransform,
) -> dict[str, tuple[float, float, float]]:
    if not coordinates:
        return {}
    labels = tuple(coordinates)
    values = np.asarray([coordinates[label] for label in labels], dtype=np.float64)
    fitted = apply_transform(values, transform.rotation, transform.translation)
    return {
        label: (float(row[0]), float(row[1]), float(row[2]))
        for label, row in zip(labels, fitted, strict=True)
    }


def _site_metrics(
    definitions: Sequence[SiteDefinition],
    reference: ParsedStructure,
    target: ParsedStructure,
    result: AnalysisResult,
) -> tuple[SiteMetrics, ...]:
    reference_chain = _single_chain(reference)
    target_chain = _single_chain(target)
    correspondence = {
        item.reference: item.target
        for item in result.correspondences
        if item.reference is not None and item.target is not None
    }
    matrix = _homogeneous_transform(result.transform)
    ligand_atoms = _reference_ligand_atoms(reference)
    return tuple(
        calculate_site_metrics(
            definition,
            reference_chain.residue_records,
            target_chain.residue_records,
            correspondence,
            target_structure_id=target_chain.structure_id,
            target_transform=matrix,
            ligand_atoms=ligand_atoms,
        )
        for definition in definitions
    )


def _reference_ligand_atoms(reference: ParsedStructure) -> dict[str, Sequence[AtomRecord]]:
    ligand_atoms: dict[str, Sequence[AtomRecord]] = {
        component.component_id: component.atoms
        for component in reference.components
        if component.is_ligand
    }
    ligand_atoms.update(
        {
            component.residue_name: component.atoms
            for component in reference.components
            if component.is_ligand and component.residue_name is not None
        }
    )
    return ligand_atoms


def _homogeneous_transform(transform: StructuralTransform | None) -> np.ndarray | None:
    if transform is None:
        return None
    matrix = np.eye(4, dtype=np.float64)
    matrix[:3, :3] = np.asarray(transform.rotation, dtype=np.float64).T
    matrix[:3, 3] = np.asarray(transform.translation, dtype=np.float64)
    return matrix


def _evidence_cards(
    result: AnalysisResult,
    reference: ProteinChain,
    msa: object | None,
    interactions: object | None,
    sites: object | None,
    definitions: Sequence[SiteDefinition],
    quality: InputQualityBundle,
    interaction_state: Availability,
    site_state: Availability,
) -> tuple[EvidenceCard, ...]:
    typed_msa = msa if isinstance(msa, MultipleSequenceAlignment) else None
    typed_interactions = interactions if isinstance(interactions, InteractionEvidence) else InteractionEvidence()
    typed_sites = tuple(sites) if isinstance(sites, tuple) and all(isinstance(item, SiteMetrics) for item in sites) else ()
    sequence_refs = {item.residue_id: item for item in _analysis_sequence(reference).residues}
    columns = {
        column.reference_residue: column
        for column in typed_msa.columns
        if column.reference_residue is not None
    } if typed_msa is not None else {}
    position_for = _position_lookup(result.correspondences)
    warnings = tuple(
        diagnostic.code
        for diagnostic in quality.reference.diagnostics + quality.target.diagnostics
        if diagnostic.severity is DiagnosticSeverity.WARNING
    )
    cards: list[EvidenceCard] = []
    for item in result.correspondences:
        if item.reference is None:
            continue
        residue_ref = sequence_refs.get(item.reference)
        if residue_ref is None:
            residue_ref = SequenceResidueRef(item.alignment_index, item.reference_one_letter or "X", item.reference)
        column = columns.get(item.reference)
        sequence = SequenceEvidence(
            reference_one_letter=item.reference_one_letter,
            target_one_letter=item.target_one_letter,
            alignment_index=item.alignment_index,
            sequence_identity=result.sequence_identity,
            conservation_fraction=column.conservation_score if column is not None else None,
            entropy_bits=column.entropy_bits if column is not None else None,
            gap_fraction=column.gap_fraction if column is not None else None,
            ambiguous_fraction=column.ambiguous_fraction if column is not None else None,
            source_refs=tuple(
                cell.residue for cell in column.cells if cell.residue is not None
            ) if column is not None else (residue_ref,),
        )
        position = _reference_position(item.reference)
        interaction_evidence = _interactions_for_position(typed_interactions, position, position_for)
        site_values = _sites_for_residue(item.reference, typed_sites, definitions, reference.residue_records)
        available = ["sequence", "structure"]
        unavailable: list[str] = []
        if interaction_state is Availability.AVAILABLE:
            available.append("interactions")
        else:
            unavailable.append("interactions")
        if site_state is Availability.AVAILABLE:
            available.append("sites")
        else:
            unavailable.append("sites")
        cards.append(
            build_evidence_card(
                residue_ref,
                target_id=result.target_id,
                sequence=sequence,
                structure=StructureEvidence(
                    item.ca_displacement_angstrom,
                    item.backbone_rmsd_angstrom,
                    item.sidechain_rmsd_angstrom,
                    item.all_heavy_atom_rmsd_angstrom,
                    available=item.target is not None and item.ca_displacement_angstrom is not None,
                ),
                interactions=interaction_evidence,
                site=SiteEvidence(site_values),
                quality=quality_for_sections(available, unavailable, warnings=warnings),
                provenance=("structlens.analysis_report", result.alignment_decision),
            )
        )
    return tuple(cards)


def _interactions_for_position(
    evidence: InteractionEvidence,
    position: str,
    position_for: Mapping[ResidueId, str],
) -> InteractionEvidence:
    differences = tuple(
        item
        for item in evidence.differences
        if position in {item.key.reference_position_a, item.key.reference_position_b}
    )
    reference = tuple(
        item
        for item in evidence.reference_interactions
        if position in {position_for.get(item.residue_a), position_for.get(item.residue_b) if item.residue_b else None}
    )
    target = tuple(
        item
        for item in evidence.target_interactions
        if position in {position_for.get(item.residue_a), position_for.get(item.residue_b) if item.residue_b else None}
    )
    return InteractionEvidence(differences, reference, target)


def _sites_for_residue(
    residue: ResidueId,
    metrics: Sequence[SiteMetrics],
    definitions: Sequence[SiteDefinition],
    reference_records: Sequence[ResidueRecord],
) -> tuple[SiteMetrics, ...]:
    relevant = {
        definition.site_id
        for definition in definitions
        if residue in {record.residue_id for record in define_site(definition, reference_records)}
    }
    return tuple(metric for metric in metrics if metric.site_id in relevant)


def _report_provenance(request: AnalysisReportRequest) -> MethodProvenance:
    analysis = request.analysis_settings
    thresholds = request.interaction_thresholds
    parameters: dict[str, FrozenJSON] = {
        "transform_direction": "target_to_reference",
        "coordinate_convention": "row_vector_target_at_rotation_plus_translation",
        "reference_selection_id": request.reference_selection.selection_id,
        "target_selection_id": request.target_selection.selection_id,
        "alignment_mode": analysis.alignment_mode.value,
        "minimum_sequence_identity": analysis.minimum_sequence_identity,
        "minimum_sequence_coverage": analysis.minimum_sequence_coverage,
        "substitution_matrix": analysis.substitution_matrix,
        "gap_open": analysis.gap_open,
        "gap_extend": analysis.gap_extend,
        "refined_rmsd": analysis.refined_rmsd,
        "refinement_cutoff_angstrom": analysis.refinement_cutoff_angstrom,
        "refinement_max_iterations": analysis.refinement_max_iterations,
        "msa_algorithm": request.msa_settings.algorithm,
        "msa_mode": request.msa_settings.mode,
        "site_count": len(request.site_definitions),
        "site_definitions": tuple(_site_definition_parameter(item) for item in request.site_definitions),
        "manual_pairs": tuple(
            {
                "reference": _residue_parameter(_canonical_residue_id(reference, "reference")),
                "target": _residue_parameter(_canonical_residue_id(target, "target")),
            }
            for reference, target in request.manual_pairs
        ),
        "interaction_scope": "selected_polymer_residue_heavy_atoms",
        "minimum_vector_magnitude_angstrom": request.minimum_vector_magnitude_angstrom,
        "maximum_vectors": request.maximum_vectors,
    }
    parameters.update(
        {f"interaction_{field.name}": float(getattr(thresholds, field.name)) for field in fields(thresholds)}
    )
    units = {
        "refinement_cutoff_angstrom": "angstrom",
        "minimum_vector_magnitude_angstrom": "angstrom",
        "site_definitions.radius_angstrom": "angstrom",
        **{
            f"interaction_{field.name}": "degree" if field.name.endswith("_degrees") else "angstrom"
            for field in fields(thresholds)
        },
    }
    return MethodProvenance(
        _METHOD_ID,
        _METHOD_VERSION,
        parameters=parameters,
        units=units,
        backend_versions={
            "biopython": _package_version("biopython"),
            "numpy": _package_version("numpy"),
            "scipy": _package_version("scipy"),
        },
        input_hashes={
            "reference_raw": request.reference_snapshot.raw_sha256,
            "reference_content": request.reference_snapshot.content_id,
            "target_raw": request.target_snapshot.raw_sha256,
            "target_content": request.target_snapshot.content_id,
        },
        analyzed_representation="asymmetric_unit_pair",
    )


def _site_definition_parameter(definition: SiteDefinition) -> Mapping[str, FrozenJSON]:
    return {
        "site_id": definition.site_id,
        "name": definition.name,
        "mode": definition.mode.value,
        "reference_residues": tuple(
            _residue_parameter(_canonical_residue_id(item)) for item in definition.reference_residues
        ),
        "center_residue": _residue_parameter(_canonical_optional_residue_id(definition.center_residue)),
        "ligand_id": definition.ligand_id,
        "radius_angstrom": definition.radius_angstrom,
    }


def _residue_parameter(residue: ResidueId | None) -> Mapping[str, FrozenJSON] | None:
    if residue is None:
        return None
    return {
        "structure_id": residue.structure_id,
        "model_id": residue.model_id,
        "chain_id": residue.chain_id,
        "auth_seq_id": residue.auth_seq_id,
        "insertion_code": residue.insertion_code,
        "residue_name": residue.residue_name,
    }


def _package_version(package: str) -> str:
    try:
        return version(package)
    except PackageNotFoundError:
        return "unavailable"


def _invalid_input_report(
    request: AnalysisReportRequest,
    input_quality: InputQualityBundle,
    provenance: MethodProvenance,
    diagnostics: tuple[Diagnostic, ...],
) -> AnalysisReport:
    state = SectionAvailability(
        input_quality=Availability.INVALID_INPUT,
        analysis=Availability.INVALID_INPUT,
        msa=Availability.INVALID_INPUT,
        interactions=Availability.INVALID_INPUT,
        sites=Availability.INVALID_INPUT,
        distance_map=Availability.INVALID_INPUT,
        displacement_vectors=Availability.INVALID_INPUT,
        evidence_cards=Availability.INVALID_INPUT,
    )
    return AnalysisReport(
        request.reference_selection,
        request.target_selection,
        input_quality,
        site_definitions=tuple(_canonical_site_definition(item) for item in request.site_definitions),
        availability=state,
        diagnostics=diagnostics,
        provenance=provenance,
    )


def _analysis_failure_report(
    request: AnalysisReportRequest,
    input_quality: InputQualityBundle,
    provenance: MethodProvenance,
    diagnostics: tuple[Diagnostic, ...],
) -> AnalysisReport:
    state = SectionAvailability(
        input_quality=Availability.AVAILABLE,
        analysis=Availability.NUMERICAL_FAILURE,
        msa=Availability.NUMERICAL_FAILURE,
        interactions=Availability.NUMERICAL_FAILURE,
        sites=Availability.NUMERICAL_FAILURE,
        distance_map=Availability.NUMERICAL_FAILURE,
        displacement_vectors=Availability.NUMERICAL_FAILURE,
        evidence_cards=Availability.NUMERICAL_FAILURE,
    )
    return AnalysisReport(
        request.reference_selection,
        request.target_selection,
        input_quality,
        site_definitions=tuple(_canonical_site_definition(item) for item in request.site_definitions),
        availability=state,
        diagnostics=diagnostics,
        provenance=provenance,
    )


def _section_failure(name: str, error: Exception) -> Diagnostic:
    return Diagnostic(
        code=f"report.{name}.failed",
        severity=DiagnosticSeverity.ERROR,
        message=f"The {name.replace('_', ' ')} section could not be calculated ({type(error).__name__}).",
        remediation="Review the input QC and the section-specific method settings.",
    )


def _merge_diagnostics(*groups: tuple[Diagnostic, ...]) -> tuple[Diagnostic, ...]:
    unique: dict[tuple[object, ...], Diagnostic] = {}
    for item in (diagnostic for group in groups for diagnostic in group):
        key = (
            item.code,
            item.severity,
            item.message,
            item.source_id,
            item.atom_id,
            item.residue_id,
            item.remediation,
        )
        unique.setdefault(key, item)
    return tuple(
        sorted(
            unique.values(),
            key=lambda item: (item.code, item.source_id or "", item.residue_id or "", item.atom_id or ""),
        )
    )


__all__ = ["ReportService"]
