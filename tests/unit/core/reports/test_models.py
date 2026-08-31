from __future__ import annotations

import json
import math
from dataclasses import FrozenInstanceError

import numpy as np
import pytest

from structlens.core.difference_maps import DistanceDifferenceMatrix, ResidueDisplacementVector
from structlens.core.evidence import (
    Availability,
    Diagnostic,
    DiagnosticSeverity,
    EvidenceCard,
    EvidenceQuality,
    InteractionEvidence,
    SequenceEvidence,
    SiteEvidence,
    StructureEvidence,
)
from structlens.core.interactions import (
    InteractionChange,
    InteractionDifference,
    InteractionRecord,
    InteractionType,
    ReferenceInteractionKey,
)
from structlens.core.models import (
    AnalysisResult,
    CorrespondenceStatus,
    MutationEvent,
    MutationKind,
    ResidueCorrespondence,
    ResidueId,
    StructuralTransform,
)
from structlens.core.msa import (
    AnalysisSequence,
    MSAColumn,
    MSAResidueCell,
    MultipleSequenceAlignment,
    SequenceResidueRef,
)
from structlens.core.parsing import InputSelection, StructureFormat
from structlens.core.provenance import MethodProvenance
from structlens.core.quality import StructureQualityReport
from structlens.core.reports import (
    AnalysisReport,
    AnalysisSnapshot,
    CorrespondenceSnapshot,
    DistanceMapSnapshot,
    InputQualityBundle,
    SectionAvailability,
)
from structlens.core.sites import SiteMetrics


def _selection(content: str, name: str, *, path: str) -> InputSelection:
    return InputSelection(content, name, StructureFormat.PDB, "1", path=path)


def _reference_residue() -> ResidueId:
    return ResidueId("reference", "1", "A", "10", None, "ALA")


def _target_residue() -> ResidueId:
    return ResidueId("target", "1", "B", "20", None, "SER")


def _correspondence() -> ResidueCorrespondence:
    return ResidueCorrespondence(
        0,
        _reference_residue(),
        _target_residue(),
        "A",
        "S",
        CorrespondenceStatus.SUBSTITUTION,
        sequence_score=2.0,
        ca_displacement_angstrom=1.25,
        mapping_source="sequence",
    )


def _analysis(correspondence: ResidueCorrespondence, provenance: dict[str, str]) -> AnalysisResult:
    return AnalysisResult(
        "reference",
        "target",
        (correspondence,),
        (),
        0.5,
        1.0,
        "sequence-guided",
        sequence_similarity=0.75,
        strict_rmsd_angstrom=1.25,
        mapped_residue_count=1,
        provenance=provenance,
        transform=StructuralTransform(translation=(1.0, 2.0, 3.0)),
    )


def _distance_map(source: np.ndarray) -> DistanceDifferenceMatrix:
    zero = np.zeros_like(source)
    return DistanceDifferenceMatrix(("A:10", "A:11"), zero, source, source, np.ones_like(source, dtype=bool))


def _report(
    analysis: AnalysisResult,
    matrix: DistanceDifferenceMatrix,
) -> AnalysisReport:
    reference = _selection("a" * 64, "reference.pdb", path="C:/private/reference.pdb")
    target = _selection("b" * 64, "target.pdb", path="C:/private/target.pdb")
    provenance = MethodProvenance(
        "structlens.analysis_report",
        "1",
        parameters={"transform_direction": "target_to_reference"},
        input_hashes={
            "reference_raw": "c" * 64,
            "reference_content": reference.content_id,
            "target_raw": "d" * 64,
            "target_content": target.content_id,
        },
        analyzed_representation="asymmetric_unit_pair",
    )
    return AnalysisReport(
        reference,
        target,
        InputQualityBundle(
            StructureQualityReport(Availability.AVAILABLE),
            StructureQualityReport(Availability.AVAILABLE),
        ),
        AnalysisSnapshot.from_result(analysis),
        distance_map=DistanceMapSnapshot.from_matrix(matrix),
        availability=SectionAvailability(
            input_quality=Availability.AVAILABLE,
            analysis=Availability.AVAILABLE,
            distance_map=Availability.AVAILABLE,
        ),
        diagnostics=(
            Diagnostic("report.test", DiagnosticSeverity.INFO, "Test diagnostic."),
        ),
        provenance=provenance,
    )


def test_correspondence_snapshot_is_a_deep_immutable_value() -> None:
    source = _correspondence()
    snapshot = CorrespondenceSnapshot.from_record(source)

    source.ca_displacement_angstrom = 99.0
    source.is_outlier = True

    assert snapshot.ca_displacement_angstrom == 1.25
    assert snapshot.is_outlier is False
    with pytest.raises(FrozenInstanceError):
        snapshot.is_outlier = True  # type: ignore[misc]


def test_report_copies_mutable_inputs_and_keeps_deterministic_identity() -> None:
    correspondence = _correspondence()
    legacy_provenance = {"engine": "test"}
    matrix_values = np.asarray(((0.0, 1.0), (1.0, 0.0)), dtype=float)
    matrix = _distance_map(matrix_values)
    report = _report(_analysis(correspondence, legacy_provenance), matrix)
    before = report.canonical_json_bytes()
    before_id = report.report_id

    correspondence.reference_one_letter = "G"
    legacy_provenance["engine"] = "mutated"
    matrix_values[0, 1] = 500.0

    assert report.analysis is not None
    assert report.analysis.correspondences[0].reference_one_letter == "A"
    assert report.analysis.legacy_provenance["engine"] == "test"
    assert report.distance_map is not None
    assert report.distance_map.delta_angstrom[0][1] == 1.0
    assert report.canonical_json_bytes() == before
    assert report.report_id == before_id


def test_report_json_is_canonical_typed_and_does_not_serialize_input_paths() -> None:
    report = _report(_analysis(_correspondence(), {}), _distance_map(np.asarray(((0.0, 1.0), (1.0, 0.0)))))

    payload = report.to_json()
    encoded = report.canonical_json_bytes()

    assert json.loads(encoded) == payload
    assert payload["report_id"] == report.report_id
    assert payload["availability"]["sites"] == "not_applicable"
    assert payload["reference_selection"]["selection_id"] == report.reference_selection.selection_id
    assert "path" not in payload["reference_selection"]
    assert b"C:/private" not in encoded
    assert encoded == report.canonical_json_bytes()


def test_report_rejects_incoherent_selection_hash_provenance() -> None:
    report = _report(_analysis(_correspondence(), {}), _distance_map(np.asarray(((0.0, 1.0), (1.0, 0.0)))))
    wrong = MethodProvenance(
        "structlens.analysis_report",
        "1",
        input_hashes={
            "reference_content": "0" * 64,
            "target_content": report.target_selection.content_id,
        },
    )

    with pytest.raises(ValueError, match="reference selection content"):
        AnalysisReport(
            report.reference_selection,
            report.target_selection,
            report.input_quality,
            report.analysis,
            availability=SectionAvailability(analysis=Availability.AVAILABLE),
            provenance=wrong,
        )


def test_report_rejects_available_section_without_its_typed_payload() -> None:
    report = _report(_analysis(_correspondence(), {}), _distance_map(np.asarray(((0.0, 1.0), (1.0, 0.0)))))

    with pytest.raises(ValueError, match="analysis availability"):
        AnalysisReport(
            report.reference_selection,
            report.target_selection,
            report.input_quality,
            availability=SectionAvailability(analysis=Availability.AVAILABLE),
        )


def test_report_rejects_payload_hidden_behind_unavailable_state() -> None:
    report = _report(_analysis(_correspondence(), {}), _distance_map(np.asarray(((0.0, 1.0), (1.0, 0.0)))))

    with pytest.raises(ValueError, match="analysis availability"):
        AnalysisReport(
            report.reference_selection,
            report.target_selection,
            report.input_quality,
            report.analysis,
        )


def test_report_serializes_rich_typed_sections_and_compatibility_aliases() -> None:
    reference_residue = _reference_residue()
    target_residue = _target_residue()
    correspondence = _correspondence()
    mutation = MutationEvent(
        0,
        MutationKind.SUBSTITUTION,
        reference_residue,
        target_residue,
        "A",
        "S",
        "ALA 10",
        "SER 20",
        "A10S",
        1,
        99,
        "nonconservative",
    )
    analysis = AnalysisResult(
        "reference",
        "target",
        (correspondence,),
        (mutation,),
        0.5,
        1.0,
        "sequence-guided",
        sequence_similarity=0.75,
        strict_rmsd_angstrom=1.25,
        refined_rmsd_angstrom=1.0,
        mapped_residue_count=1,
        refined_residue_count=1,
        excluded_alignment_indices=(4,),
        tm_score=0.9,
        provenance={"engine": "test"},
        transform=StructuralTransform(translation=(1.0, 2.0, 3.0)),
        method_provenance=MethodProvenance(
            "alignment.method",
            "2",
            input_hashes={"alignment_input": "e" * 64},
        ),
    )
    residue_ref = SequenceResidueRef(0, "A", reference_residue)
    target_ref = SequenceResidueRef(0, "S", target_residue)
    msa = MultipleSequenceAlignment(
        sequences=(
            AnalysisSequence("reference", "A", "A", (residue_ref,), "structure"),
            AnalysisSequence("target", "B", "S", (target_ref,), "structure"),
        ),
        aligned_rows=(("reference", "A"), ("target", "S")),
        columns=(
            MSAColumn(
                index=0,
                reference_label="A:10",
                reference_residue=reference_residue,
                cells=(
                    MSAResidueCell("reference", 0, residue_ref, "A"),
                    MSAResidueCell("target", 0, target_ref, "S"),
                ),
                non_gap_count=2,
                gap_fraction=0.0,
                ambiguous_fraction=0.0,
                conservation_score=0.5,
                entropy_bits=0.7,
            ),
        ),
        reference_structure_id="reference",
        algorithm="muscle5",
        provenance=("manual",),
    )
    reference_record = InteractionRecord(
        "reference",
        InteractionType.HYDROPHOBIC,
        reference_residue,
        None,
        "CA",
        None,
        4.0,
    )
    difference = InteractionDifference(
        ReferenceInteractionKey(InteractionType.HYDROPHOBIC, "A:10"),
        InteractionChange.LOST,
        reference_record=reference_record,
    )
    site = SiteMetrics(
        "active",
        "target",
        1,
        1.0,
        global_frame_backbone_rmsd_angstrom=0.4,
        site_fitted_backbone_rmsd_angstrom=0.2,
        centroid_displacement_angstrom=0.5,
        radius_of_gyration_angstrom=1.2,
        atomic_envelope_volume_angstrom3=25.0,
        sasa_angstrom2=33.0,
        polar_residue_fraction=1.0,
        charged_residue_fraction=0.0,
    )
    vector = ResidueDisplacementVector(
        "A:10",
        reference_residue,
        target_residue,
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 1.0),
        (0.0, 0.0, 1.0),
        1.0,
    )
    card = EvidenceCard(
        residue_ref,
        target_id="target",
        sequence=SequenceEvidence(
            reference_one_letter="A",
            target_one_letter="S",
            alignment_index=0,
            sequence_identity=0.5,
            source_refs=(residue_ref, target_ref),
        ),
        structure=StructureEvidence(ca_displacement_angstrom=1.0, available=True),
        interactions=InteractionEvidence((difference,), (reference_record,), ()),
        site=SiteEvidence((site,)),
        quality=EvidenceQuality(
            overall_status="available",
            available_sections=("sequence", "structure", "sites"),
            source_count=3,
        ),
        provenance=("structlens.analysis_report", "sequence-guided"),
    )
    report = AnalysisReport(
        _selection("a" * 64, "reference.pdb", path="C:/private/reference.pdb"),
        _selection("b" * 64, "target.pdb", path="C:/private/target.pdb"),
        InputQualityBundle(
            StructureQualityReport(Availability.AVAILABLE),
            StructureQualityReport(Availability.AVAILABLE),
        ),
        analysis=AnalysisSnapshot.from_result(analysis),
        msa=msa,
        interactions=InteractionEvidence((difference,), (reference_record,), ()),
        sites=SiteEvidence((site,)),
        distance_map=DistanceMapSnapshot.from_matrix(_distance_map(np.asarray(((0.0, 1.0), (1.0, 0.0))))),
        displacement_vectors=(vector,),
        evidence_cards=(card,),
        availability=SectionAvailability(
            input_quality="available",
            analysis="available",
            msa="available",
            interactions="available",
            sites="available",
            distance_map="available",
            displacement_vectors="available",
            evidence_cards="available",
        ),
        provenance=MethodProvenance(
            "structlens.analysis_report",
            "1",
            input_hashes={
                "reference_selection_content": "a" * 64,
                "target-selection-content": "b" * 64,
            },
        ),
    )

    payload = report.to_json()

    assert payload["analysis"]["mutations"][0]["canonical_notation"] == "A10S"
    assert payload["analysis"]["method_provenance"]["method_id"] == "alignment.method"
    assert payload["msa"]["columns"][0]["cells"][1]["character"] == "S"
    assert payload["interactions"]["differences"][0]["target_record"] is None
    assert payload["sites"]["metrics"][0]["atomic_envelope_volume_angstrom3"] == 25.0
    assert payload["displacement_vectors"][0]["magnitude_angstrom"] == 1.0
    assert payload["evidence_cards"][0]["quality"]["available_sections"] == ["sequence", "structure", "sites"]
    assert report.availability.input_quality is Availability.AVAILABLE
    assert report.availability.vectors is Availability.AVAILABLE
    assert report.vectors[0].magnitude_angstrom == 1.0


def test_report_accepts_direct_tuple_sections_and_absent_provenance() -> None:
    reference = _selection("a" * 64, "reference.pdb", path="C:/private/reference.pdb")
    target = _selection("b" * 64, "target.pdb", path="C:/private/target.pdb")
    residue = _reference_residue()
    interaction = InteractionRecord(
        "reference",
        InteractionType.HYDROPHOBIC,
        residue,
        None,
        "CA",
        None,
        4.0,
    )
    site = SiteMetrics("active", "reference", 1, 1.0)
    report = AnalysisReport(
        reference,
        target,
        InputQualityBundle(
            StructureQualityReport(Availability.AVAILABLE),
            StructureQualityReport(Availability.AVAILABLE),
        ),
        interactions=(interaction,),
        sites=(site,),
        availability=SectionAvailability(
            interactions=Availability.AVAILABLE,
            sites=Availability.AVAILABLE,
        ),
    )

    payload = report.to_json()

    assert payload["interactions"][0]["interaction_type"] == "hydrophobic"
    assert payload["sites"][0]["site_id"] == "active"
    assert payload["provenance"] is None


@pytest.mark.parametrize(
    ("factory", "message"),
    [
        (
            lambda: CorrespondenceSnapshot(
                alignment_index=-1,
                reference=None,
                target=None,
                reference_one_letter=None,
                target_one_letter=None,
                status="substitution",
            ),
            "alignment_index must be non-negative",
        ),
        (lambda: CorrespondenceSnapshot.from_record(object()), "record must be a ResidueCorrespondence"),
        (
            lambda: AnalysisSnapshot(
                "reference",
                "target",
                ("bad",),
                (),
                0.5,
                1.0,
                "sequence-guided",
            ),
            "correspondences must contain CorrespondenceSnapshot values",
        ),
        (
            lambda: AnalysisSnapshot(
                "reference",
                "target",
                (),
                ("bad",),
                0.5,
                1.0,
                "sequence-guided",
            ),
            "mutations must contain MutationEvent values",
        ),
        (
            lambda: AnalysisSnapshot(
                "reference",
                "target",
                (),
                (),
                0.5,
                1.0,
                "sequence-guided",
                mapped_residue_count=-1,
            ),
            "mapped_residue_count must be non-negative or None",
        ),
        (
            lambda: AnalysisSnapshot(
                "reference",
                "target",
                (),
                (),
                0.5,
                1.0,
                "sequence-guided",
                method_provenance="bad",
            ),
            "method_provenance must be MethodProvenance or None",
        ),
        (lambda: DistanceMapSnapshot.from_matrix(object()), "matrix must be a DistanceDifferenceMatrix"),
        (
            lambda: DistanceMapSnapshot(("A:10",), ((0.0, 1.0),), ((0.0, 1.0),), ((0.0, 1.0),), ((True, False),)),
            "reference_distances_angstrom shape must match reference_positions",
        ),
        (
            lambda: DistanceMapSnapshot(("A:10",), ((math.nan,),), ((0.0,),), ((0.0,),), ((True,),)),
            "reference_distances_angstrom must contain finite values",
        ),
        (
            lambda: DistanceMapSnapshot(("A:10",), ((0.0,),), ((0.0,),), ((0.0,),), ((True, False),)),
            "valid_mask shape must match reference_positions",
        ),
        (
            lambda: InputQualityBundle("bad", StructureQualityReport(Availability.AVAILABLE)),
            "reference and target must be StructureQualityReport values",
        ),
        (lambda: SectionAvailability(input_quality="bad"), "unknown input_quality"),
        (
            lambda: AnalysisReport(
                _selection("a" * 64, "reference.pdb", path="C:/private/reference.pdb"),
                _selection("b" * 64, "target.pdb", path="C:/private/target.pdb"),
                InputQualityBundle(
                    StructureQualityReport(Availability.AVAILABLE),
                    StructureQualityReport(Availability.AVAILABLE),
                ),
                interactions=("bad",),
            ),
            "interactions must contain typed interaction values",
        ),
        (
            lambda: AnalysisReport(
                _selection("a" * 64, "reference.pdb", path="C:/private/reference.pdb"),
                _selection("b" * 64, "target.pdb", path="C:/private/target.pdb"),
                InputQualityBundle(
                    StructureQualityReport(Availability.AVAILABLE),
                    StructureQualityReport(Availability.AVAILABLE),
                ),
                sites=("bad",),
            ),
            "sites must contain SiteMetrics values",
        ),
        (
            lambda: AnalysisReport(
                _selection("a" * 64, "reference.pdb", path="C:/private/reference.pdb"),
                _selection("b" * 64, "target.pdb", path="C:/private/target.pdb"),
                InputQualityBundle(
                    StructureQualityReport(Availability.AVAILABLE),
                    StructureQualityReport(Availability.AVAILABLE),
                ),
                distance_map="bad",
            ),
            "distance_map must be a DistanceMapSnapshot or None",
        ),
        (
            lambda: AnalysisReport(
                _selection("a" * 64, "reference.pdb", path="C:/private/reference.pdb"),
                _selection("b" * 64, "target.pdb", path="C:/private/target.pdb"),
                InputQualityBundle(
                    StructureQualityReport(Availability.AVAILABLE),
                    StructureQualityReport(Availability.AVAILABLE),
                ),
                displacement_vectors=("bad",),
            ),
            "vectors must contain ResidueDisplacementVector values",
        ),
        (
            lambda: AnalysisReport(
                _selection("a" * 64, "reference.pdb", path="C:/private/reference.pdb"),
                _selection("b" * 64, "target.pdb", path="C:/private/target.pdb"),
                InputQualityBundle(
                    StructureQualityReport(Availability.AVAILABLE),
                    StructureQualityReport(Availability.AVAILABLE),
                ),
                evidence_cards=("bad",),
            ),
            "evidence_cards must contain EvidenceCard values",
        ),
        (
            lambda: AnalysisReport(
                _selection("a" * 64, "reference.pdb", path="C:/private/reference.pdb"),
                _selection("b" * 64, "target.pdb", path="C:/private/target.pdb"),
                InputQualityBundle(
                    StructureQualityReport(Availability.AVAILABLE),
                    StructureQualityReport(Availability.AVAILABLE),
                ),
                availability="bad",
            ),
            "availability must be SectionAvailability",
        ),
        (
            lambda: AnalysisReport(
                _selection("a" * 64, "reference.pdb", path="C:/private/reference.pdb"),
                _selection("b" * 64, "target.pdb", path="C:/private/target.pdb"),
                InputQualityBundle(
                    StructureQualityReport(Availability.AVAILABLE),
                    StructureQualityReport(Availability.AVAILABLE),
                ),
                diagnostics=("bad",),
            ),
            "diagnostics must contain Diagnostic values",
        ),
        (
            lambda: AnalysisReport(
                _selection("a" * 64, "reference.pdb", path="C:/private/reference.pdb"),
                _selection("b" * 64, "target.pdb", path="C:/private/target.pdb"),
                InputQualityBundle(
                    StructureQualityReport(Availability.AVAILABLE),
                    StructureQualityReport(Availability.AVAILABLE),
                ),
                provenance="bad",
            ),
            "provenance must be MethodProvenance or None",
        ),
    ],
)
def test_report_models_reject_invalid_values(factory: object, message: str) -> None:
    with pytest.raises((TypeError, ValueError), match=message):
        factory()
