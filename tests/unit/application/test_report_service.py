from __future__ import annotations

from pathlib import Path

import pytest

from structlens.application.analysis_service import AnalysisService
from structlens.application.dto import AnalysisReportRequest
from structlens.application.report_service import ReportService
from structlens.application.visualization_service import VisualizationService
from structlens.core.evidence import Availability
from structlens.core.interactions import InteractionThresholds
from structlens.core.models import AlignmentMode, AnalysisSettings, ResidueId
from structlens.core.msa import MSASettings
from structlens.core.parsing import InputSelection, StructureFormat, capture_snapshot, load_structure_evidence
from structlens.core.sites import SiteDefinition, SiteDefinitionMode
from structlens.plugin.visualization.renderer import VisualizationState

_FIXTURES = Path(__file__).parents[2] / "fixtures" / "parsing"


def _request(tmp_path: Path, *, invalid_reference: bool = False) -> AnalysisReportRequest:
    source = (_FIXTURES / "enriched.pdb").read_bytes()
    reference_path = tmp_path / "reference.pdb"
    target_path = tmp_path / "target.pdb"
    if invalid_reference:
        source = (
            b"ATOM      1  CA  ALA A   1       1.000   2.000   3.000  2.00 10.00           C  \nEND\n"
        )
    reference_path.write_bytes(source)
    target_path.write_bytes((_FIXTURES / "enriched.pdb").read_bytes())
    reference = capture_snapshot(reference_path)
    target = capture_snapshot(target_path)
    reference_selection = InputSelection(
        reference.content_id,
        reference.display_name,
        StructureFormat.PDB,
        "1" if not invalid_reference else "0",
        author_chain_ids=("A",),
        path=str(reference_path),
    )
    target_selection = InputSelection(
        target.content_id,
        target.display_name,
        StructureFormat.PDB,
        "1",
        author_chain_ids=("A",),
        path=str(target_path),
    )
    sites = ()
    if not invalid_reference:
        sites = (
            SiteDefinition(
                "active",
                "Active site",
                SiteDefinitionMode.KEY_RESIDUES,
                (ResidueId("reference", "1", "A", "100", None, "ALA"),),
            ),
        )
    return AnalysisReportRequest(
        reference,
        target,
        reference_selection,
        target_selection,
        analysis_settings=AnalysisSettings(),
        msa_settings=MSASettings(),
        interaction_thresholds=InteractionThresholds(),
        site_definitions=sites,
    )


def test_normal_report_service_populates_every_applicable_v03_section(tmp_path: Path) -> None:
    request = _request(tmp_path)

    report = ReportService().analyze(request)

    assert report.availability.input_quality is Availability.AVAILABLE
    assert report.availability.analysis is Availability.AVAILABLE
    assert report.availability.msa is Availability.AVAILABLE
    assert report.availability.interactions is Availability.AVAILABLE
    assert report.availability.sites is Availability.AVAILABLE
    assert report.availability.distance_map is Availability.AVAILABLE
    assert report.availability.displacement_vectors is Availability.AVAILABLE
    assert report.availability.evidence_cards is Availability.AVAILABLE
    assert report.analysis is not None and report.analysis.correspondences
    assert report.msa is not None and report.msa.columns
    assert report.distance_map is not None
    assert report.input_quality.reference is not None
    assert report.input_quality.target is not None
    assert report.sites and report.evidence_cards
    assert report.provenance.parameters["transform_direction"] == "target_to_reference"
    assert report.provenance.input_hashes == {
        "reference_raw": request.reference_snapshot.raw_sha256,
        "reference_content": request.reference_snapshot.content_id,
        "target_raw": request.target_snapshot.raw_sha256,
        "target_content": request.target_snapshot.content_id,
    }
    serialized = report.canonical_json_bytes()
    assert str(tmp_path).encode() not in serialized


def test_invalid_raw_input_keeps_typed_qc_and_never_runs_downstream_analysis(tmp_path: Path) -> None:
    class ForbiddenAnalysis:
        def analyze(self, *_: object, **__: object) -> object:
            raise AssertionError("downstream analysis must not run after invalid input QC")

    report = ReportService(analysis_service=ForbiddenAnalysis()).analyze(
        _request(tmp_path, invalid_reference=True)
    )

    assert report.availability.input_quality is Availability.INVALID_INPUT
    assert report.availability.analysis is Availability.INVALID_INPUT
    assert report.availability.msa is Availability.INVALID_INPUT
    assert report.analysis is None
    assert report.input_quality.reference is not None
    assert report.input_quality.reference.availability is Availability.INVALID_INPUT
    assert "occupancy.invalid" in {
        item.code for item in report.input_quality.reference.diagnostics
    }


def test_absent_optional_sites_are_not_applicable_not_an_empty_measurement(tmp_path: Path) -> None:
    request = _request(tmp_path)
    request = AnalysisReportRequest(
        request.reference_snapshot,
        request.target_snapshot,
        request.reference_selection,
        request.target_selection,
        analysis_settings=request.analysis_settings,
        msa_settings=request.msa_settings,
        interaction_thresholds=request.interaction_thresholds,
    )

    report = ReportService().analyze(request)

    assert report.availability.sites is Availability.NOT_APPLICABLE
    assert report.sites == ()


def test_request_rejects_selection_snapshot_hash_mismatch(tmp_path: Path) -> None:
    request = _request(tmp_path)
    wrong = InputSelection(
        "f" * 64,
        request.reference_selection.display_name,
        StructureFormat.PDB,
        "1",
    )

    with pytest.raises(ValueError, match="reference selection content"):
        AnalysisReportRequest(
            request.reference_snapshot,
            request.target_snapshot,
            wrong,
            request.target_selection,
        )


def test_optional_section_failure_is_contained_to_that_section(tmp_path: Path) -> None:
    class FailingMSA:
        def align(self, *_: object, **__: object) -> object:
            raise RuntimeError("alignment backend failed")

    report = ReportService(msa_engine=FailingMSA()).analyze(_request(tmp_path))

    # A presentation/report consumer must be able to use the independent
    # structural measurements even when the optional MSA backend fails.
    assert report.availability.analysis is Availability.AVAILABLE
    assert report.availability.msa is Availability.NUMERICAL_FAILURE
    assert report.availability.interactions is Availability.AVAILABLE
    assert report.availability.distance_map is Availability.AVAILABLE
    assert report.availability.displacement_vectors is Availability.AVAILABLE
    assert report.availability.evidence_cards is Availability.AVAILABLE
    assert report.analysis is not None
    assert report.interactions is not None
    assert report.distance_map is not None
    assert any(item.code == "report.msa.failed" for item in report.diagnostics)


def test_analysis_failure_is_typed_and_prevents_dependent_sections(tmp_path: Path) -> None:
    class FailingAnalysis:
        def analyze(self, *_: object, **__: object) -> object:
            raise RuntimeError("mapping backend failed")

    report = ReportService(analysis_service=FailingAnalysis()).analyze(_request(tmp_path))

    assert report.availability.input_quality is Availability.AVAILABLE
    assert report.availability.analysis is Availability.NUMERICAL_FAILURE
    assert report.availability.msa is Availability.NUMERICAL_FAILURE
    assert report.availability.interactions is Availability.NUMERICAL_FAILURE
    assert report.availability.sites is Availability.NUMERICAL_FAILURE
    assert report.availability.distance_map is Availability.NUMERICAL_FAILURE
    assert report.availability.displacement_vectors is Availability.NUMERICAL_FAILURE
    assert report.availability.evidence_cards is Availability.NUMERICAL_FAILURE
    assert report.analysis is None
    assert any(item.code == "report.analysis.failed" for item in report.diagnostics)


def test_identical_report_rerun_has_stable_identity_and_bytes(tmp_path: Path) -> None:
    request = _request(tmp_path)
    service = ReportService()

    first = service.analyze(request)
    second = service.analyze(request)

    assert first.report_id == second.report_id
    assert first.canonical_json_bytes() == second.canonical_json_bytes()


def test_pairwise_report_rejects_a_selection_with_multiple_protein_chains(tmp_path: Path) -> None:
    source = _two_chain_pdb()
    reference_path = tmp_path / "multi-reference.pdb"
    target_path = tmp_path / "multi-target.pdb"
    reference_path.write_bytes(source)
    target_path.write_bytes(source)
    reference = capture_snapshot(reference_path)
    target = capture_snapshot(target_path)
    reference_selection = InputSelection(
        reference.content_id,
        reference.display_name,
        StructureFormat.PDB,
        "1",
        author_chain_ids=("A", "B"),
        path=str(reference_path),
    )
    target_selection = InputSelection(
        target.content_id,
        target.display_name,
        StructureFormat.PDB,
        "1",
        author_chain_ids=("A", "B"),
        path=str(target_path),
    )

    class ForbiddenAnalysis:
        def analyze(self, *_: object, **__: object) -> object:
            raise AssertionError("multi-chain selections must stop before analysis")

    report = ReportService(analysis_service=ForbiddenAnalysis()).analyze(
        AnalysisReportRequest(reference, target, reference_selection, target_selection)
    )

    assert report.availability.input_quality is Availability.INVALID_INPUT
    assert report.availability.analysis is Availability.INVALID_INPUT
    assert report.analysis is None
    assert {item.code for item in report.diagnostics} >= {
        "report.input.reference.chain_count",
        "report.input.target.chain_count",
    }


def test_empty_interaction_and_vector_results_still_mean_measurement_available(tmp_path: Path) -> None:
    report = ReportService().analyze(_request(tmp_path))

    # The absence of a detected feature is not the same as an unavailable or
    # failed measurement. This distinction prevents downstream code from
    # treating scientifically valid negatives as missing values.
    assert report.availability.interactions is Availability.AVAILABLE
    assert report.interactions is not None
    assert not report.interactions.reference_interactions
    assert not report.interactions.target_interactions
    assert not report.interactions.differences
    assert report.availability.displacement_vectors is Availability.AVAILABLE
    assert report.displacement_vectors == ()


def test_visualization_selection_consumes_the_canonical_report(tmp_path: Path) -> None:
    report = ReportService().analyze(_request(tmp_path))

    selected = VisualizationService().select(report, VisualizationState())

    assert selected == report.analysis.correspondences


def test_visualization_selection_preserves_legacy_analysis_result_flow(tmp_path: Path) -> None:
    request = _request(tmp_path)
    reference = load_structure_evidence(request.reference_snapshot, selection=request.reference_selection)
    target = load_structure_evidence(request.target_snapshot, selection=request.target_selection)
    result = AnalysisService().analyze(
        reference.protein_structure,
        target.protein_structure,
        request.analysis_settings,
        reference_chain_id=reference.protein_structure.chains[0].chain_id,
        target_chain_id=target.protein_structure.chains[0].chain_id,
    )

    selected = VisualizationService().select(result, VisualizationState())

    assert selected == result.correspondences


def test_visualization_selection_handles_reports_without_analysis(tmp_path: Path) -> None:
    report = ReportService().analyze(_request(tmp_path, invalid_reference=True))

    selected = VisualizationService().select(report, VisualizationState())

    assert report.analysis is None
    assert selected == ()


def test_report_identity_does_not_change_when_identical_sources_are_renamed(tmp_path: Path) -> None:
    source = (_FIXTURES / "enriched.pdb").read_bytes()

    def renamed_request(prefix: str) -> AnalysisReportRequest:
        reference_path = tmp_path / f"{prefix}_reference.pdb"
        target_path = tmp_path / f"{prefix}_target.pdb"
        reference_path.write_bytes(source)
        target_path.write_bytes(source)
        reference = capture_snapshot(reference_path)
        target = capture_snapshot(target_path)
        return AnalysisReportRequest(
            reference,
            target,
            InputSelection(reference.content_id, reference.display_name, StructureFormat.PDB, "1", ("A",)),
            InputSelection(target.content_id, target.display_name, StructureFormat.PDB, "1", ("A",)),
        )

    first = ReportService().analyze(renamed_request("first"))
    second = ReportService().analyze(renamed_request("second"))

    assert first.report_id == second.report_id
    assert first.canonical_json_bytes() == second.canonical_json_bytes()


def test_unresolved_site_is_not_a_successful_zero_measurement(tmp_path: Path) -> None:
    request = _request(tmp_path)
    unknown = SiteDefinition(
        "missing-ligand",
        "Missing ligand site",
        SiteDefinitionMode.LIGAND_RADIUS,
        ligand_id="DOES_NOT_EXIST",
        radius_angstrom=5.0,
    )
    request = AnalysisReportRequest(
        request.reference_snapshot,
        request.target_snapshot,
        request.reference_selection,
        request.target_selection,
        site_definitions=(unknown,),
    )

    report = ReportService().analyze(request)

    assert report.availability.sites is Availability.NOT_DETECTED
    assert report.sites == ()
    assert any(item.code == "report.sites.missing-ligand.not_detected" for item in report.diagnostics)


def test_evidence_cards_preserve_failed_interaction_availability(tmp_path: Path) -> None:
    class FailingInteractions:
        def detect(self, *_: object, **__: object) -> object:
            raise RuntimeError("interaction backend failed")

        def compare(self, *_: object, **__: object) -> object:
            raise AssertionError("comparison cannot run after detection failure")

    report = ReportService(interaction_service=FailingInteractions()).analyze(_request(tmp_path))

    assert report.availability.interactions is Availability.NUMERICAL_FAILURE
    assert report.evidence_cards
    assert all("interactions" in card.quality.unavailable_sections for card in report.evidence_cards)
    assert all("interactions" not in card.quality.available_sections for card in report.evidence_cards)


def test_site_definition_is_preserved_in_report_and_provenance(tmp_path: Path) -> None:
    report = ReportService().analyze(_request(tmp_path))

    assert report.site_definitions[0].site_id == "active"
    serialized = report.provenance.parameters["site_definitions"]
    assert serialized[0]["site_id"] == "active"
    assert serialized[0]["mode"] == "key_residues"


def test_manual_mapping_runs_through_the_same_canonical_report_path(tmp_path: Path) -> None:
    request = _request(tmp_path)
    manual = AnalysisReportRequest(
        request.reference_snapshot,
        request.target_snapshot,
        request.reference_selection,
        request.target_selection,
        analysis_settings=AnalysisSettings(alignment_mode=AlignmentMode.MANUAL),
        manual_pairs=((
            ResidueId("reference", "1", "A", "100", None, "ALA"),
            ResidueId("target", "1", "A", "100", None, "ALA"),
        ),),
    )

    report = ReportService().analyze(manual)

    assert report.availability.analysis is Availability.AVAILABLE
    assert report.analysis is not None
    assert len(report.analysis.correspondences) == 1
    assert report.provenance.parameters["manual_pairs"]


def _two_chain_pdb() -> bytes:
    """Return a tiny valid two-chain PDB for exact-chain selection tests."""

    def atom(serial: int, name: str, residue: str, chain: str, sequence: int, x: float, element: str) -> str:
        return (
            f"ATOM  {serial:5d} {name:^4s} {residue:>3s} {chain}{sequence:4d}    "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{10.0:6.2f}          {element:>2s}\n"
        )

    lines = ["HEADER    TWO CHAIN TEST\n", "MODEL        1\n"]
    serial = 1
    for chain, offset in (("A", 0.0), ("B", 20.0)):
        for sequence, residue in ((1, "ALA"), (2, "GLY")):
            for name, delta, element in (("N", 0.0, "N"), ("CA", 1.3, "C"), ("C", 2.6, "C"), ("O", 3.6, "O")):
                lines.append(atom(serial, name, residue, chain, sequence, offset + delta + (sequence - 1) * 3.8, element))
                serial += 1
        lines.append(f"TER   {serial:5d}      GLY {chain}{2:4d}    \n")
        serial += 1
    lines.extend(("ENDMDL\n", "END\n"))
    return "".join(lines).encode("ascii")
