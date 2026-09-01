"""Small, hand-checked fixtures shared by Task 13 round-four regressions."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from structlens.application.dto import AnalysisReportRequest
from structlens.application.report_service import ReportService
from structlens.core.difference_maps import ResidueDisplacementVector
from structlens.core.evidence import (
    Availability,
    Diagnostic,
    DiagnosticSeverity,
    InteractionEvidence,
)
from structlens.core.interactions import (
    InteractionChange,
    InteractionDifference,
    InteractionRecord,
    InteractionType,
    ReferenceInteractionKey,
)
from structlens.core.models import ResidueId
from structlens.core.parsing import InputSelection, StructureFormat
from structlens.core.reports import AnalysisReport
from structlens.core.sites import SiteMetrics
from tests.unit.application.test_report_service import _request as _service_request

FIXTURES = Path(__file__).parents[2] / "fixtures" / "parsing"


def request(tmp_path: Path) -> AnalysisReportRequest:
    """Create a canonical request whose source paths are disposable."""

    return _service_request(tmp_path)


def rich_report(tmp_path: Path) -> tuple[AnalysisReportRequest, AnalysisReport]:
    """Return a report with one value in every section used by Qt assertions."""

    current_request = request(tmp_path)
    report = ReportService().analyze(current_request)
    assert report.analysis is not None
    correspondence = report.analysis.correspondences[0]
    assert correspondence.reference is not None and correspondence.target is not None

    reference_residue = correspondence.reference
    target_residue = correspondence.target
    record = InteractionRecord(
        "target",
        InteractionType.HBOND_GEOMETRIC,
        target_residue,
        None,
        "N",
        "O",
        2.70,
        evidence_mode="heavy_atom_geometry",
    )
    interaction = InteractionDifference(
        ReferenceInteractionKey(InteractionType.HBOND_GEOMETRIC, _residue_label(reference_residue)),
        InteractionChange.GAINED,
        target_record=record,
    )
    vector = ResidueDisplacementVector(
        reference_position=_residue_label(reference_residue),
        reference_residue=reference_residue,
        target_residue=target_residue,
        start_xyz=(0.0, 0.0, 0.0),
        end_xyz=(3.0, 4.0, 0.0),
        vector_xyz=(3.0, 4.0, 0.0),
        magnitude_angstrom=5.0,
    )
    existing_site = report.sites[0] if report.sites else SiteMetrics("round4", "target", 0, 0.0)
    site = replace(existing_site, site_id="round4-site", atomic_envelope_volume_angstrom3=12.5)
    diagnostic = Diagnostic(
        "round4.integration.diagnostic",
        DiagnosticSeverity.WARNING,
        "A controlled round-four diagnostic remains visible.",
        remediation="Inspect the fixture evidence.",
    )
    availability = replace(
        report.availability,
        sites=Availability.AVAILABLE,
        interactions=Availability.AVAILABLE,
        displacement_vectors=Availability.AVAILABLE,
    )
    rich = replace(
        report,
        interactions=InteractionEvidence(differences=(interaction,), target_interactions=(record,)),
        sites=(site,),
        displacement_vectors=(vector,),
        diagnostics=(diagnostic,),
        availability=availability,
    )
    return current_request, rich


def gapped_card(report: AnalysisReport) -> AnalysisReport:
    """Set a correspondence alignment index apart from its sequence index."""

    assert report.evidence_cards
    card = report.evidence_cards[0]
    shifted = replace(card, sequence=replace(card.sequence, alignment_index=17))
    availability = replace(report.availability, evidence_cards=Availability.AVAILABLE)
    return replace(report, evidence_cards=(shifted,), availability=availability)


def multi_model_chain_pdb() -> bytes:
    """Two models and two chains, so a reset to the first item is observable."""

    def atom(serial: int, name: str, residue: str, chain: str, seq: int, x: float) -> str:
        return (
            f"ATOM  {serial:5d} {name:^4s} {residue:>3s} {chain}{seq:4d}    "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{10.0:6.2f}          "
            f"{name[0]:>2s}\n"
        )

    lines = ["HEADER    ROUND FOUR MODEL CHAIN FIXTURE\n"]
    serial = 1
    for model in (1, 2):
        lines.append(f"MODEL        {model}\n")
        for chain, offset in (("A", 0.0), ("B", 20.0)):
            for seq, residue in ((1, "ALA"), (2, "GLY")):
                for name, delta in (("N", 0.0), ("CA", 1.3), ("C", 2.6), ("O", 3.6)):
                    lines.append(atom(serial, name, residue, chain, seq, offset + delta + seq * 3.8))
                    serial += 1
            lines.append(f"TER   {serial:5d}      GLY {chain}{2:4d}\n")
            serial += 1
        lines.append("ENDMDL\n")
    lines.append("END\n")
    return "".join(lines).encode("ascii")


def non_first_request_paths(tmp_path: Path) -> tuple[Path, Path]:
    source = multi_model_chain_pdb()
    reference = tmp_path / "multi-reference.pdb"
    target = tmp_path / "multi-target.pdb"
    reference.write_bytes(source)
    target.write_bytes(source)
    return reference, target


def selection(snapshot: object, *, model: str, chain: str, path: Path) -> InputSelection:
    """Build a hand-checked selection for direct snapshot-boundary tests."""

    assert hasattr(snapshot, "content_id") and hasattr(snapshot, "display_name")
    assert hasattr(snapshot, "logical_format")
    return InputSelection(
        snapshot.content_id,  # type: ignore[attr-defined]
        snapshot.display_name,  # type: ignore[attr-defined]
        StructureFormat(snapshot.logical_format),  # type: ignore[attr-defined]
        model,
        author_chain_ids=(chain,),
        path=str(path),
    )


def _residue_label(residue: ResidueId) -> str:
    return f"{residue.chain_id}:{residue.auth_seq_id}{residue.insertion_code or ''}"


__all__ = [
    "FIXTURES",
    "gapped_card",
    "multi_model_chain_pdb",
    "non_first_request_paths",
    "rich_report",
    "request",
    "selection",
]
