from __future__ import annotations

from pathlib import Path

from structlens.application.analysis_service import AnalysisService
from structlens.application.dto import AnalysisReportRequest
from structlens.application.report_service import ReportService
from structlens.core.evidence import Availability
from structlens.core.models import (
    AlignmentMode,
    AnalysisSettings,
    CorrespondenceStatus,
    ResidueCorrespondence,
)
from structlens.core.parsing import InputSelection, StructureFormat, capture_snapshot
from structlens.integrations.usalign.adapter import USAlignAlignmentResult


def test_structural_report_passes_only_the_selected_model_and_chain(tmp_path: Path) -> None:
    source_path = tmp_path / "multi-model.pdb"
    source_path.write_bytes(_multi_model_two_chain_pdb())
    snapshot = capture_snapshot(source_path)
    selection = InputSelection(
        snapshot.content_id,
        snapshot.display_name,
        StructureFormat.PDB,
        "2",
        author_chain_ids=("B",),
        path=str(source_path),
    )

    class CapturingAdapter:
        def align(self, reference: object, target: object, settings: object) -> USAlignAlignmentResult:
            for chain in (reference, target):
                assert chain.model_id == "2"  # type: ignore[attr-defined]
                assert chain.chain_id == "B"  # type: ignore[attr-defined]
                assert tuple(
                    atom.coordinate[0]
                    for residue in chain.residue_records  # type: ignore[attr-defined]
                    for atom in residue.atoms
                    if atom.name == "CA"
                ) == (120.0, 123.80000305175781)
            correspondences = tuple(
                ResidueCorrespondence(
                    index,
                    reference.residue_records[index].residue_id,  # type: ignore[attr-defined]
                    target.residue_records[index].residue_id,  # type: ignore[attr-defined]
                    reference.residue_records[index].one_letter,  # type: ignore[attr-defined]
                    target.residue_records[index].one_letter,  # type: ignore[attr-defined]
                    CorrespondenceStatus.CONSERVED,
                    mapping_source="US-align",
                )
                for index in range(len(reference.residue_records))  # type: ignore[attr-defined]
            )
            return USAlignAlignmentResult(correspondences, 1.0, None, "test", {})

    source_path.write_bytes(b"REPLACED AFTER SNAPSHOT\n")
    report = ReportService(analysis_service=AnalysisService(CapturingAdapter())).analyze(
        AnalysisReportRequest(
            snapshot,
            snapshot,
            selection,
            selection,
            analysis_settings=AnalysisSettings(alignment_mode=AlignmentMode.STRUCTURE),
        )
    )

    assert report.availability.analysis is Availability.AVAILABLE
    assert report.analysis is not None
    assert str(source_path).encode() not in report.canonical_json_bytes()


def _multi_model_two_chain_pdb() -> bytes:
    lines = ["HEADER    MULTI MODEL SELECTION TEST\n"]
    serial = 1
    for model, model_offset in ((1, 0.0), (2, 100.0)):
        lines.append(f"MODEL     {model:4d}\n")
        for chain, chain_offset in (("A", 0.0), ("B", 20.0)):
            for residue_index, residue_name in ((1, "ALA"), (2, "GLY")):
                origin = model_offset + chain_offset + (residue_index - 1) * 3.8
                for atom_name, delta, element in (
                    ("N", -1.3, "N"),
                    ("CA", 0.0, "C"),
                    ("C", 1.3, "C"),
                    ("O", 1.3, "O"),
                ):
                    y = 1.2 if atom_name == "O" else 0.0
                    lines.append(
                        f"ATOM  {serial:5d} {atom_name:^4s} {residue_name:>3s} {chain}{residue_index:4d}    "
                        f"{origin + delta:8.3f}{y:8.3f}{0.0:8.3f}{1.0:6.2f}{10.0:6.2f}          "
                        f"{element:>2s}\n"
                    )
                    serial += 1
            lines.append("TER\n")
        lines.append("ENDMDL\n")
    lines.append("END\n")
    return "".join(lines).encode("ascii")
