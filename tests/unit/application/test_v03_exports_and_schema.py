from __future__ import annotations

from pathlib import Path

from openpyxl import load_workbook

from structlens.application.export_service import export_v03_xlsx
from structlens.application.project_state import ProjectState
from structlens.core.evidence import EvidenceCard, quality_for_sections
from structlens.core.models import ResidueId


def test_v03_project_schema_round_trip_preserves_authoritative_metadata() -> None:
    state = ProjectState(
        msa_settings={"algorithm": "muscle5"},
        interaction_thresholds={"hbond_distance_angstrom": 3.5},
        site_definitions=({"site_id": "active", "mode": "key_residues"},),
        distance_difference_metadata={"valid_pair_count": 4},
        evidence_sources={"A:10": ["msa", "structure"]},
    )
    restored = ProjectState.from_json(state.to_json())
    assert restored.schema_version == "3.0"
    assert restored.msa_settings["algorithm"] == "muscle5"
    assert restored.site_definitions[0]["site_id"] == "active"


def test_v03_xlsx_has_required_sheets_and_blank_unavailable(tmp_path: Path) -> None:
    path = tmp_path / "results.xlsx"
    card = EvidenceCard(ResidueId("ref", "1", "A", "10", None, "ALA"), "target", quality=quality_for_sections(["sequence"], ["structure"]))
    export_v03_xlsx(path, evidence_cards=(card,), provenance=("fixture",))
    workbook = load_workbook(path)
    assert {"MSA", "Conservation", "Amino Acid Frequencies", "Interaction Differences", "Sites", "Site Metrics", "Residue Evidence", "Evidence Quality"}.issubset(workbook.sheetnames)
    assert workbook["Residue Evidence"]["E2"].value is None
