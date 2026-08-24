from pathlib import Path

from structlens.application.analysis_service import AnalysisService
from structlens.core.models import ComparisonMode
from structlens.core.parsing import load_structure


def test_reference_vs_many_keeps_requested_target_ids() -> None:
    structure = load_structure(Path("tests/fixtures/parsing/numbering_altloc.pdb"))
    analysis = AnalysisService().analyze_reference_vs_many(
        structure,
        {"target_a": structure, "target_b": structure},
    )
    assert analysis.comparison_mode is ComparisonMode.REFERENCE_VS_MANY
    assert analysis.target_ids == ("target_a", "target_b")


def test_all_vs_all_calculates_unique_pairs_once() -> None:
    structure = load_structure(Path("tests/fixtures/parsing/numbering_altloc.pdb"))
    analysis = AnalysisService().analyze_all_vs_all(
        {"a": structure, "b": structure, "c": structure}
    )
    assert len(analysis.pair_results) == 3
    assert analysis.matrices["sequence_identity"].value("c", "a") == 1.0
