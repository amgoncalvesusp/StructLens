from structlens.core.models import ResidueId
from structlens.integrations.pymol.selections import residue_selection, selection_name


def test_selection_names_are_namespaced_and_safe() -> None:
    name = selection_name("Project A", "Target/1", "mutations")
    assert name.startswith("structlens_")
    assert "/" not in name
    assert " " not in name


def test_residue_selection_uses_chain_and_insertion_code() -> None:
    residue = ResidueId("target", "1", "B", "100", "A", "GLY")
    assert residue_selection(residue) == "chain B and resi 100A"
