"""Public-package smoke test for the foundational domain models."""

from structlens.core.models import (
    AlignmentMode,
    CorrespondenceStatus,
    MutationKind,
    ResidueId,
)


def test_domain_models_are_exported_from_core_models() -> None:
    residue = ResidueId("fixture", "1", "A", "1", None, "ALA")

    assert residue.residue_name == "ALA"
    assert AlignmentMode.AUTO.value == "auto"
    assert CorrespondenceStatus.CONSERVED.value == "conserved"
    assert MutationKind.SUBSTITUTION.value == "substitution"
