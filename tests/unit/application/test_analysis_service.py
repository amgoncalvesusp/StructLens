from __future__ import annotations

from threading import Event

import numpy as np
import pytest

from structlens.application.analysis_service import AnalysisService
from structlens.core.errors import AnalysisCancelledError
from structlens.core.models import (
    AlignmentMode,
    AnalysisSettings,
    AtomRecord,
    ProteinChain,
    ProteinStructure,
    ResidueId,
    ResidueNumbering,
    ResidueRecord,
)


def _structure(structure_id: str, offset: float = 0.0) -> ProteinStructure:
    records = []
    for index, (name, one_letter) in enumerate((("ALA", "A"), ("GLY", "G")), 1):
        rid = ResidueId(structure_id, "1", "A", str(index), None, name)
        atoms = tuple(
            AtomRecord(atom, "C", (index + offset, i + offset, 0.0))
            for i, atom in enumerate(("N", "CA", "C", "O"))
        )
        records.append(
            ResidueRecord(
                rid,
                ResidueNumbering(str(index), str(index), None),
                name,
                one_letter,
                atoms,
            )
        )
    chain = ProteinChain(
        structure_id,
        "1",
        "A",
        tuple(record.residue_id for record in records),
        "AG",
        tuple(records),
    )
    return ProteinStructure(structure_id, (chain,))


def test_analysis_service_builds_authoritative_map_and_zero_self_rmsd() -> None:
    result = AnalysisService().analyze(_structure("ref"), _structure("target"))

    assert result.alignment_decision.startswith("sequence-guided")
    assert result.sequence_identity == 1.0
    assert result.mutation_count == 0
    assert result.strict_rmsd_angstrom == 0.0
    assert len(result.correspondences) == 2
    assert all(
        correspondence.backbone_rmsd_angstrom == 0.0
        for correspondence in result.correspondences
    )


def test_analysis_service_is_rigid_body_invariant() -> None:
    target = _structure("target", offset=4.0)
    result = AnalysisService().analyze(
        _structure("ref"), target, AnalysisSettings(refined_rmsd=True)
    )

    assert np.isclose(result.strict_rmsd_angstrom, 0.0)
    assert result.refined_rmsd_angstrom == result.strict_rmsd_angstrom


def test_manual_mapping_coverage_uses_full_chain_denominators() -> None:
    reference = _structure("ref")
    target = _structure("target")
    reference_chain = reference.chains[0]
    target_chain = target.chains[0]

    result = AnalysisService().analyze(
        reference,
        target,
        AnalysisSettings(alignment_mode=AlignmentMode.MANUAL),
        manual_pairs=[(reference_chain.residue_records[0].residue_id, target_chain.residue_records[0].residue_id)],
    )

    assert result.sequence_identity == 1.0
    assert result.sequence_coverage == 0.5


def test_analysis_honors_cancellation_before_work_starts() -> None:
    cancel_event = Event()
    cancel_event.set()

    with pytest.raises(AnalysisCancelledError):
        AnalysisService().analyze(
            _structure("ref"),
            _structure("target"),
            cancel_event=cancel_event,
        )


def _structure_with_nonstandard(structure_id: str) -> ProteinStructure:
    """A chain carrying a selenomethionine, which has no canonical one-letter code."""

    residues = (("ALA", "A"), ("MSE", None), ("GLY", "G"))
    records = []
    for index, (name, one_letter) in enumerate(residues, 1):
        rid = ResidueId(structure_id, "1", "A", str(index), None, name)
        atoms = tuple(
            AtomRecord(atom, "C", (float(index), float(i), 0.0))
            for i, atom in enumerate(("N", "CA", "C", "O"))
        )
        records.append(
            ResidueRecord(
                rid,
                ResidueNumbering(str(index), str(index), None),
                name,
                one_letter,
                atoms,
                is_standard=one_letter is not None,
            )
        )
    chain = ProteinChain(
        structure_id,
        "1",
        "A",
        tuple(record.residue_id for record in records),
        "".join(record.one_letter or "X" for record in records),
        tuple(records),
    )
    return ProteinStructure(structure_id, (chain,))


def test_nonstandard_residue_is_not_counted_as_a_self_mutation() -> None:
    """Comparing a structure to itself must report no mutations, MSE included."""

    result = AnalysisService().analyze(
        _structure_with_nonstandard("ref"), _structure_with_nonstandard("target")
    )

    assert result.mutation_count == 0
    nonstandard = [event for event in result.mutations if event.kind.value == "nonstandard"]
    assert nonstandard, "the MSE position must still be reported as non-standard"


def test_sequence_comparison_stores_the_authoritative_transform() -> None:
    """Site metrics need the transform that produced the structural numbers."""

    result = AnalysisService().analyze(_structure("ref"), _structure("target", offset=4.0))

    assert result.transform is not None
    rotation = np.asarray(result.transform.rotation, dtype=float)
    translation = np.asarray(result.transform.translation, dtype=float)
    assert rotation.shape == (3, 3)
    assert np.linalg.det(rotation) == pytest.approx(1.0)

    reference = _structure("ref").chains[0].residue_records
    target = _structure("target", offset=4.0).chains[0].residue_records
    reference_ca = np.asarray(
        [next(a.coordinate for a in r.atoms if a.name == "CA") for r in reference], dtype=float
    )
    target_ca = np.asarray(
        [next(a.coordinate for a in r.atoms if a.name == "CA") for r in target], dtype=float
    )
    # The stored transform must be the one that actually places the target onto
    # the reference, not an unrelated or identity placeholder.
    fitted = target_ca @ rotation + translation
    assert np.allclose(fitted, reference_ca, atol=1e-6)
