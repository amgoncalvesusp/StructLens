from __future__ import annotations

import numpy as np

from structlens.application.analysis_service import AnalysisService
from structlens.core.models import (
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


def test_analysis_service_is_rigid_body_invariant() -> None:
    target = _structure("target", offset=4.0)
    result = AnalysisService().analyze(
        _structure("ref"), target, AnalysisSettings(refined_rmsd=True)
    )

    assert np.isclose(result.strict_rmsd_angstrom, 0.0)
    assert result.refined_rmsd_angstrom == result.strict_rmsd_angstrom
