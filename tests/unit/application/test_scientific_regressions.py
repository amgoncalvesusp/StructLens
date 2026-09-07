from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from zipfile import ZipFile

import numpy as np
import pytest

from structlens.application.analysis_service import AnalysisService
from structlens.core.alignment.superposition import superpose
from structlens.core.geometry.kabsch import apply_transform
from structlens.core.models import (
    AlignmentMode,
    AnalysisSettings,
    AtomRecord,
    ProteinChain,
    ResidueId,
    ResidueNumbering,
    ResidueRecord,
    StructuralTransform,
)
from structlens.integrations.pymol_bundle import write_pymol_bundle


def _chain(identifier: str, sequence: str, coordinates: np.ndarray | None = None) -> ProteinChain:
    records = tuple(
        ResidueRecord(
            ResidueId(identifier, "1", "A", str(index + 1), None, letter),
            ResidueNumbering(str(index + 1), str(index + 1), None),
            letter,
            letter,
            (AtomRecord("CA", "C", tuple(coordinates[index])),) if coordinates is not None else (),
        )
        for index, letter in enumerate(sequence)
    )
    return ProteinChain(identifier, "1", "A", tuple(record.residue_id for record in records), sequence, records)


def test_multiple_alignment_anchors_on_reference_residue_after_target_insertion() -> None:
    reference = _chain("ref", "ACDEFGHIK")
    same = _chain("a", "ACDEFGHIK")
    inserted = _chain("b", "ACWDEFGHIK")
    analysis = AnalysisService().analyze_multiple_structure_alignment(
        reference, {"a": same, "b": inserted}, AnalysisSettings(alignment_mode=AlignmentMode.SEQUENCE)
    )

    for position, ref_residue in zip(
        (position for position in analysis.aligned_positions if position.reference_residue is not None),
        reference.residues,
        strict=True,
    ):
        assert position.reference_residue == ref_residue
        assert position.residues["a"].residue_name == ref_residue.residue_name
        assert position.residues["b"].residue_name == ref_residue.residue_name
        assert position.coverage == 1.0
    insertion = next(position for position in analysis.aligned_positions if position.reference_residue is None)
    assert insertion.residues["a"] is None
    assert insertion.residues["b"] == inserted.residues[2]
    assert insertion.coverage == pytest.approx(1 / 3)
    assert insertion.per_structure_deviation_angstrom["ref"] is None
    assert insertion.ca_positional_variability_angstrom is None


@pytest.mark.parametrize("sequence", ["WACDEFGHIK", "ACWDEFGHIK", "ACDEFGHIKW"])
def test_multiple_alignment_does_not_claim_equivalence_between_target_insertions(sequence: str) -> None:
    analysis = AnalysisService().analyze_multiple_structure_alignment(
        _chain("ref", "ACDEFGHIK"),
        {"a": _chain("a", sequence), "b": _chain("b", sequence)},
        AnalysisSettings(alignment_mode=AlignmentMode.SEQUENCE),
    )
    insertions = [position for position in analysis.aligned_positions if position.reference_residue is None]
    assert len(insertions) == 2
    assert [position.alignment_index for position in analysis.aligned_positions] == list(range(11))
    for position in insertions:
        assert sum(residue is not None for residue in position.residues.values()) == 1
        assert position.coverage == pytest.approx(1 / 3)


def test_multiple_alignment_keeps_reference_positions_across_target_deletions() -> None:
    reference = _chain("ref", "ACDEFGHIK")
    analysis = AnalysisService().analyze_multiple_structure_alignment(
        reference,
        {"deleted": _chain("deleted", "ACEFGHIK"), "same": _chain("same", reference.sequence)},
        AnalysisSettings(alignment_mode=AlignmentMode.SEQUENCE),
    )
    assert tuple(position.reference_residue for position in analysis.aligned_positions) == reference.residues
    deletion = analysis.aligned_positions[2]
    assert deletion.residues["deleted"] is None
    assert deletion.residues["same"].residue_name == "D"
    assert deletion.coverage == pytest.approx(2 / 3)


@pytest.mark.parametrize("cutoff", [1.0, 0.000001])
def test_refinement_final_fit_matches_reported_inliers_after_one_iteration(cutoff: float) -> None:
    reference = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [20, 0, 0]], dtype=float)
    target = reference + np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 10, 0]])
    result = AnalysisService().analyze(
        _chain("ref", "ACDE", reference),
        _chain("target", "ACDE", target),
        AnalysisSettings(refined_rmsd=True, refinement_cutoff_angstrom=cutoff, refinement_max_iterations=1),
    )
    included = [index for index in range(4) if index not in result.excluded_alignment_indices]
    assert result.refined_residue_count == len(included)
    assert tuple(item.alignment_index for item in result.correspondences if item.is_outlier) == result.excluded_alignment_indices
    if included:
        expected = superpose(reference[included], target[included])
        assert result.refined_rmsd_angstrom == pytest.approx(expected.strict_rmsd_angstrom)
    else:
        assert result.refined_rmsd_angstrom is None


@pytest.mark.parametrize("mode", [AlignmentMode.SEQUENCE, AlignmentMode.STRUCTURE])
def test_result_and_bundle_transform_reproduce_the_reported_strict_frame(mode: AlignmentMode, tmp_path: Path) -> None:
    reference_coordinates = np.array([[0, 0, 0], [1, 0, 0], [2, 1, 0], [20, 0, 1]], dtype=float)
    rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], dtype=float)
    target_coordinates = reference_coordinates @ rotation + np.array([10, 20, 30])
    target_coordinates = target_coordinates + np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 5, 0]])
    chains = []
    for identifier, coordinates in (("ref", reference_coordinates), ("target", target_coordinates)):
        path = tmp_path / f"{identifier}.pdb"
        path.write_text("\n".join(
            f"ATOM  {index:5d}  CA  ALA A{index:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C"
            for index, (x, y, z) in enumerate(coordinates, 1)
        ) + "\nEND\n", encoding="ascii")
        chains.append(replace(_chain(identifier, "AAAA", coordinates), source_path=str(path)))
    reference, target = chains
    mapping = AnalysisService().analyze(reference, target).correspondences
    adapter_result = SimpleNamespace(
        correspondences=mapping, tm_score=0.9, executable_version="test", metadata={},
        transform=StructuralTransform(translation=(100.0, 200.0, 300.0)),
    )
    adapter = SimpleNamespace(align=lambda *_args: adapter_result)
    result = AnalysisService(adapter).analyze(
        reference, target, AnalysisSettings(alignment_mode=mode, refined_rmsd=True, refinement_max_iterations=1)
    )
    assert result.transform is not None
    fitted = apply_transform(target_coordinates, result.transform.rotation, result.transform.translation)
    distances = np.linalg.norm(reference_coordinates - fitted, axis=1)
    assert np.sqrt(np.mean(distances ** 2)) == pytest.approx(result.strict_rmsd_angstrom)
    assert distances == pytest.approx([row.ca_displacement_angstrom for row in result.correspondences])
    bundle = write_pymol_bundle(
        tmp_path / "analysis.structlens-pymol", reference=reference, targets={"target": target}, analysis=result
    )
    with ZipFile(bundle) as archive:
        exported = json.loads(archive.read("transforms/transforms.json"))["target"]
    exported_fit = apply_transform(target_coordinates, exported["rotation"], exported["translation"])
    assert exported_fit == pytest.approx(fitted)


def test_analysis_without_coordinates_does_not_fabricate_a_transform() -> None:
    result = AnalysisService().analyze(_chain("ref", "ACDE"), _chain("target", "ACDE"))
    assert result.transform is None
