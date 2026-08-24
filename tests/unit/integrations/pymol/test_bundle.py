import zipfile
from pathlib import Path

import pytest

from structlens.core.models import AnalysisResult, ProteinChain
from structlens.integrations.pymol_bundle import validate_pymol_bundle, write_pymol_bundle


def _empty_result() -> AnalysisResult:
    return AnalysisResult(
        reference_id="ref",
        target_id="target",
        correspondences=(),
        mutations=(),
        sequence_identity=0.0,
        sequence_coverage=0.0,
        alignment_decision="sequence-guided",
    )


def test_bundle_round_trip_is_data_only(tmp_path: Path) -> None:
    ref_path = tmp_path / "ref.pdb"
    target_path = tmp_path / "target.pdb"
    ref_path.write_text("ATOM\n", encoding="utf-8")
    target_path.write_text("ATOM\n", encoding="utf-8")
    reference = ProteinChain("ref", "1", "A", source_path=str(ref_path))
    target = ProteinChain("target", "1", "A", source_path=str(target_path))
    bundle = write_pymol_bundle(
        tmp_path / "analysis.structlens-pymol",
        reference=reference,
        targets=(target,),
        analysis=_empty_result(),
    )
    manifest = validate_pymol_bundle(bundle)
    assert manifest["schema_version"] == "1.0"
    with zipfile.ZipFile(bundle) as archive:
        assert all(not name.endswith((".py", ".exe")) for name in archive.namelist())


def test_bundle_rejects_path_traversal(tmp_path: Path) -> None:
    bundle = tmp_path / "unsafe.structlens-pymol"
    with zipfile.ZipFile(bundle, "w") as archive:
        archive.writestr("../evil.py", "raise SystemExit")
    with pytest.raises(Exception, match="Unsafe|path"):
        validate_pymol_bundle(bundle)
