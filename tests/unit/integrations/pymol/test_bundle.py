import json
import zipfile
from pathlib import Path

import pytest

from structlens.core.errors import BundleValidationError
from structlens.core.metrics.sequence_metrics import SequenceAlignmentMetrics
from structlens.core.models import (
    AnalysisResult,
    ProteinChain,
    ReferenceVsManyAnalysis,
    TargetAnalysis,
)
from structlens.core.provenance import MethodProvenance
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


def test_bundle_preserves_full_typed_and_legacy_provenance(tmp_path: Path) -> None:
    ref_path = tmp_path / "ref.pdb"
    target_path = tmp_path / "target.pdb"
    ref_path.write_text("ATOM\n", encoding="utf-8")
    target_path.write_text("ATOM\n", encoding="utf-8")
    typed = MethodProvenance(
        "structlens.compare",
        "0.4.0",
        {"grid_spacing": 0.5, "nested": {"enabled": True}},
        {"grid_spacing": "angstrom"},
        {"scipy": "1.17.1"},
        {"target": "a" * 64},
    )
    analysis = AnalysisResult(
        reference_id="ref",
        target_id="target",
        correspondences=(),
        mutations=(),
        sequence_identity=1.0,
        sequence_coverage=1.0,
        alignment_decision="test",
        provenance={"backend": "legacy"},
        method_provenance=typed,
    )
    bundle = write_pymol_bundle(
        tmp_path / "typed.structlens-pymol",
        reference=ProteinChain("ref", "1", "A", source_path=str(ref_path)),
        targets=(ProteinChain("target", "1", "A", source_path=str(target_path)),),
        analysis=analysis,
    )

    with zipfile.ZipFile(bundle) as archive:
        payload = json.loads(archive.read("provenance.json"))
    assert payload["legacy"] == {"backend": "legacy"}
    assert payload["method_provenance"] == typed.to_json()


def test_reference_many_bundle_keys_typed_provenance_by_target(tmp_path: Path) -> None:
    ref_path = tmp_path / "ref.pdb"
    target_path = tmp_path / "target.pdb"
    ref_path.write_text("ATOM\n", encoding="utf-8")
    target_path.write_text("ATOM\n", encoding="utf-8")
    typed = MethodProvenance(
        "structlens.compare",
        "0.4.0",
        {},
        {},
        input_hashes={"target": "c" * 64},
    )
    target_analysis = TargetAnalysis(
        "target",
        (),
        (),
        SequenceAlignmentMetrics(1.0, 1.0, 1.0, 1.0, 1.0, 0),
        method_provenance=typed,
    )
    bundle = write_pymol_bundle(
        tmp_path / "many.structlens-pymol",
        reference=ProteinChain("ref", "1", "A", source_path=str(ref_path)),
        targets={"target": ProteinChain("target", "1", "A", source_path=str(target_path))},
        analysis=ReferenceVsManyAnalysis("ref", {"target": target_analysis}),
    )

    with zipfile.ZipFile(bundle) as archive:
        payload = json.loads(archive.read("provenance.json"))
    assert payload["method_provenance"]["targets"]["target"] == typed.to_json()


def test_bundle_validation_rejects_unknown_reference_many_provenance_target(tmp_path: Path) -> None:
    ref_path = tmp_path / "ref.pdb"
    target_path = tmp_path / "target.pdb"
    ref_path.write_text("ATOM\n", encoding="utf-8")
    target_path.write_text("ATOM\n", encoding="utf-8")
    typed = MethodProvenance(
        "structlens.compare",
        "0.4.0",
        {},
        {},
        input_hashes={"target": "c" * 64},
    )
    target_analysis = TargetAnalysis(
        "target",
        (),
        (),
        SequenceAlignmentMetrics(1.0, 1.0, 1.0, 1.0, 1.0, 0),
        method_provenance=typed,
    )
    bundle = write_pymol_bundle(
        tmp_path / "unknown-target.structlens-pymol",
        reference=ProteinChain("ref", "1", "A", source_path=str(ref_path)),
        targets={"target": ProteinChain("target", "1", "A", source_path=str(target_path))},
        analysis=ReferenceVsManyAnalysis("ref", {"target": target_analysis}),
    )
    rewritten = tmp_path / "unknown-target-rewritten.structlens-pymol"
    with zipfile.ZipFile(bundle) as source, zipfile.ZipFile(rewritten, "w") as destination:
        for name in source.namelist():
            content = source.read(name)
            if name == "provenance.json":
                payload = json.loads(content)
                method_payload = payload["method_provenance"]
                method_payload["targets"]["other_target"] = method_payload["targets"].pop("target")
                content = json.dumps(payload).encode("utf-8")
            destination.writestr(name, content)

    with pytest.raises(BundleValidationError, match="target"):
        validate_pymol_bundle(rewritten)


def test_bundle_rejects_conflicting_explicit_typed_provenance(tmp_path: Path) -> None:
    ref_path = tmp_path / "ref.pdb"
    target_path = tmp_path / "target.pdb"
    ref_path.write_text("ATOM\n", encoding="utf-8")
    target_path.write_text("ATOM\n", encoding="utf-8")
    embedded = MethodProvenance("structlens.compare", "0.4.0", {}, {}, input_hashes={"target": "a" * 64})
    explicit = MethodProvenance("structlens.compare", "0.4.0", {}, {}, input_hashes={"target": "b" * 64})
    analysis = AnalysisResult(
        "ref", "target", (), (), 1.0, 1.0, "test", method_provenance=embedded
    )

    with pytest.raises(BundleValidationError, match="conflict"):
        write_pymol_bundle(
            tmp_path / "conflict.structlens-pymol",
            reference=ProteinChain("ref", "1", "A", source_path=str(ref_path)),
            targets=(ProteinChain("target", "1", "A", source_path=str(target_path)),),
            analysis=analysis,
            method_provenance=explicit,
        )


def test_bundle_validation_rejects_tampered_typed_provenance(tmp_path: Path) -> None:
    ref_path = tmp_path / "ref.pdb"
    target_path = tmp_path / "target.pdb"
    ref_path.write_text("ATOM\n", encoding="utf-8")
    target_path.write_text("ATOM\n", encoding="utf-8")
    typed = MethodProvenance("structlens.compare", "0.4.0", {}, {}, input_hashes={"target": "a" * 64})
    bundle = write_pymol_bundle(
        tmp_path / "tampered.structlens-pymol",
        reference=ProteinChain("ref", "1", "A", source_path=str(ref_path)),
        targets=(ProteinChain("target", "1", "A", source_path=str(target_path)),),
        analysis=AnalysisResult("ref", "target", (), (), 1.0, 1.0, "test", method_provenance=typed),
    )
    rewritten = tmp_path / "tampered-rewritten.structlens-pymol"
    with zipfile.ZipFile(bundle) as source, zipfile.ZipFile(rewritten, "w") as destination:
        for name in source.namelist():
            content = source.read(name)
            if name == "provenance.json":
                payload = json.loads(content)
                payload["method_provenance"]["artifact_id"] = "0" * 64
                content = json.dumps(payload).encode("utf-8")
            destination.writestr(name, content)

    with pytest.raises(BundleValidationError, match="provenance"):
        validate_pymol_bundle(rewritten)


def test_bundle_validation_rejects_null_typed_provenance(tmp_path: Path) -> None:
    ref_path = tmp_path / "ref.pdb"
    target_path = tmp_path / "target.pdb"
    ref_path.write_text("ATOM\n", encoding="utf-8")
    target_path.write_text("ATOM\n", encoding="utf-8")
    typed = MethodProvenance("structlens.compare", "0.4.0", {}, {}, input_hashes={"target": "a" * 64})
    bundle = write_pymol_bundle(
        tmp_path / "null-typed.structlens-pymol",
        reference=ProteinChain("ref", "1", "A", source_path=str(ref_path)),
        targets=(ProteinChain("target", "1", "A", source_path=str(target_path)),),
        analysis=AnalysisResult("ref", "target", (), (), 1.0, 1.0, "test", method_provenance=typed),
    )
    rewritten = tmp_path / "null-typed-rewritten.structlens-pymol"
    with zipfile.ZipFile(bundle) as source, zipfile.ZipFile(rewritten, "w") as destination:
        for name in source.namelist():
            content = source.read(name)
            if name == "provenance.json":
                payload = json.loads(content)
                payload["method_provenance"] = None
                content = json.dumps(payload).encode("utf-8")
            destination.writestr(name, content)

    with pytest.raises(BundleValidationError, match="provenance"):
        validate_pymol_bundle(rewritten)
