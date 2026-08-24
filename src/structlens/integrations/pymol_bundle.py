"""Versioned, data-only ``.structlens-pymol`` interchange bundles."""

from __future__ import annotations

import gzip
import json
import os
import tempfile
import zipfile
from collections.abc import Mapping, Sequence
from dataclasses import asdict
from pathlib import Path
from typing import Any

from structlens import __version__
from structlens.core.errors import (
    BundleCompatibilityError,
    BundleValidationError,
    UnsafeBundleError,
)
from structlens.core.models import (
    AnalysisResult,
    ComparisonMode,
    ProteinChain,
    ProteinStructure,
    ReferenceVsManyAnalysis,
    ResidueId,
)

BUNDLE_FORMAT = "structlens-pymol"
BUNDLE_SCHEMA_VERSION = "1.0"
SUPPORTED_BUNDLE_SCHEMA_MAJOR = 1


def write_pymol_bundle(
    output_path: str | Path,
    *,
    reference: ProteinStructure | ProteinChain,
    targets: Mapping[str, ProteinStructure | ProteinChain]
    | Sequence[ProteinStructure | ProteinChain],
    analysis: AnalysisResult | ReferenceVsManyAnalysis,
    provenance: Mapping[str, str] | None = None,
    visualization: Mapping[str, Any] | None = None,
    msa_summary: Mapping[str, Any] | None = None,
    conservation: Mapping[str, Any] | None = None,
    interactions: Mapping[str, Any] | None = None,
    sites: Mapping[str, Any] | None = None,
    evidence: Mapping[str, Any] | None = None,
    vectors: Mapping[str, Any] | None = None,
) -> Path:
    """Write and validate a deterministic, atomic analysis bundle."""

    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    target_items = _target_items(targets)
    reference_id = _structure_id(reference)
    target_ids = tuple(identifier for identifier, _ in target_items)
    if not reference_id or any(not identifier for identifier in target_ids):
        raise BundleValidationError("Structure identifiers must not be empty")
    if reference_id in target_ids:
        raise BundleValidationError("Reference structure cannot also be a target")
    mode = (
        analysis.comparison_mode.value
        if isinstance(analysis, ReferenceVsManyAnalysis)
        else ComparisonMode.PAIRWISE.value
    )
    structure_payloads: dict[str, tuple[str, bytes]] = {}
    structure_payloads[reference_id] = _read_structure(reference)
    for target_id, target in target_items:
        structure_payloads[target_id] = _read_structure(target)
    analysis_files = _analysis_payloads(analysis)
    manifest: dict[str, Any] = {
        "format": BUNDLE_FORMAT,
        "schema_version": BUNDLE_SCHEMA_VERSION,
        "created_by": {"application": "StructLens", "version": __version__},
        "plugin_compatibility": {"minimum_version": "0.1.0"},
        "analysis_id": _analysis_id(reference_id, target_ids, mode),
        "comparison_mode": mode,
        "reference_id": reference_id,
        "target_ids": list(target_ids),
        "structures": {
            identifier: f"structures/{'reference' if identifier == reference_id else f'target_{index:04d}'}.{suffix}"
            for index, (identifier, (_, payload)) in enumerate(structure_payloads.items())
            for suffix in [_suffix_for_payload(identifier, payload, structure_payloads[identifier][0])]
        },
    }
    # Rebuild the structure paths using the source suffix, keeping IDs deterministic.
    structure_entries: dict[str, bytes] = {}
    for index, (identifier, (suffix, payload)) in enumerate(structure_payloads.items()):
        filename = "reference" if identifier == reference_id else f"target_{index:04d}"
        entry = f"structures/{filename}.{suffix}"
        manifest["structures"][identifier] = entry
        structure_entries[entry] = payload
    files: dict[str, bytes] = {
        "manifest.json": _json_bytes(manifest),
        "provenance.json": _json_bytes(dict(provenance or {})),
        "analysis/summary.json": _json_bytes(_summary_payload(analysis)),
        "analysis/correspondence.json": _json_bytes(analysis_files["correspondence"]),
        "analysis/mutations.json": _json_bytes(analysis_files["mutations"]),
        "analysis/key_residues.json": _json_bytes(analysis_files["key_residues"]),
        "analysis/outliers.json": _json_bytes(analysis_files["outliers"]),
        "analysis/metrics.json": _json_bytes(analysis_files["metrics"]),
        "transforms/transforms.json": _json_bytes(analysis_files["transforms"]),
        "visualization/default_view.json": _json_bytes(dict(visualization or {})),
        "visualization/presets.json": _json_bytes(_default_presets()),
        **structure_entries,
    }
    optional_payloads = {
        "analysis/msa_summary.json": msa_summary,
        "analysis/conservation.json": conservation,
        "analysis/interactions.json": interactions,
        "analysis/sites.json": sites,
        "analysis/evidence.json": evidence,
        "visualization/vectors.json": vectors,
    }
    files.update({name: _json_bytes(dict(payload)) for name, payload in optional_payloads.items() if payload is not None})
    _write_atomic(destination, files)
    validate_pymol_bundle(destination)
    return destination


def validate_pymol_bundle(path: str | Path) -> dict[str, Any]:
    """Validate ZIP safety, schema compatibility, and referenced data files."""

    bundle_path = Path(path)
    if not bundle_path.is_file():
        raise BundleValidationError(f"Bundle does not exist: {bundle_path}")
    try:
        with zipfile.ZipFile(bundle_path) as archive:
            _validate_entries(archive)
            try:
                manifest = json.loads(archive.read("manifest.json"))
            except (KeyError, json.JSONDecodeError) as error:
                raise BundleValidationError("Bundle manifest.json is missing or malformed") from error
            if not isinstance(manifest, dict):
                raise BundleValidationError("Bundle manifest must be a JSON object")
            _validate_manifest(manifest)
            for identifier, entry in manifest["structures"].items():
                if entry not in archive.namelist():
                    raise BundleValidationError(
                        f"Structure '{identifier}' references missing bundle entry '{entry}'"
                    )
            for required in (
                "provenance.json",
                "analysis/summary.json",
                "analysis/correspondence.json",
                "analysis/mutations.json",
                "analysis/key_residues.json",
                "analysis/outliers.json",
                "analysis/metrics.json",
                "transforms/transforms.json",
            ):
                if required not in archive.namelist():
                    raise BundleValidationError(f"Bundle is missing required entry '{required}'")
            return manifest
    except zipfile.BadZipFile as error:
        raise BundleValidationError("The .structlens-pymol file is not a valid ZIP archive") from error


def _validate_entries(archive: zipfile.ZipFile) -> None:
    names = archive.namelist()
    if len(names) != len(set(names)):
        raise UnsafeBundleError("Bundle contains duplicate ZIP entries")
    for name in names:
        path = Path(name)
        if path.is_absolute() or ".." in path.parts or "\\" in name:
            raise UnsafeBundleError(f"Unsafe bundle path: {name}")
        if path.suffix.lower() in {".py", ".pyc", ".exe", ".dll", ".so", ".sh", ".bat"}:
            raise UnsafeBundleError(f"Executable payloads are not allowed: {name}")


def _validate_manifest(manifest: Mapping[str, Any]) -> None:
    if manifest.get("format") != BUNDLE_FORMAT:
        raise BundleValidationError("Unsupported StructLens bundle format")
    schema = str(manifest.get("schema_version", ""))
    try:
        major = int(schema.split(".", 1)[0])
    except (ValueError, IndexError) as error:
        raise BundleValidationError("Bundle schema_version must be MAJOR.MINOR") from error
    if major > SUPPORTED_BUNDLE_SCHEMA_MAJOR:
        raise BundleCompatibilityError(f"Bundle schema major {major} is newer than supported major 1")
    required = {"analysis_id", "comparison_mode", "reference_id", "target_ids", "structures"}
    missing = sorted(required - set(manifest))
    if missing:
        raise BundleValidationError(f"Bundle manifest is missing: {', '.join(missing)}")
    target_ids = manifest["target_ids"]
    structures = manifest["structures"]
    if not isinstance(target_ids, list) or not isinstance(structures, dict):
        raise BundleValidationError("Manifest target_ids and structures have invalid types")
    if len(target_ids) != len(set(target_ids)):
        raise BundleValidationError("Manifest target_ids must be unique")
    if manifest["reference_id"] in target_ids:
        raise BundleValidationError("Manifest reference_id cannot also be a target")
    if len(set(structures)) != len(structures):
        raise BundleValidationError("Manifest contains duplicate structure identifiers")
    if manifest["reference_id"] not in structures or any(target not in structures for target in target_ids):
        raise BundleValidationError("Manifest references an unknown structure identifier")


def _write_atomic(destination: Path, files: Mapping[str, bytes]) -> None:
    with tempfile.NamedTemporaryFile(prefix=".structlens-", suffix=".tmp", dir=destination.parent, delete=False) as handle:
        temporary = Path(handle.name)
    try:
        with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for name in sorted(files):
                info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
                info.compress_type = zipfile.ZIP_DEFLATED
                info.external_attr = 0o600 << 16
                archive.writestr(info, files[name])
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def _read_structure(structure: ProteinStructure | ProteinChain) -> tuple[str, bytes]:
    source = structure.source_path
    if source is None and isinstance(structure, ProteinStructure):
        source = structure.source_path
    if source is None and isinstance(structure, ProteinChain):
        source = structure.source_path
    if source is None:
        raise BundleValidationError(f"Structure '{_structure_id(structure)}' has no source file to export")
    path = Path(source)
    if not path.is_file():
        raise BundleValidationError(f"Structure source does not exist: {path}")
    if path.suffix.lower() == ".gz":
        suffix = path.name[:-3].rsplit(".", 1)[-1].lower()
        with gzip.open(path, "rb") as handle:
            return suffix, handle.read()
    return path.suffix.lower().lstrip(".") or "pdb", path.read_bytes()


def _target_items(
    targets: Mapping[str, ProteinStructure | ProteinChain]
    | Sequence[ProteinStructure | ProteinChain],
) -> tuple[tuple[str, ProteinStructure | ProteinChain], ...]:
    if isinstance(targets, Mapping):
        items = tuple((str(identifier), target) for identifier, target in targets.items())
    else:
        items = tuple((_structure_id(target), target) for target in targets)
    if len({identifier for identifier, _ in items}) != len(items):
        raise BundleValidationError("Target structure identifiers must be unique")
    return items


def _structure_id(structure: ProteinStructure | ProteinChain) -> str:
    return structure.structure_id


def _suffix_for_payload(identifier: str, payload: bytes, source: str) -> str:
    return source.rsplit(".", 1)[-1].lower() or "pdb"


def _analysis_id(reference_id: str, target_ids: Sequence[str], mode: str) -> str:
    import hashlib

    seed = "|".join((reference_id, mode, *target_ids)).encode("utf-8")
    return hashlib.sha256(seed).hexdigest()[:16]


def _summary_payload(analysis: AnalysisResult | ReferenceVsManyAnalysis) -> dict[str, Any]:
    if isinstance(analysis, AnalysisResult):
        return {"comparison_mode": "pairwise", "target_ids": [analysis.target_id], "metrics": _result_metrics(analysis)}
    return {
        "comparison_mode": analysis.comparison_mode.value,
        "reference_id": analysis.reference_id,
        "target_ids": list(analysis.target_ids),
        "metrics": {target_id: _target_metrics(target) for target_id, target in analysis.targets.items()},
    }


def _analysis_payloads(analysis: AnalysisResult | ReferenceVsManyAnalysis) -> dict[str, Any]:
    results = [analysis] if isinstance(analysis, AnalysisResult) else list(analysis.targets.values())
    correspondences: list[dict[str, Any]] = []
    mutations: list[dict[str, Any]] = []
    key_residues: list[dict[str, Any]] = []
    outliers: list[dict[str, Any]] = []
    metrics: dict[str, Any] = {}
    transforms: dict[str, Any] = {}
    for result in results:
        target_id = result.target_id
        rows = result.correspondences if isinstance(result, AnalysisResult) else result.correspondence
        events = result.mutations
        correspondences.extend({"target_id": target_id, **_correspondence_payload(row)} for row in rows)
        mutations.extend({"target_id": target_id, **asdict(event)} for event in events)
        key_residues.extend(
            {"target_id": target_id, **_correspondence_payload(row)} for row in rows if row.is_key_residue
        )
        outliers.extend(
            {"target_id": target_id, **_correspondence_payload(row)} for row in rows if row.is_outlier
        )
        metrics[target_id] = _result_metrics(result) if isinstance(result, AnalysisResult) else _target_metrics(result)
        if isinstance(result, AnalysisResult):
            transform = result.transform
            transforms[target_id] = (
                {}
                if transform is None
                else {
                    "rotation": [list(row) for row in transform.rotation],
                    "translation": list(transform.translation),
                }
            )
        else:
            transforms[target_id] = asdict(result.transform)
    return {
        "correspondence": correspondences,
        "mutations": mutations,
        "key_residues": key_residues,
        "outliers": outliers,
        "metrics": metrics,
        "transforms": transforms,
    }


def _correspondence_payload(row: Any) -> dict[str, Any]:
    payload = asdict(row)
    payload["reference"] = _residue_payload(row.reference)
    payload["target"] = _residue_payload(row.target)
    payload["status"] = row.status.value
    return payload


def _residue_payload(residue: ResidueId | None) -> dict[str, Any] | None:
    return None if residue is None else asdict(residue)


def _result_metrics(result: AnalysisResult) -> dict[str, Any]:
    return {
        "sequence_identity": result.sequence_identity,
        "sequence_similarity": result.sequence_similarity,
        "sequence_coverage": result.sequence_coverage,
        "strict_ca_rmsd_angstrom": result.strict_rmsd_angstrom,
        "refined_ca_rmsd_angstrom": result.refined_rmsd_angstrom,
        "tm_score": result.tm_score,
    }


def _target_metrics(result: Any) -> dict[str, Any]:
    structural = result.structural_metrics
    return {
        "sequence_identity": result.sequence_metrics.identity,
        "sequence_similarity": result.sequence_metrics.similarity,
        "sequence_coverage": result.sequence_metrics.coverage,
        "strict_ca_rmsd_angstrom": None if structural is None else structural.strict_ca_rmsd_angstrom,
    }


def _default_presets() -> dict[str, Any]:
    return {
        "Overview": {"highlight": "none", "representation": "cartoon"},
        "Mutation Focus": {"highlight": "mutations", "representation": "sticks"},
        "Key Residues": {"highlight": "key_residues", "representation": "sticks"},
        "Mutated Key Residues": {"highlight": "mutated_key_residues", "representation": "sticks+spheres"},
        "Structural Deviation": {"highlight": "high_ca_displacement", "representation": "sticks"},
        "Outliers": {"highlight": "outliers", "representation": "sticks+spheres"},
        "Publication": {"highlight": "key_residues", "representation": "cartoon+sticks"},
    }


def _json_bytes(payload: Any) -> bytes:
    return (json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")


__all__ = [
    "BUNDLE_FORMAT",
    "BUNDLE_SCHEMA_VERSION",
    "SUPPORTED_BUNDLE_SCHEMA_MAJOR",
    "validate_pymol_bundle",
    "write_pymol_bundle",
]
