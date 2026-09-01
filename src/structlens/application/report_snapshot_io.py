"""Evidence persistence and role-specific binding for canonical reports."""

from __future__ import annotations

import hashlib
import json
import os
import re
import zlib
from collections.abc import Mapping, Sequence
from pathlib import Path

from structlens.core.parsing import InputSelection, SourceSnapshot
from structlens.core.parsing.limits import SnapshotError
from structlens.core.parsing.snapshot import _GzipDecoder
from structlens.core.reports import AnalysisReport
from structlens.core.reports.safe_io import atomic_write_bytes, read_bounded_bytes, validate_safe_ancestors

_SHA256 = re.compile(r"^[0-9a-f]{64}$")
MAX_SNAPSHOT_BYTES = 100 * 1024 * 1024
_METADATA_FIELDS = frozenset(
    {
        "content_id",
        "raw_sha256",
        "decompressed_sha256",
        "display_name",
        "logical_format",
        "is_gzip",
        "storage_encoding",
    }
)


def _digest(value: str, name: str) -> str:
    if not isinstance(value, str) or not _SHA256.fullmatch(value.lower()):
        raise ValueError(f"{name} must be a lowercase SHA-256 hexadecimal digest")
    return value.lower()


def snapshot_path(content_id: str, directory: str | Path, raw_sha256: str | None = None) -> Path:
    """Return a legacy logical path or an exact logical-plus-raw path."""

    logical = _digest(content_id, "content_id")
    suffix = "" if raw_sha256 is None else f".{_digest(raw_sha256, 'raw_sha256')}"
    return Path(directory) / f"{logical}{suffix}.snapshot"


def _metadata_path(target: Path) -> Path:
    return target.with_suffix(".json")


def _stored_bytes(snapshot: SourceSnapshot) -> tuple[bytes, str]:
    if snapshot.raw_bytes is not None:
        return snapshot.raw_bytes, "raw"
    return snapshot.decompressed_bytes, "decompressed"


def _metadata_bytes(snapshot: SourceSnapshot, encoding: str) -> bytes:
    return json.dumps(
        {
            "content_id": snapshot.content_id,
            "raw_sha256": snapshot.raw_sha256,
            "decompressed_sha256": snapshot.decompressed_sha256,
            "display_name": snapshot.display_name,
            "logical_format": snapshot.logical_format,
            "is_gzip": snapshot.is_gzip,
            "storage_encoding": encoding,
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")


def _cleanup_partial(paths: tuple[Path, ...]) -> None:
    for artifact in paths:
        try:
            info = artifact.lstat()
            if artifact.is_symlink() or not artifact.is_file() or info.st_nlink != 1:
                continue
            os.unlink(artifact)
        except FileNotFoundError:
            pass


def _write_artifact(snapshot: SourceSnapshot, target: Path) -> None:
    stored, encoding = _stored_bytes(snapshot)
    metadata = _metadata_path(target)
    try:
        atomic_write_bytes(target, stored, max_bytes=MAX_SNAPSHOT_BYTES, label="snapshot", idempotent=True)
        atomic_write_bytes(
            metadata,
            _metadata_bytes(snapshot, encoding),
            max_bytes=MAX_SNAPSHOT_BYTES,
            label="snapshot metadata",
            idempotent=True,
        )
    except Exception:
        _cleanup_partial((target, metadata))
        raise


def _artifact_candidates(content_id: str, root: Path) -> tuple[Path, ...]:
    legacy = snapshot_path(content_id, root)
    candidates: list[Path] = [legacy] if legacy.exists() or _metadata_path(legacy).exists() else []
    for candidate in root.glob(f"{content_id}.*.snapshot"):
        if candidate not in candidates:
            candidates.append(candidate)
        if len(candidates) > 3:
            break
    return tuple(candidates)


def store_snapshot(snapshot: SourceSnapshot, directory: str | Path) -> Path:
    if not isinstance(snapshot, SourceSnapshot):
        raise TypeError("snapshot must be a SourceSnapshot")
    raw_size = len(snapshot.raw_bytes) if snapshot.raw_bytes is not None else 0
    if max(len(snapshot.decompressed_bytes), raw_size) > MAX_SNAPSHOT_BYTES:
        raise ValueError("snapshot exceeds the maximum supported size")
    root = Path(directory)
    validate_safe_ancestors(root, label="snapshot directory")
    if root.exists() and root.is_symlink():
        raise ValueError("snapshot directory must not be a link")
    root.mkdir(parents=True, exist_ok=True)
    validate_safe_ancestors(root, label="snapshot directory")

    legacy = snapshot_path(snapshot.content_id, root)
    exact = snapshot_path(snapshot.content_id, root, snapshot.raw_sha256)
    exact_metadata = _metadata_path(exact)
    if exact.exists() or exact_metadata.exists():
        existing = _load_artifact(exact, snapshot.content_id, snapshot.raw_sha256)
        if existing != snapshot:
            raise ValueError("existing persisted snapshot conflicts with the requested artifact")
        return exact

    if legacy.exists() or _metadata_path(legacy).exists():
        try:
            existing = _load_artifact(legacy, snapshot.content_id, None)
        except (OSError, ValueError) as exc:
            raise ValueError("existing persisted snapshot is invalid or conflicting") from exc
        if existing == snapshot:
            return legacy
        if existing.raw_sha256 == snapshot.raw_sha256 or snapshot.raw_bytes is None:
            raise ValueError("existing persisted snapshot conflicts with the requested artifact")
        _write_artifact(existing, snapshot_path(existing.content_id, root, existing.raw_sha256))
        _write_artifact(snapshot, exact)
        return exact

    if _artifact_candidates(snapshot.content_id, root):
        _write_artifact(snapshot, exact)
        return exact
    _write_artifact(snapshot, legacy)
    return legacy


def _logical_bytes(stored: bytes, *, is_gzip: bool, encoding: str) -> bytes:
    if encoding == "decompressed" or not is_gzip:
        return stored
    output = bytearray()
    decoder = _GzipDecoder(MAX_SNAPSHOT_BYTES)
    try:
        for start in range(0, len(stored), 64 * 1024):
            decoder.feed(stored[start : start + 64 * 1024], output)
        decoder.finish(output)
    except (SnapshotError, zlib.error) as exc:
        raise ValueError("persisted raw snapshot is not valid bounded gzip evidence") from exc
    return bytes(output)


def _load_artifact(target: Path, content_id: str, expected_raw: str | None) -> SourceSnapshot:
    metadata_path = _metadata_path(target)
    if not target.is_file() or not metadata_path.is_file() or target.is_symlink() or metadata_path.is_symlink():
        raise ValueError(f"persisted snapshot is missing for {content_id}")
    stored = read_bounded_bytes(target, max_bytes=MAX_SNAPSHOT_BYTES, label="persisted snapshot")
    try:
        metadata = json.loads(
            read_bounded_bytes(metadata_path, max_bytes=MAX_SNAPSHOT_BYTES, label="snapshot metadata")
        )
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError("persisted snapshot metadata is invalid") from exc
    if not isinstance(metadata, Mapping):
        raise ValueError("persisted snapshot metadata must be an object")
    unknown = set(metadata).difference(_METADATA_FIELDS)
    required = _METADATA_FIELDS.difference({"storage_encoding"})
    if unknown or required.difference(metadata):
        raise ValueError("persisted snapshot metadata is invalid")
    try:
        fields = ("content_id", "raw_sha256", "decompressed_sha256", "display_name", "logical_format")
        if any(not isinstance(metadata[field], str) for field in fields) or not isinstance(metadata["is_gzip"], bool):
            raise ValueError("persisted snapshot metadata is invalid")
        logical_id = _digest(str(metadata["content_id"]), "content_id")
        raw_id = _digest(str(metadata["raw_sha256"]), "raw_sha256")
        if logical_id != content_id or (expected_raw is not None and raw_id != expected_raw):
            raise ValueError("persisted snapshot content ID or raw identity is invalid")
        encoding = metadata.get("storage_encoding", "decompressed")
        if encoding not in {"raw", "decompressed"}:
            raise ValueError("persisted snapshot metadata is invalid")
        if encoding == "raw" and hashlib.sha256(stored).hexdigest() != raw_id:
            raise ValueError("persisted raw snapshot hash mismatch")
        logical = _logical_bytes(stored, is_gzip=bool(metadata["is_gzip"]), encoding=str(encoding))
        if hashlib.sha256(logical).hexdigest() != logical_id:
            raise ValueError("persisted snapshot hash mismatch")
        raw_bytes = stored if encoding == "raw" else None
        return SourceSnapshot(
            logical,
            raw_id,
            logical_id,
            str(metadata["decompressed_sha256"]),
            str(metadata["display_name"]),
            str(metadata["logical_format"]),
            bool(metadata["is_gzip"]),
            raw_bytes,
        )
    except (KeyError, TypeError, ValueError) as exc:
        if isinstance(exc, ValueError) and "snapshot" in str(exc):
            raise
        raise ValueError("persisted snapshot metadata is invalid") from exc


def load_snapshot(content_id: str, directory: str | Path, raw_sha256: str | None = None) -> SourceSnapshot:
    root = Path(directory)
    validate_safe_ancestors(root, label="snapshot directory")
    if root.is_symlink():
        raise ValueError("snapshot directory must not be a link")
    logical = _digest(content_id, "content_id")
    expected_raw = None if raw_sha256 is None else _digest(raw_sha256, "raw_sha256")
    if expected_raw is not None:
        exact = snapshot_path(logical, root, expected_raw)
        if exact.exists() or _metadata_path(exact).exists():
            return _load_artifact(exact, logical, expected_raw)
        legacy = snapshot_path(logical, root)
        if legacy.exists() or _metadata_path(legacy).exists():
            return _load_artifact(legacy, logical, expected_raw)
        raise ValueError(f"persisted snapshot is missing for {logical}")

    candidates = _artifact_candidates(logical, root)
    if not candidates:
        raise ValueError(f"persisted snapshot is missing for {logical}")
    loaded = tuple(_load_artifact(path, logical, None) for path in candidates)
    by_raw = {item.raw_sha256: item for item in loaded}
    if len(by_raw) != 1:
        raise ValueError(f"persisted snapshot is ambiguous; raw SHA-256 is required for {logical}")
    return next(iter(by_raw.values()))


def verify_snapshot(snapshot: SourceSnapshot, selection: InputSelection) -> None:
    if not isinstance(snapshot, SourceSnapshot) or not isinstance(selection, InputSelection):
        raise TypeError("snapshot and selection must be typed values")
    if snapshot.content_id != selection.content_id or snapshot.logical_format != selection.format.value:
        raise ValueError("source snapshot does not match the report selection content")


def _capture_source_path(path: str | Path) -> SourceSnapshot:
    source = Path(path)
    validate_safe_ancestors(source, label="source path")
    try:
        info = source.lstat()
    except FileNotFoundError:
        raise
    if source.is_symlink() or not source.is_file() or info.st_nlink != 1:
        raise ValueError("source path must be a regular, non-link file")
    return SourceSnapshot.from_path(source)


def verify_source_path(path: str | Path, selection: InputSelection) -> SourceSnapshot:
    snapshot = _capture_source_path(path)
    verify_snapshot(snapshot, selection)
    return snapshot


def _provenance_hash(hashes: Mapping[str, str], role: str, kind: str) -> str | None:
    for key, value in hashes.items():
        if key.casefold().replace("-", "_") in {f"{role}_{kind}", f"{role}_{kind}_sha256", f"{role}_source_{kind}"}:
            return value
    return None


def verify_report_inputs(
    report: AnalysisReport,
    snapshots: Sequence[SourceSnapshot],
    source_paths: Sequence[str | Path],
    snapshot_dir: str | Path | None,
) -> None:
    if len(snapshots) > 2 or len(source_paths) > 2:
        raise ValueError("at most reference and target source evidence may be supplied")
    if not snapshots and not source_paths and snapshot_dir is None:
        raise ValueError("complete reference and target source evidence is required")
    roles = ("reference", "target")
    selections = (report.reference_selection, report.target_selection)
    hashes = {} if report.provenance is None else report.provenance.input_hashes
    expected: dict[str, tuple[str, str]] = {}
    for role, selection in zip(roles, selections, strict=True):
        content_hash = _provenance_hash(hashes, role, "content")
        raw_hash = _provenance_hash(hashes, role, "raw")
        if content_hash is None or raw_hash is None:
            raise ValueError(f"report provenance is missing {role} raw/content evidence hashes")
        expected[role] = (_digest(content_hash, f"{role} content hash"), _digest(raw_hash, f"{role} raw hash"))
        if expected[role][0] != selection.content_id:
            raise ValueError(f"{role} report selection does not match report provenance content hash")

    evidence: dict[str, SourceSnapshot] = {}
    for snapshot_index, snapshot in enumerate(snapshots):
        if not isinstance(snapshot, SourceSnapshot):
            raise TypeError("snapshots must contain SourceSnapshot values")
        matching_roles = [role for role in roles if expected[role] == (snapshot.content_id, snapshot.raw_sha256)]
        if len(matching_roles) > 1:
            named_roles = [
                role for role in matching_roles if selections[roles.index(role)].display_name == snapshot.display_name
            ]
            if named_roles:
                matching_roles = named_roles
            elif snapshot_index < len(roles) and roles[snapshot_index] in matching_roles:
                matching_roles = [roles[snapshot_index]]
        if not matching_roles:
            if any(selection.content_id == snapshot.content_id for selection in selections):
                raise ValueError("multiple conflicting snapshots supplied for one source role")
            raise ValueError("source snapshot raw/content identity does not match a report provenance role")
        for role in matching_roles:
            if role in evidence and evidence[role] != snapshot:
                raise ValueError("multiple conflicting snapshots supplied for one source role")
            evidence[role] = snapshot
    for source_path in source_paths:
        captured = _capture_source_path(source_path)
        matching_roles = [role for role in roles if expected[role] == (captured.content_id, captured.raw_sha256)]
        if not matching_roles:
            raise ValueError("source path does not match report provenance raw/content hashes")
        for role in matching_roles:
            verify_snapshot(captured, selections[roles.index(role)])
            if role in evidence and evidence[role] != captured:
                raise ValueError("source path does not match the supplied typed snapshot")
            evidence[role] = captured
    if snapshot_dir is not None:
        for role, selection in zip(roles, selections, strict=True):
            stored = load_snapshot(selection.content_id, snapshot_dir, expected[role][1])
            if role in evidence and evidence[role] != stored:
                raise ValueError("persisted snapshot conflicts with supplied source evidence")
            evidence[role] = stored
    if set(evidence) != set(roles):
        raise ValueError("complete reference and target source evidence is required")
    for role, selection in zip(roles, selections, strict=True):
        snapshot = evidence[role]
        verify_snapshot(snapshot, selection)
        if (snapshot.content_id, snapshot.raw_sha256) != expected[role]:
            raise ValueError(f"{role} snapshot does not match report provenance raw/content hashes")


__all__ = [
    "MAX_SNAPSHOT_BYTES",
    "load_snapshot",
    "snapshot_path",
    "store_snapshot",
    "verify_report_inputs",
    "verify_snapshot",
    "verify_source_path",
]
