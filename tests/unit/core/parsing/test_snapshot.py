"""Tests for bounded, content-addressed source snapshots."""

from __future__ import annotations

import gzip
import hashlib
from pathlib import Path

import pytest

from structlens.core.parsing.limits import ParseLimits
from structlens.core.parsing.snapshot import (
    SnapshotError,
    SnapshotLimitError,
    SourceSnapshot,
    capture_snapshot,
)


def test_parse_limits_are_frozen_and_require_positive_values() -> None:
    limits = ParseLimits(
        max_raw_bytes=10,
        max_decompressed_bytes=20,
        max_models=2,
        max_chains=3,
        max_residues=4,
        max_atoms=5,
    )

    assert limits.max_raw_bytes == 10
    with pytest.raises((AttributeError, TypeError)):
        limits.max_atoms = 6  # type: ignore[misc]
    with pytest.raises(ValueError, match="max_raw_bytes"):
        ParseLimits(max_raw_bytes=0)
    with pytest.raises(ValueError, match="max_atoms"):
        ParseLimits(max_atoms=-1)


def test_capture_snapshot_keeps_decompressed_bytes_and_hashes(tmp_path: Path) -> None:
    source = b"ATOM      1  CA  ALA A   1      0.000   0.000   0.000\n"
    path = tmp_path / "sample.pdb"
    path.write_bytes(source)

    snapshot = capture_snapshot(path)

    assert isinstance(snapshot, SourceSnapshot)
    assert snapshot.decompressed_bytes == source
    assert snapshot.content_id == hashlib.sha256(source).hexdigest()
    assert snapshot.decompressed_sha256 == snapshot.content_id
    assert snapshot.raw_sha256 == hashlib.sha256(source).hexdigest()
    assert snapshot.display_name == "sample.pdb"
    assert snapshot.logical_format == "pdb"
    assert snapshot.is_gzip is False


def test_gzip_and_plain_input_share_content_id_but_not_raw_hash(tmp_path: Path) -> None:
    source = b"data_demo\n_atom_site.id 1\n"
    plain = tmp_path / "same-name.cif"
    compressed = tmp_path / "same-name.cif.gz"
    plain.write_bytes(source)
    with gzip.open(compressed, "wb") as handle:
        handle.write(source)

    plain_snapshot = SourceSnapshot.from_path(plain)
    gzip_snapshot = SourceSnapshot.from_path(compressed)

    assert gzip_snapshot.decompressed_bytes == source
    assert gzip_snapshot.content_id == plain_snapshot.content_id
    assert gzip_snapshot.raw_sha256 != plain_snapshot.raw_sha256
    assert gzip_snapshot.logical_format == "mmcif"
    assert gzip_snapshot.is_gzip is True


def test_same_basename_with_different_content_does_not_collide(tmp_path: Path) -> None:
    first_dir = tmp_path / "one"
    second_dir = tmp_path / "two"
    first_dir.mkdir()
    second_dir.mkdir()
    first = first_dir / "structure.pdb"
    second = second_dir / "structure.pdb"
    first.write_bytes(b"first")
    second.write_bytes(b"second")

    first_snapshot = capture_snapshot(first)
    second_snapshot = capture_snapshot(second)

    assert first_snapshot.display_name == second_snapshot.display_name
    assert first_snapshot.content_id != second_snapshot.content_id


def test_snapshot_is_independent_from_later_path_replacement(tmp_path: Path) -> None:
    path = tmp_path / "mutable.pdb"
    original = b"original bytes"
    replacement = b"replacement bytes"
    path.write_bytes(original)

    snapshot = capture_snapshot(path)
    path.write_bytes(replacement)

    assert snapshot.decompressed_bytes == original
    assert snapshot.content_id == hashlib.sha256(original).hexdigest()
    assert snapshot.raw_sha256 == hashlib.sha256(original).hexdigest()


def test_raw_byte_limit_is_checked_before_gzip_expansion(tmp_path: Path) -> None:
    path = tmp_path / "expanded.pdb.gz"
    with gzip.open(path, "wb") as handle:
        handle.write(b"A" * 256)

    with pytest.raises(SnapshotLimitError, match="raw bytes") as error:
        capture_snapshot(
            path,
            ParseLimits(max_raw_bytes=10, max_decompressed_bytes=10_000),
        )

    assert error.value.limit_name == "raw bytes"


def test_decompressed_byte_limit_rejects_expansion_without_unbounded_read(
    tmp_path: Path,
) -> None:
    path = tmp_path / "expanded.pdb.gz"
    with gzip.open(path, "wb") as handle:
        handle.write(b"A" * 256)

    with pytest.raises(SnapshotLimitError, match="decompressed bytes"):
        capture_snapshot(
            path,
            ParseLimits(max_raw_bytes=10_000, max_decompressed_bytes=32),
        )


def test_truncated_gzip_is_rejected_with_snapshot_error(tmp_path: Path) -> None:
    path = tmp_path / "truncated.pdb.gz"
    path.write_bytes(gzip.compress(b"ATOM\n")[:-2])

    with pytest.raises(SnapshotError, match="truncated gzip"):
        capture_snapshot(path)


def test_malformed_gzip_is_rejected_with_snapshot_error(tmp_path: Path) -> None:
    path = tmp_path / "malformed.pdb.gz"
    path.write_bytes(b"not a gzip stream")

    with pytest.raises(SnapshotError, match="unable to capture"):
        capture_snapshot(path)


def test_gzip_with_trailing_garbage_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "trailing-garbage.pdb.gz"
    path.write_bytes(gzip.compress(b"ATOM\n") + b"garbage")

    with pytest.raises(SnapshotError, match="unable to capture"):
        capture_snapshot(path)


def test_gzip_with_zero_padding_is_accepted_until_eof(tmp_path: Path) -> None:
    path = tmp_path / "zero-padding.pdb.gz"
    path.write_bytes(gzip.compress(b"ATOM\n") + b"\x00" * 32)

    snapshot = capture_snapshot(path)

    assert snapshot.decompressed_bytes == b"ATOM\n"


def test_gzip_with_zero_padding_then_nonzero_garbage_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "zero-padding-garbage.pdb.gz"
    path.write_bytes(gzip.compress(b"ATOM\n") + b"\x00" * 32 + b"x")

    with pytest.raises(SnapshotError):
        capture_snapshot(path)


def test_concatenated_gzip_members_are_explicitly_concatenated(tmp_path: Path) -> None:
    path = tmp_path / "concatenated.pdb.gz"
    path.write_bytes(gzip.compress(b"first\n") + gzip.compress(b"second\n"))

    snapshot = capture_snapshot(path)

    assert snapshot.decompressed_bytes == b"first\nsecond\n"


def test_count_limits_have_clear_diagnostics() -> None:
    limits = ParseLimits(
        max_raw_bytes=10,
        max_decompressed_bytes=10,
        max_models=1,
        max_chains=2,
        max_residues=3,
        max_atoms=4,
    )

    with pytest.raises(SnapshotLimitError, match="models") as error:
        limits.validate_counts(models=2, chains=1, residues=1, atoms=1)
    assert error.value.limit_name == "models"
    assert error.value.observed == 2


def test_snapshot_rejects_unknown_logical_format(tmp_path: Path) -> None:
    path = tmp_path / "input.txt"
    path.write_bytes(b"not a structure")

    with pytest.raises(ValueError, match="unsupported logical format"):
        capture_snapshot(path)
