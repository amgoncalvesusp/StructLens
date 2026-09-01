from __future__ import annotations

import os
from pathlib import Path

import pytest

import structlens.core.reports.safe_io as safe_io
from structlens.core.reports.safe_io import atomic_write_bytes, read_bounded_bytes, validate_safe_ancestors


def test_ancestor_validation_handles_missing_non_directory_and_linked_ancestors(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    validate_safe_ancestors(tmp_path / "missing" / "nested" / "artifact.bin")

    non_directory = tmp_path / "not-a-directory"
    non_directory.write_bytes(b"file")
    with pytest.raises(ValueError, match="directory"):
        validate_safe_ancestors(non_directory / "artifact.bin")

    linked = tmp_path / "linked"
    linked.mkdir()
    original_is_link = safe_io._is_link

    def pretend_link(path: Path, info: os.stat_result) -> bool:
        return path == linked or original_is_link(path, info)

    monkeypatch.setattr(safe_io, "_is_link", pretend_link)
    with pytest.raises(ValueError, match="link"):
        validate_safe_ancestors(linked / "artifact.bin")


def test_ancestor_state_accepts_missing_and_rejects_link_or_identity_change(tmp_path: Path) -> None:
    missing = tmp_path / "missing"
    safe_io._assert_ancestor_state(((missing, None),), label="artifact")

    changed = tmp_path / "changed"
    changed.mkdir()
    state = safe_io._ancestor_state(changed / "artifact.bin", label="artifact")
    changed.rmdir()
    changed.write_bytes(b"not a directory")
    with pytest.raises(ValueError, match="changed|link"):
        safe_io._assert_ancestor_state(state, label="artifact")

    changed.unlink()
    changed.mkdir()
    with pytest.raises(ValueError, match="changed"):
        safe_io._assert_ancestor_state(state, label="artifact")


def test_bounded_read_rejects_open_failure_and_source_change(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    source = tmp_path / "source.bin"
    source.write_bytes(b"payload")

    def fail_open(*_args: object, **_kwargs: object) -> int:
        raise OSError("injected open failure")

    monkeypatch.setattr(safe_io.os, "open", fail_open)
    with pytest.raises(ValueError, match="opened safely"):
        read_bounded_bytes(source, max_bytes=100)

    monkeypatch.undo()
    original_is_link = safe_io._is_link
    original_open = safe_io.os.open

    def open_then_report_link(*args: object, **kwargs: object) -> int:
        descriptor = original_open(*args, **kwargs)
        monkeypatch.setattr(
            safe_io,
            "_is_link",
            lambda path, info: path == source or original_is_link(path, info),
        )
        return descriptor

    monkeypatch.setattr(safe_io.os, "open", open_then_report_link)
    with pytest.raises(ValueError, match="changed or is a link"):
        read_bounded_bytes(source, max_bytes=100)


def test_atomic_write_detects_target_disappearance_and_change(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    original_regular = safe_io._regular

    disappeared = tmp_path / "disappeared.bin"
    disappeared.write_bytes(b"old")
    calls = 0

    def disappear_on_replace(path: Path, *, label: str) -> os.stat_result | None:
        nonlocal calls
        calls += 1
        if calls == 2:
            return None
        return original_regular(path, label=label)

    monkeypatch.setattr(safe_io, "_regular", disappear_on_replace)
    with pytest.raises(FileExistsError, match="disappeared"):
        atomic_write_bytes(disappeared, b"new", max_bytes=100, replace_existing=True)

    monkeypatch.setattr(safe_io, "_regular", original_regular)
    changed = tmp_path / "changed.bin"
    changed.write_bytes(b"old")
    calls = 0

    def change_on_replace(path: Path, *, label: str) -> os.stat_result | None:
        nonlocal calls
        calls += 1
        if calls == 2:
            path.write_bytes(b"changed")
        return original_regular(path, label=label)

    monkeypatch.setattr(safe_io, "_regular", change_on_replace)
    with pytest.raises(FileExistsError, match="changed while writing"):
        atomic_write_bytes(changed, b"new", max_bytes=100, replace_existing=True)
