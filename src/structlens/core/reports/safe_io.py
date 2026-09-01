"""Bounded, link-resistant atomic file primitives for report artifacts."""

from __future__ import annotations

import os
import stat
import tempfile
import uuid
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

_AncestorIdentity = tuple[Path, tuple[int, int, int, int, int] | None]


class SafeReadLimitError(ValueError):
    """A bounded descriptor read observed more bytes than its configured cap."""

    def __init__(self, label: str, limit: int, observed: int) -> None:
        self.limit = limit
        self.observed = observed
        super().__init__(f"{label} exceeds the maximum supported size")


def _is_link(path: Path, info: os.stat_result) -> bool:
    reparse = getattr(stat, "FILE_ATTRIBUTE_REPARSE_POINT", 0x400)
    return stat.S_ISLNK(info.st_mode) or path.is_symlink() or bool(getattr(info, "st_file_attributes", 0) & reparse)


def validate_safe_ancestors(path: str | Path, *, label: str = "file") -> None:
    """Reject links/reparse points in every existing ancestor of ``path``."""

    _ancestor_state(path, label=label)


def _ancestor_identity(info: os.stat_result) -> tuple[int, int, int, int, int]:
    reparse = getattr(stat, "FILE_ATTRIBUTE_REPARSE_POINT", 0x400)
    return (
        info.st_dev,
        info.st_ino,
        info.st_nlink,
        info.st_mode,
        int(getattr(info, "st_file_attributes", 0) & reparse),
    )


def _ancestor_state(path: str | Path, *, label: str) -> tuple[_AncestorIdentity, ...]:
    absolute = Path(os.path.abspath(path))
    state: list[_AncestorIdentity] = []
    for ancestor in absolute.parents:
        try:
            info = ancestor.lstat()
        except FileNotFoundError:
            state.append((ancestor, None))
            continue
        if _is_link(ancestor, info):
            raise ValueError(f"{label} ancestor must not be a link")
        if not stat.S_ISDIR(info.st_mode):
            raise ValueError(f"{label} ancestor must be a directory")
        state.append((ancestor, _ancestor_identity(info)))
    return tuple(state)


def _assert_ancestor_state(state: tuple[_AncestorIdentity, ...], *, label: str) -> None:
    """Fail closed if an ancestor's identity or link state changed."""

    for ancestor, expected in state:
        try:
            info = ancestor.lstat()
        except FileNotFoundError:
            actual = None
        else:
            if _is_link(ancestor, info) or not stat.S_ISDIR(info.st_mode):
                raise ValueError(f"{label} ancestor changed or became a link")
            actual = _ancestor_identity(info)
        if actual != expected:
            raise ValueError(f"{label} ancestor changed while reading or writing")


def _regular(path: Path, *, label: str) -> os.stat_result | None:
    validate_safe_ancestors(path, label=label)
    try:
        info = path.lstat()
    except FileNotFoundError:
        return None
    if _is_link(path, info) or not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
        raise ValueError(f"{label} must be a regular, non-link file")
    return info


def _read_identity(info: os.stat_result) -> tuple[int, int, int, int, int, int, int]:
    return (
        info.st_dev,
        info.st_ino,
        info.st_mode,
        info.st_nlink,
        info.st_size,
        info.st_mtime_ns,
        info.st_ctime_ns,
    )


def read_bounded_bytes(path: str | Path, *, max_bytes: int, label: str = "file") -> bytes:
    source = Path(path)
    ancestor_state = _ancestor_state(source, label=label)
    info = _regular(source, label=label)
    if info is None:
        raise FileNotFoundError(source)
    if info.st_size > max_bytes:
        raise SafeReadLimitError(label, max_bytes, info.st_size)
    flags = os.O_RDONLY | getattr(os, "O_BINARY", 0) | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(source, flags)
    except OSError as exc:
        raise ValueError(f"{label} cannot be opened safely") from exc
    try:
        opened = os.fstat(descriptor)
        current = source.lstat()
        if (
            _is_link(source, current)
            or not stat.S_ISREG(opened.st_mode)
            or opened.st_nlink != 1
            or _read_identity(opened) != _read_identity(current)
            or _read_identity(opened) != _read_identity(info)
        ):
            raise ValueError(f"{label} changed or is a link")
        _assert_ancestor_state(ancestor_state, label=label)
        chunks: list[bytes] = []
        total = 0
        while total <= max_bytes:
            chunk = os.read(descriptor, min(1024 * 1024, max_bytes + 1 - total))
            if not chunk:
                break
            chunks.append(chunk)
            total += len(chunk)
        if total > max_bytes:
            raise SafeReadLimitError(label, max_bytes, total)
        finished = os.fstat(descriptor)
        current = source.lstat()
        if (
            not stat.S_ISREG(finished.st_mode)
            or finished.st_nlink != 1
            or _is_link(source, current)
            or _read_identity(finished) != _read_identity(opened)
            or _read_identity(current) != _read_identity(opened)
        ):
            raise ValueError(f"{label} changed while reading")
        _assert_ancestor_state(ancestor_state, label=label)
        return b"".join(chunks)
    finally:
        os.close(descriptor)


@dataclass(frozen=True, slots=True)
class _StableDirectory:
    path: Path
    descriptor: int | None = None
    windows_handles: tuple[int, ...] = ()


def _lock_windows_directories(parent: Path, *, label: str) -> tuple[int, ...]:
    import ctypes
    from ctypes import wintypes

    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    create_file = kernel32.CreateFileW
    create_file.argtypes = (
        wintypes.LPCWSTR,
        wintypes.DWORD,
        wintypes.DWORD,
        wintypes.LPVOID,
        wintypes.DWORD,
        wintypes.DWORD,
        wintypes.HANDLE,
    )
    create_file.restype = wintypes.HANDLE
    close_handle = kernel32.CloseHandle
    close_handle.argtypes = (wintypes.HANDLE,)
    close_handle.restype = wintypes.BOOL
    invalid = ctypes.c_void_p(-1).value
    handles: list[int] = []
    paths = tuple(reversed((parent, *parent.parents)))
    try:
        for directory in paths:
            handle = create_file(
                str(directory),
                0x80000000,
                0x00000001 | 0x00000002,
                None,
                3,
                0x02000000 | 0x00200000,
                None,
            )
            value = int(handle) if handle is not None else 0
            if value in (0, invalid):
                raise OSError(ctypes.get_last_error(), f"cannot lock {directory}")
            handles.append(value)
    except OSError as exc:
        for handle in reversed(handles):
            close_handle(wintypes.HANDLE(handle))
        raise ValueError(f"{label} directory cannot be held safely") from exc
    return tuple(handles)


def _unlock_windows_directories(handles: tuple[int, ...]) -> None:
    import ctypes
    from ctypes import wintypes

    close_handle = ctypes.WinDLL("kernel32", use_last_error=True).CloseHandle
    close_handle.argtypes = (wintypes.HANDLE,)
    close_handle.restype = wintypes.BOOL
    for handle in reversed(handles):
        close_handle(wintypes.HANDLE(handle))


@contextmanager
def _stable_posix_directory(parent: Path, *, label: str) -> Iterator[_StableDirectory]:  # pragma: no cover
    """Hold one POSIX directory descriptor for descriptor-relative writes."""

    flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(parent, flags)
    except OSError as exc:
        raise ValueError(f"{label} directory cannot be held safely") from exc
    try:
        opened = os.fstat(descriptor)
        current = parent.lstat()
        if not stat.S_ISDIR(opened.st_mode) or _ancestor_identity(opened) != _ancestor_identity(current):
            raise ValueError(f"{label} directory changed while opening")
        yield _StableDirectory(parent, descriptor=descriptor)
    finally:
        os.close(descriptor)


@contextmanager
def _stable_directory(parent: Path, *, label: str) -> Iterator[_StableDirectory]:
    """Hold the resolved parent for every create, fsync and replace operation."""

    if os.name != "nt":  # pragma: no cover - exercised on POSIX hosts
        with _stable_posix_directory(parent, label=label) as directory:
            yield directory
        return
    handles = _lock_windows_directories(parent, label=label)
    try:
        yield _StableDirectory(parent, windows_handles=handles)
    finally:
        _unlock_windows_directories(handles)


def _temporary_file(directory: _StableDirectory, target: Path) -> tuple[int, str]:
    if directory.descriptor is None:
        descriptor, temporary = tempfile.mkstemp(prefix=f".{target.name}.", dir=str(directory.path))
        return descriptor, temporary
    for _attempt in range(10):  # pragma: no cover - POSIX descriptor path
        name = f".{target.name}.{uuid.uuid4().hex}"
        flags = os.O_RDWR | os.O_CREAT | os.O_EXCL | getattr(os, "O_BINARY", 0)
        try:
            return os.open(name, flags, 0o600, dir_fd=directory.descriptor), name
        except FileExistsError:
            continue
    raise FileExistsError(f"cannot allocate a private temporary file for {target.name}")  # pragma: no cover


def _replace_temporary(directory: _StableDirectory, temporary: str, target: Path) -> None:
    if directory.descriptor is None:
        os.replace(temporary, target)
        return
    os.replace(  # pragma: no cover - POSIX descriptor path
        temporary, target.name, src_dir_fd=directory.descriptor, dst_dir_fd=directory.descriptor
    )
    os.fsync(directory.descriptor)  # pragma: no cover


def _unlink_temporary(directory: _StableDirectory, temporary: str) -> None:
    try:
        if directory.descriptor is None:
            os.unlink(temporary)
        else:  # pragma: no cover - POSIX descriptor path
            os.unlink(temporary, dir_fd=directory.descriptor)
    except FileNotFoundError:
        pass


def atomic_write_bytes(
    path: str | Path,
    data: bytes,
    *,
    max_bytes: int,
    label: str = "file",
    idempotent: bool = False,
    replace_existing: bool = False,
) -> None:
    if not isinstance(data, bytes):
        raise TypeError("atomic write data must be bytes")
    if len(data) > max_bytes:
        raise ValueError(f"{label} exceeds the maximum supported size")
    target = Path(os.path.abspath(path))
    parent = target.parent
    validate_safe_ancestors(parent, label=label)
    parent.mkdir(parents=True, exist_ok=True)
    validate_safe_ancestors(parent, label=label)
    ancestor_state = _ancestor_state(target, label=label)
    with _stable_directory(parent, label=label) as directory:
        _assert_ancestor_state(ancestor_state, label=label)
        existing = _regular(target, label=label)
        existed_before_write = existing is not None
        if existing is not None:
            if idempotent and read_bounded_bytes(target, max_bytes=max_bytes, label=label) == data:
                return
            if idempotent or not replace_existing:
                raise FileExistsError(f"conflicting {label} changed or already exists")
        descriptor, temporary = _temporary_file(directory, target)
        try:
            with os.fdopen(descriptor, "wb") as handle:
                handle.write(data)
                handle.flush()
                os.fsync(handle.fileno())
            _assert_ancestor_state(ancestor_state, label=label)
            current = _regular(target, label=label)
            if current is not None and not existed_before_write:
                raise FileExistsError(f"{label} appeared while writing")
            if current is None and existed_before_write:
                raise FileExistsError(f"{label} disappeared while writing")
            if current is not None and existing is not None and _read_identity(current) != _read_identity(existing):
                raise FileExistsError(f"{label} changed while writing")
            _assert_ancestor_state(ancestor_state, label=label)
            _replace_temporary(directory, temporary, target)
            temporary = ""
        finally:
            if temporary:
                _unlink_temporary(directory, temporary)


def atomic_write_text(
    path: str | Path,
    text: str,
    *,
    max_bytes: int,
    label: str = "file",
    idempotent: bool = False,
    replace_existing: bool = False,
) -> None:
    if not isinstance(text, str):
        raise TypeError("atomic write text must be a string")
    atomic_write_bytes(
        path,
        text.encode("utf-8"),
        max_bytes=max_bytes,
        label=label,
        idempotent=idempotent,
        replace_existing=replace_existing,
    )


__all__ = [
    "SafeReadLimitError",
    "atomic_write_bytes",
    "atomic_write_text",
    "read_bounded_bytes",
    "validate_safe_ancestors",
]
