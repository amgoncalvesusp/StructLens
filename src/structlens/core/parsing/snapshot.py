"""Bounded, content-addressed source capture for structure parsers."""

from __future__ import annotations

import hashlib
import zlib
from dataclasses import dataclass
from pathlib import Path

from .limits import ParseLimits, SnapshotError, SnapshotLimitError

_READ_CHUNK_SIZE = 64 * 1024
_GZIP_WBITS = 16 + zlib.MAX_WBITS
_SHA256_HEX_LENGTH = hashlib.sha256().digest_size * 2


@dataclass(frozen=True, slots=True)
class SourceSnapshot:
    """Immutable bytes and hashes captured before any structure parsing.

    ``content_id`` is the SHA-256 digest of ``decompressed_bytes``.  The raw
    digest deliberately remains separate because gzip containers can differ
    while carrying the same logical structure source.
    """

    decompressed_bytes: bytes
    raw_sha256: str
    content_id: str
    decompressed_sha256: str
    display_name: str
    logical_format: str
    is_gzip: bool
    raw_bytes: bytes | None = None

    def __post_init__(self) -> None:
        content = bytes(self.decompressed_bytes)
        object.__setattr__(self, "decompressed_bytes", content)
        if not self.display_name:
            raise ValueError("display_name must not be empty")
        if self.logical_format not in {"pdb", "mmcif"}:
            raise ValueError("logical_format must be 'pdb' or 'mmcif'")
        if not isinstance(self.is_gzip, bool):
            raise TypeError("is_gzip must be a boolean")
        _validate_sha256("raw_sha256", self.raw_sha256)
        _validate_sha256("content_id", self.content_id)
        _validate_sha256("decompressed_sha256", self.decompressed_sha256)
        expected_content_hash = hashlib.sha256(content).hexdigest()
        if self.content_id != expected_content_hash:
            raise ValueError("content_id does not match decompressed_bytes")
        if self.decompressed_sha256 != expected_content_hash:
            raise ValueError("decompressed_sha256 does not match decompressed_bytes")
        raw = self.raw_bytes
        if raw is not None:
            raw = bytes(raw)
            if hashlib.sha256(raw).hexdigest() != self.raw_sha256:
                raise ValueError("raw_sha256 does not match raw_bytes")
            object.__setattr__(self, "raw_bytes", raw)
        elif self.raw_sha256 == expected_content_hash:
            object.__setattr__(self, "raw_bytes", content)

    @classmethod
    def from_path(
        cls,
        path: str | Path,
        limits: ParseLimits | None = None,
    ) -> SourceSnapshot:
        """Capture ``path`` once using bounded reads."""

        return capture_snapshot(path, limits)

    @property
    def raw_hash(self) -> str:
        """Compatibility spelling for the raw SHA-256 digest."""

        return self.raw_sha256

    @property
    def content_sha256(self) -> str:
        """SHA-256 digest of the logical, decompressed content."""

        return self.content_id

    @property
    def decompressed_hash(self) -> str:
        """Compatibility spelling for the decompressed SHA-256 digest."""

        return self.decompressed_sha256

    @property
    def content(self) -> bytes:
        """Compatibility view of the immutable logical source bytes."""

        return self.decompressed_bytes

    @property
    def format(self) -> str:
        """Logical structure format (``pdb`` or ``mmcif``)."""

        return self.logical_format


def capture_snapshot(
    path: str | Path,
    limits: ParseLimits | None = None,
) -> SourceSnapshot:
    """Read a supported PDB/mmCIF source in bounded chunks.

    For gzip input, the compressed byte cap is checked as bytes are read and
    the decompressed cap is checked as output is produced.  No parser receives
    a path or an unbounded file-like object from this function.
    """

    from structlens.core.reports.safe_io import SafeReadLimitError, read_bounded_bytes

    source_path = Path(path)
    logical_format, is_gzip = _logical_format(source_path)
    active_limits = limits if limits is not None else ParseLimits()
    raw_hash = hashlib.sha256()
    content = bytearray()

    try:
        raw = read_bounded_bytes(source_path, max_bytes=active_limits.max_raw_bytes, label="source")
        raw_hash.update(raw)
        decoder = _GzipDecoder(active_limits.max_decompressed_bytes) if is_gzip else None
        for start in range(0, len(raw), _READ_CHUNK_SIZE):
            chunk = raw[start : start + _READ_CHUNK_SIZE]
            if decoder is None:
                _append_decompressed(content, chunk, active_limits)
            else:
                decoder.feed(chunk, content)
        if decoder is not None:
            decoder.finish(content)
    except SafeReadLimitError as error:
        raise SnapshotLimitError("raw bytes", error.limit, error.observed) from error
    except SnapshotError:
        raise
    except (OSError, zlib.error) as error:
        raise SnapshotError(f"unable to capture {source_path.name}: {error}") from error

    decompressed = bytes(content)
    content_hash = hashlib.sha256(decompressed).hexdigest()
    return SourceSnapshot(
        decompressed_bytes=decompressed,
        raw_sha256=raw_hash.hexdigest(),
        content_id=content_hash,
        decompressed_sha256=content_hash,
        display_name=source_path.name,
        logical_format=logical_format,
        is_gzip=is_gzip,
        raw_bytes=raw,
    )


class _GzipDecoder:
    """Incremental gzip decoder that never asks zlib for unlimited output."""

    def __init__(self, max_decompressed_bytes: int) -> None:
        self._max_decompressed_bytes = max_decompressed_bytes
        self._decoder = zlib.decompressobj(_GZIP_WBITS)

    def feed(self, data: bytes, output: bytearray) -> None:
        pending = data
        while pending:
            if self._decoder.eof:
                pending = pending.lstrip(b"\x00")
                if not pending:
                    return
                self._decoder = zlib.decompressobj(_GZIP_WBITS)
            decoded = self._decoder.decompress(
                pending,
                self._remaining_output(output) + 1,
            )
            _append_decompressed(output, decoded, self._max_decompressed_bytes)
            if self._decoder.unused_data:
                pending = self._decoder.unused_data
                continue
            pending = self._decoder.unconsumed_tail
            if pending:
                continue
            self._drain(output)

    def finish(self, output: bytearray) -> None:
        if not self._decoder.eof:
            decoded = self._decoder.flush(self._remaining_output(output) + 1)
            _append_decompressed(output, decoded, self._max_decompressed_bytes)
            raise SnapshotError("invalid or truncated gzip source")

    def _drain(self, output: bytearray) -> None:
        while True:
            decoded = self._decoder.decompress(
                b"",
                self._remaining_output(output) + 1,
            )
            _append_decompressed(output, decoded, self._max_decompressed_bytes)
            if not decoded:
                return

    def _remaining_output(self, output: bytearray) -> int:
        return self._max_decompressed_bytes - len(output)


def _append_decompressed(
    output: bytearray,
    decoded: bytes,
    limits: ParseLimits | int,
) -> None:
    max_size = limits.max_decompressed_bytes if isinstance(limits, ParseLimits) else limits
    next_size = len(output) + len(decoded)
    if next_size > max_size:
        raise SnapshotLimitError("decompressed bytes", max_size, next_size)
    output.extend(decoded)


def _logical_format(path: Path) -> tuple[str, bool]:
    name = path.name.lower()
    is_gzip = name.endswith(".gz")
    if is_gzip:
        name = name[:-3]
    if name.endswith((".pdb", ".ent")):
        return "pdb", is_gzip
    if name.endswith((".cif", ".mmcif")):
        return "mmcif", is_gzip
    raise ValueError(f"unsupported logical format for source: {path.name}")


def _validate_sha256(name: str, value: str) -> None:
    if (
        not isinstance(value, str)
        or len(value) != _SHA256_HEX_LENGTH
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise ValueError(f"{name} must be a lowercase SHA-256 hexadecimal digest")


__all__ = [
    "SnapshotError",
    "SnapshotLimitError",
    "SourceSnapshot",
    "capture_snapshot",
]
