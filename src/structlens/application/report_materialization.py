"""Ephemeral source files required by structural-analysis backends."""

from __future__ import annotations

from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory

from structlens.core.models import ProteinStructure
from structlens.core.parsing import ParsedStructure, SourceSnapshot


@contextmanager
def materialized_analysis_inputs(
    reference: ParsedStructure,
    target: ParsedStructure,
    reference_snapshot: SourceSnapshot,
    target_snapshot: SourceSnapshot,
) -> Iterator[tuple[ProteinStructure, ProteinStructure]]:
    """Expose immutable snapshots as short-lived files for structural tools.

    The parser boundary intentionally accepts snapshots without paths.  A
    structural backend such as US-align still needs a filesystem path, so this
    boundary materializes the exact decompressed bytes into a private temporary
    directory for the duration of one analysis call.  The returned structures
    are copies with ephemeral paths; the parsed structures and report inputs
    remain unchanged and no temporary path can enter report identity.
    """

    with TemporaryDirectory(prefix="structlens-analysis-") as temporary_directory:
        root = Path(temporary_directory)
        reference_path = _write_snapshot(reference_snapshot, root, "reference")
        target_path = _write_snapshot(target_snapshot, root, "target")
        yield (
            _with_source_path(reference.protein_structure, reference_path),
            _with_source_path(target.protein_structure, target_path),
        )


def _write_snapshot(snapshot: SourceSnapshot, directory: Path, role: str) -> Path:
    suffix = ".pdb" if snapshot.logical_format == "pdb" else ".cif"
    path = directory / f"{role}{suffix}"
    path.write_bytes(snapshot.decompressed_bytes)
    return path


def _with_source_path(structure: ProteinStructure, path: Path) -> ProteinStructure:
    source_path = str(path)
    chains = tuple(replace(chain, source_path=source_path) for chain in structure.chains)
    return replace(structure, chains=chains, source_path=source_path)


__all__ = ["materialized_analysis_inputs"]
