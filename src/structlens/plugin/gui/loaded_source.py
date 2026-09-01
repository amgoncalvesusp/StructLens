"""Immutable capture/parse boundary for one GUI structure source."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from structlens.core.models import ProteinStructure
from structlens.core.parsing import (
    InputSelection,
    SourceSnapshot,
    StructureFormat,
    capture_snapshot,
    load_structure_legacy_snapshot,
)


@dataclass(frozen=True, slots=True)
class LoadedSource:
    """One path capture and its structure parsed from exactly those bytes."""

    path: Path
    snapshot: SourceSnapshot
    structure: ProteinStructure
    pymol_object_name: str | None = None

    @classmethod
    def capture(cls, path: str | Path, *, pymol_object_name: str | None = None) -> LoadedSource:
        source_path = Path(path)
        snapshot = capture_snapshot(source_path)
        structure = load_structure_legacy_snapshot(snapshot, path=source_path)
        return cls(source_path, snapshot, structure, pymol_object_name)

    @classmethod
    def from_snapshot(
        cls,
        path: str | Path,
        snapshot: SourceSnapshot,
        *,
        pymol_object_name: str | None = None,
    ) -> LoadedSource:
        source_path = Path(path)
        structure = load_structure_legacy_snapshot(snapshot, path=source_path)
        return cls(source_path, snapshot, structure, pymol_object_name)

    def selection(self, *, model_id: str, chain_id: str) -> InputSelection:
        return InputSelection(
            self.snapshot.content_id,
            self.snapshot.display_name,
            StructureFormat(self.snapshot.logical_format),
            model_id,
            author_chain_ids=(chain_id,),
            path=str(self.path),
        )


__all__ = ["LoadedSource"]
