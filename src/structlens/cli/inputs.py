"""Boundary-safe structure capture and selection helpers for the CLI."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from structlens.core.parsing import (
    AltlocPolicy,
    AssemblyScope,
    InputSelection,
    ParsedStructure,
    SourceSnapshot,
    StructureFormat,
    capture_snapshot,
    load_structure_evidence,
)


@dataclass(frozen=True, slots=True)
class CapturedInput:
    """One bounded source together with its validated selection and parse."""

    snapshot: SourceSnapshot
    selection: InputSelection
    parsed: ParsedStructure


def capture_input(
    path: str | Path,
    *,
    model_id: str = "1",
    author_chain_ids: Sequence[str] = (),
    label_chain_ids: Sequence[str] = (),
    altloc_policy: str | AltlocPolicy = AltlocPolicy.HIGHEST_OCCUPANCY,
    assembly_scope: str | AssemblyScope = AssemblyScope.ASYMMETRIC_UNIT,
) -> CapturedInput:
    """Capture and validate one supported structure before any science runs."""

    source_path = _path(path)
    model = _nonempty(model_id, "model")
    authors = _chain_values(author_chain_ids, "author chain")
    labels = _chain_values(label_chain_ids, "label chain")
    snapshot = capture_snapshot(source_path)
    # Biopython assigns implicit PDB models the serial ``0``.  Treat the
    # documented CLI default ``1`` as the first implicit PDB model so legacy
    # PDB invocations remain usable while explicit model selection stays strict.
    if model == "1" and snapshot.logical_format == "pdb" and not _has_explicit_model(snapshot):
        model = "0"
    selection = InputSelection(
        content_id=snapshot.content_id,
        display_name=snapshot.display_name,
        format=StructureFormat(snapshot.logical_format),
        model_id=model,
        author_chain_ids=authors,
        label_chain_ids=labels,
        altloc_policy=altloc_policy,
        assembly_scope=assembly_scope,
        path=str(source_path),
    )
    parsed = load_structure_evidence(snapshot, selection=selection)
    return CapturedInput(snapshot, parsed.selection, parsed)


def _path(value: str | Path) -> Path:
    candidate = Path(value)
    if not str(candidate).strip():
        raise ValueError("structure path must not be empty")
    if not candidate.exists():
        raise ValueError(f"structure path does not exist: {candidate}")
    if not candidate.is_file():
        raise ValueError(f"structure path is not a regular file: {candidate}")
    return candidate


def _nonempty(value: object, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must not be empty")
    return value.strip()


def _chain_values(values: Sequence[str], label: str) -> tuple[str, ...]:
    normalized = tuple(_nonempty(value, label) for value in values)
    if len(set(normalized)) != len(normalized):
        raise ValueError(f"{label} values must not repeat")
    return normalized


def _has_explicit_model(snapshot: SourceSnapshot) -> bool:
    """Distinguish PDB's implicit model 0 from explicit MODEL serials."""

    return any(line.startswith("MODEL") for line in snapshot.decompressed_bytes.decode("utf-8").splitlines())


__all__ = ["CapturedInput", "capture_input"]
