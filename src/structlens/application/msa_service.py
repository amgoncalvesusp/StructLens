"""Application service for reproducible MSA row/ResidueId mapping."""

from __future__ import annotations

import difflib
import math
import subprocess
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory

from structlens.core.msa import (
    AnalysisSequence,
    MSAColumn,
    MSAResidueCell,
    MSASettings,
    MultipleSequenceAlignment,
    SequenceResidueRef,
)

_CANONICAL = frozenset("ACDEFGHIKLMNPQRSTVWY")


def _internal_ids(count: int) -> tuple[str, ...]:
    return tuple(f"SLSEQ{index:06d}" for index in range(1, count + 1))


def _reference_label(ref: SequenceResidueRef | None, previous: str | None, insertion_index: int) -> str:
    if ref is not None and ref.residue_id is not None:
        residue = ref.residue_id
        chain = residue.chain_id or "?"
        return f"{chain}:{residue.auth_seq_id}{residue.insertion_code or ''}"
    if previous:
        return f"{previous}+{insertion_index}"
    return f"N+{insertion_index}"


def _fallback_rows(sequences: Sequence[AnalysisSequence]) -> tuple[tuple[str, str], ...]:
    """Deterministically align every row to the first row without shell calls."""

    if not sequences:
        raise ValueError("at least one sequence is required")
    reference = sequences[0].sequence
    rows: list[list[str]] = [list(reference)]
    for sequence in sequences[1:]:
        matcher = difflib.SequenceMatcher(a=reference, b=sequence.sequence, autojunk=False)
        target_row: list[str] = []
        ref_row = rows[0]
        # Rebuild the existing columns when an insertion is observed. This is
        # intentionally deterministic and is used only when MUSCLE is absent.
        for tag, i1, i2, j1, j2 in matcher.get_opcodes():
            if tag == "equal":
                target_row.extend(sequence.sequence[j1:j2])
            elif tag == "delete":
                target_row.extend("-" for _ in range(i1, i2))
            elif tag == "replace":
                width = max(i2 - i1, j2 - j1)
                target_row.extend(sequence.sequence[j1:j2])
                target_row.extend("-" for _ in range(width - (j2 - j1)))
            elif tag == "insert":
                target_row.extend(sequence.sequence[j1:j2])
        # A progressive fallback can only safely pad the reference row here;
        # columns remain explicit and never receive fabricated ResidueIds.
        if len(target_row) > len(ref_row):
            rows[0].extend("-" for _ in range(len(target_row) - len(ref_row)))
        elif len(target_row) < len(ref_row):
            target_row.extend("-" for _ in range(len(ref_row) - len(target_row)))
        rows.append(target_row)
    width = max(len(row) for row in rows)
    rows = [row + ["-"] * (width - len(row)) for row in rows]
    ids = _internal_ids(len(sequences))
    return tuple((internal_id, "".join(row)) for internal_id, row in zip(ids, rows, strict=True))


def _columns(sequences: Sequence[AnalysisSequence], rows: Sequence[tuple[str, str]]) -> tuple[MSAColumn, ...]:
    ref = sequences[0]
    width = len(rows[0][1])
    refs = iter(ref.residues)
    previous: str | None = None
    insertion_index = 0
    output: list[MSAColumn] = []
    for column_index in range(width):
        ref_char = rows[0][1][column_index]
        ref_residue = next(refs, None) if ref_char != "-" else None
        if ref_residue is not None:
            label = _reference_label(ref_residue, previous, insertion_index)
            previous = label
            insertion_index = 0
        else:
            insertion_index += 1
            label = _reference_label(None, previous, insertion_index)
        cells: list[MSAResidueCell] = []
        for sequence, (_, aligned) in zip(sequences, rows, strict=True):
            character = aligned[column_index]
            # Map a non-gap row character to its source sequence index by
            # counting non-gap symbols to the left of this alignment column.
            residue = None
            if character != "-":
                sequence_index = sum(char != "-" for char in aligned[:column_index])
                if sequence_index < len(sequence.residues):
                    residue = sequence.residues[sequence_index]
            cells.append(MSAResidueCell(sequence.structure_id, column_index, residue, character))
        non_gap = sum(cell.character != "-" for cell in cells)
        gap_fraction = 1.0 - (non_gap / len(cells))
        ambiguous = sum(cell.character.upper() in {"X", "B", "Z", "J"} for cell in cells)
        valid = [cell.character.upper() for cell in cells if cell.character.upper() in _CANONICAL]
        entropy = None
        conservation = None
        if len(valid) >= 2:
            frequencies = {symbol: valid.count(symbol) / len(valid) for symbol in set(valid)}
            entropy = -sum(value * math.log2(value) for value in frequencies.values())
            conservation = 1.0 - entropy / math.log2(20)
        output.append(
            MSAColumn(
                column_index,
                label,
                ref_residue.residue_id if ref_residue is not None else None,
                tuple(cells),
                non_gap,
                gap_fraction,
                ambiguous / len(cells),
                conservation,
                entropy,
            )
        )
    return tuple(output)


@dataclass(frozen=True, slots=True)
class MuscleAlignmentEngine:
    executable: str | Path | None = None

    def align(self, sequences: Sequence[AnalysisSequence], settings: MSASettings) -> MultipleSequenceAlignment:
        if not sequences:
            raise ValueError("at least one sequence is required")
        rows = _fallback_rows(sequences)
        algorithm = "fallback"
        provenance = ("deterministic fallback alignment",)
        if self.executable is not None:
            executable = Path(self.executable)
            if not executable.exists():
                raise FileNotFoundError(f"MUSCLE executable not found: {executable}")
            with TemporaryDirectory(prefix="structlens-msa-") as directory:
                input_path = Path(directory) / "input.fasta"
                output_path = Path(directory) / "output.afa"
                input_path.write_text(
                    "".join(f">{identifier}\n{sequence.sequence}\n" for identifier, sequence in zip(_internal_ids(len(sequences)), sequences, strict=True)),
                    encoding="utf-8",
                )
                subprocess.run(
                    [str(executable), "-align", str(input_path), "-output", str(output_path)],
                    check=True,
                    shell=False,
                    capture_output=True,
                    text=True,
                )
                rows = parse_alignment(output_path.read_text(encoding="utf-8"), sequences)
                algorithm = settings.algorithm
                provenance = (f"MUSCLE executable: {executable}",)
        return MultipleSequenceAlignment(tuple(sequences), tuple(rows), _columns(sequences, rows), sequences[0].structure_id, algorithm, provenance)


def parse_alignment(text: str, sequences: Sequence[AnalysisSequence]) -> tuple[tuple[str, str], ...]:
    """Parse FASTA alignment and reconnect generated IDs to source sequences."""

    expected = _internal_ids(len(sequences))
    parsed: dict[str, str] = {}
    current: str | None = None
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            current = line[1:].split()[0]
            if current in parsed:
                raise ValueError(f"duplicate alignment ID: {current}")
            parsed[current] = ""
        elif current is not None:
            parsed[current] += line
    if set(parsed) != set(expected):
        missing = sorted(set(expected) - set(parsed))
        unexpected = sorted(set(parsed) - set(expected))
        raise ValueError(f"alignment IDs mismatch; missing={missing}, unexpected={unexpected}")
    widths = {len(parsed[identifier]) for identifier in expected}
    if len(widths) != 1:
        raise ValueError("alignment rows must have equal length")
    return tuple((identifier, parsed[identifier]) for identifier in expected)


def align_sequences(
    sequences: Sequence[AnalysisSequence],
    settings: MSASettings | None = None,
    executable: str | Path | None = None,
) -> MultipleSequenceAlignment:
    return MuscleAlignmentEngine(executable).align(sequences, settings or MSASettings())


__all__ = ["MuscleAlignmentEngine", "align_sequences", "parse_alignment"]
