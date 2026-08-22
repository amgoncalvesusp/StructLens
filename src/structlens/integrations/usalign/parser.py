"""Parse the stable, human-readable sections of US-align output."""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import cast

from .executable import USAlignOutputError


@dataclass(frozen=True, slots=True)
class USAlignAlignedPair:
    """One column of US-align's gapped sequence rendering."""

    reference_index: int | None
    target_index: int | None
    reference_one_letter: str | None
    target_one_letter: str | None


@dataclass(frozen=True, slots=True)
class USAlignTransform:
    """Rigid transform emitted by US-align, when that output section exists."""

    translation: tuple[float, float, float]
    rotation: tuple[tuple[float, float, float], ...]


@dataclass(frozen=True, slots=True)
class USAlignParsedOutput:
    """Structured, non-persistent projection of a US-align stdout report."""

    tm_score: float
    aligned_pairs: tuple[USAlignAlignedPair, ...]
    transform: USAlignTransform | None
    version: str | None
    metadata: dict[str, str]


_TM_SCORE = re.compile(r"^TM-score=\s*([0-9]+(?:\.[0-9]+)?)", re.MULTILINE)
_VERSION = re.compile(r"US-align\s*\(\s*Version\s+([^\s)]+)", re.IGNORECASE)
_NAME = re.compile(r"^Name of Structure_(\d+):\s*(.+?)\s*$", re.MULTILINE)
_SEQUENCE = re.compile(r"^[A-Za-z-]+$")
_MARKER = re.compile(r"^[\s:.|]+$")
_MATRIX_ROW = re.compile(
    r"^\s*[012]\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+"
    r"([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*$"
)


def parse_usalign_output(output: str) -> USAlignParsedOutput:
    """Parse TM-score, alignment columns, transform, and provenance metadata."""

    tm_match = _TM_SCORE.search(output)
    if tm_match is None:
        raise USAlignOutputError("US-align output did not contain a TM-score.")
    aligned_pairs = _parse_alignment_pairs(output.splitlines())
    if not aligned_pairs:
        raise USAlignOutputError(
            "US-align output did not contain aligned residue pairs."
        )

    version_match = _VERSION.search(output)
    metadata = {f"structure_{number}": name for number, name in _NAME.findall(output)}
    return USAlignParsedOutput(
        tm_score=float(tm_match.group(1)),
        aligned_pairs=tuple(aligned_pairs),
        transform=_parse_transform(output.splitlines()),
        version=version_match.group(1) if version_match else None,
        metadata=metadata,
    )


def _parse_alignment_pairs(lines: list[str]) -> list[USAlignAlignedPair]:
    for index in range(len(lines) - 2):
        reference = lines[index].strip()
        marker = lines[index + 1].strip()
        target = lines[index + 2].strip()
        if (
            len(reference) == len(target)
            and len(reference) > 0
            and _SEQUENCE.fullmatch(reference)
            and _MARKER.fullmatch(marker)
        ):
            return _pairs_from_gapped_sequences(reference, target)
    return []


def _pairs_from_gapped_sequences(
    reference: str, target: str
) -> list[USAlignAlignedPair]:
    pairs: list[USAlignAlignedPair] = []
    reference_index = 0
    target_index = 0
    for reference_letter, target_letter in zip(
        reference.upper(), target.upper(), strict=True
    ):
        current_reference = None if reference_letter == "-" else reference_index
        current_target = None if target_letter == "-" else target_index
        if current_reference is not None:
            reference_index += 1
        if current_target is not None:
            target_index += 1
        if current_reference is not None or current_target is not None:
            pairs.append(
                USAlignAlignedPair(
                    current_reference,
                    current_target,
                    None if current_reference is None else reference_letter,
                    None if current_target is None else target_letter,
                )
            )
    return pairs


def _parse_transform(lines: list[str]) -> USAlignTransform | None:
    rows: list[tuple[float, float, float, float]] = []
    for line in lines:
        match = _MATRIX_ROW.match(line)
        if match is not None:
            values = tuple(float(value) for value in match.groups())
            rows.append(cast(tuple[float, float, float, float], values))
    if len(rows) < 3:
        return None
    matrix = rows[-3:]
    return USAlignTransform(
        translation=(matrix[0][0], matrix[1][0], matrix[2][0]),
        rotation=tuple((row[1], row[2], row[3]) for row in matrix),
    )


__all__ = [
    "USAlignAlignedPair",
    "USAlignParsedOutput",
    "USAlignTransform",
    "parse_usalign_output",
]
