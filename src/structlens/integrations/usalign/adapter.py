"""Structural-alignment engine backed by a locally installed US-align binary."""

from __future__ import annotations

import subprocess
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

from structlens.core.models import (
    CorrespondenceStatus,
    ProteinChain,
    ResidueCorrespondence,
    StructuralAlignmentSettings,
)

from .executable import USAlignExecutionError, USAlignOutputError, discover_executable
from .parser import USAlignParsedOutput, USAlignTransform, parse_usalign_output


@dataclass(frozen=True, slots=True)
class USAlignAlignmentResult:
    """Adapter-local result ready to be promoted into an analysis result."""

    correspondences: tuple[ResidueCorrespondence, ...]
    tm_score: float
    transform: USAlignTransform | None
    executable_version: str | None
    metadata: Mapping[str, str]


class USAlignAdapter:
    """Run US-align without shell interpolation and map output to domain records."""

    def __init__(
        self,
        *,
        executable: str | Path | None = None,
        structure_paths: Mapping[str, str | Path],
    ) -> None:
        self._configured_executable = executable
        self._structure_paths = {
            structure_id: Path(path) for structure_id, path in structure_paths.items()
        }

    def align(
        self,
        reference: ProteinChain,
        target: ProteinChain,
        settings: StructuralAlignmentSettings,
    ) -> USAlignAlignmentResult:
        """Run alignment for two chains whose source files are registered locally."""

        reference_path = self._source_path(reference)
        target_path = self._source_path(target)
        executable = discover_executable(
            self._configured_executable
            if self._configured_executable is not None
            else settings.executable
        )
        try:
            completed = subprocess.run(
                [str(executable), str(reference_path), str(target_path)],
                capture_output=True,
                check=False,
                shell=False,
                text=True,
                timeout=settings.timeout_seconds,
            )
        except subprocess.TimeoutExpired as error:
            raise USAlignExecutionError(
                f"US-align timed out after {settings.timeout_seconds} seconds."
            ) from error
        except OSError as error:
            raise USAlignExecutionError(f"Could not run US-align: {error}") from error
        if completed.returncode != 0:
            details = (
                completed.stderr.strip() or completed.stdout.strip() or "no output"
            )
            raise USAlignExecutionError(
                f"US-align exited with code {completed.returncode}: {details}"
            )
        parsed = parse_usalign_output(completed.stdout)
        return self._result_from_parsed(reference, target, parsed)

    def _source_path(self, chain: ProteinChain) -> Path:
        try:
            return self._structure_paths[chain.structure_id]
        except KeyError as error:
            raise USAlignExecutionError(
                f"No source path is registered for structure '{chain.structure_id}'."
            ) from error

    @staticmethod
    def _result_from_parsed(
        reference: ProteinChain,
        target: ProteinChain,
        parsed: USAlignParsedOutput,
    ) -> USAlignAlignmentResult:
        correspondences: list[ResidueCorrespondence] = []
        for alignment_index, pair in enumerate(parsed.aligned_pairs):
            try:
                reference_residue = (
                    None
                    if pair.reference_index is None
                    else reference.residues[pair.reference_index]
                )
                target_residue = (
                    None
                    if pair.target_index is None
                    else target.residues[pair.target_index]
                )
            except IndexError as error:
                raise USAlignOutputError(
                    "US-align alignment columns exceed the supplied chain "
                    "residue count."
                ) from error
            correspondences.append(
                ResidueCorrespondence(
                    alignment_index=alignment_index,
                    reference=reference_residue,
                    target=target_residue,
                    reference_one_letter=pair.reference_one_letter,
                    target_one_letter=pair.target_one_letter,
                    status=_status_for(
                        pair.reference_one_letter, pair.target_one_letter
                    ),
                    mapping_source="US-align",
                )
            )
        return USAlignAlignmentResult(
            correspondences=tuple(correspondences),
            tm_score=parsed.tm_score,
            transform=parsed.transform,
            executable_version=parsed.version,
            metadata=parsed.metadata,
        )


def _status_for(
    reference_one_letter: str | None, target_one_letter: str | None
) -> CorrespondenceStatus:
    if reference_one_letter is None:
        return CorrespondenceStatus.INSERTION
    if target_one_letter is None:
        return CorrespondenceStatus.DELETION
    if reference_one_letter == target_one_letter:
        return CorrespondenceStatus.CONSERVED
    return CorrespondenceStatus.SUBSTITUTION


__all__ = ["USAlignAdapter", "USAlignAlignmentResult"]
