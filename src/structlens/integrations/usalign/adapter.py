"""Structural-alignment engine backed by a locally installed US-align binary."""

from __future__ import annotations

import subprocess
import tempfile
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

from structlens.core.models import (
    CorrespondenceStatus,
    ProteinChain,
    ResidueCorrespondence,
    StructuralAlignmentSettings,
)

from .executable import USAlignExecutionError, USAlignOutputError, resolve_backend
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
    ) -> None:
        self._configured_executable = executable

    def align(
        self,
        reference: ProteinChain,
        target: ProteinChain,
        settings: StructuralAlignmentSettings,
    ) -> USAlignAlignmentResult:
        """Run alignment for two chains whose source files are registered locally."""

        configured = self._configured_executable
        if configured is None and settings.executable not in {"", "USalign", "US-align"}:
            configured = settings.executable
        backend = resolve_backend(configured)
        with tempfile.TemporaryDirectory(prefix="structlens-usalign-") as temporary_directory:
            directory = Path(temporary_directory)
            reference_path = _write_selected_trace(reference, directory, "reference")
            target_path = _write_selected_trace(target, directory, "target")
            try:
                completed = subprocess.run(
                    [str(backend.path), str(reference_path), str(target_path)],
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
        result = self._result_from_parsed(reference, target, parsed)
        return USAlignAlignmentResult(
            correspondences=result.correspondences,
            tm_score=result.tm_score,
            transform=result.transform,
            executable_version=result.executable_version or backend.version,
            metadata={
                **dict(result.metadata),
                "backend": "US-align",
                "binary_source": backend.source,
                "platform": backend.platform,
                "command_options": "selected-model-chain C-alpha traces; reference then target",
            },
        )

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


def _write_selected_trace(chain: ProteinChain, directory: Path, role: str) -> Path:
    path = directory / f"{role}.pdb"
    path.write_bytes(_selected_ca_pdb(chain))
    return path


def _selected_ca_pdb(chain: ProteinChain) -> bytes:
    """Serialize only one validated chain's ordered C-alpha trace for US-align."""

    if len(chain.residue_records) > 9_999:
        raise USAlignExecutionError("Selected chain exceeds the temporary PDB residue limit.")
    lines: list[str] = []
    for index, residue in enumerate(chain.residue_records, start=1):
        alpha_carbons = tuple(atom for atom in residue.atoms if atom.name.upper() == "CA")
        if len(alpha_carbons) != 1:
            raise USAlignExecutionError(
                f"Selected residue {residue.residue_id.auth_seq_id!r} must contain exactly one C-alpha atom."
            )
        x, y, z = alpha_carbons[0].coordinate
        if any(value < -999.999 or value > 9_999.999 for value in (x, y, z)):
            raise USAlignExecutionError("Selected C-alpha coordinate exceeds the temporary PDB field range.")
        lines.append(
            f"ATOM  {index:5d}  CA  {residue.residue_name:>3.3s} A{index:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}           C  \n"
        )
    if not lines:
        raise USAlignExecutionError("Selected chain has no residues for structural alignment.")
    lines.extend(("TER\n", "END\n"))
    return "".join(lines).encode("ascii")
