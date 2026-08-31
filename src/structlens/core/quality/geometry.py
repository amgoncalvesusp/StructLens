"""Conservative, coordinate-only geometry screens for normalized chains.

These checks are sanity screens for coordinate quality.  They are deliberately
not crystallographic validation: a reported distance only says that the
available coordinates fall inside (or outside) an explicit screening interval.
Pairs are considered only when source numbering provides evidence that the two
residues are adjacent in the same structure, model, and chain.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import TypeAlias

from structlens.core.evidence import Diagnostic, DiagnosticSeverity
from structlens.core.models import AtomRecord, ProteinChain, ResidueRecord

# These deliberately permissive screening limits flag coordinate discontinuity,
# not refinement-geometry outliers. Exact peptide geometry needs a versioned
# restraint library and is outside this coordinate-only gate.
DEFAULT_MIN_CN_DISTANCE_ANGSTROM = 0.0
DEFAULT_MAX_CN_DISTANCE_ANGSTROM = 2.0
DEFAULT_MIN_CA_DISTANCE_ANGSTROM = 2.5
DEFAULT_MAX_CA_DISTANCE_ANGSTROM = 4.5

ChainInput: TypeAlias = ProteinChain | Sequence[ProteinChain]


def _chains(value: ChainInput) -> tuple[ProteinChain, ...]:
    """Normalize one chain or an explicit collection without reordering it."""

    if isinstance(value, ProteinChain):
        return (value,)
    try:
        chains = tuple(value)
    except TypeError as exc:
        raise TypeError("chains must be a ProteinChain or a sequence of ProteinChain values") from exc
    if any(not isinstance(chain, ProteinChain) for chain in chains):
        raise TypeError("chains must contain ProteinChain values")
    return chains


def _validate_limits(minimum: float, maximum: float, *, name: str) -> tuple[float, float]:
    try:
        lower = float(minimum)
        upper = float(maximum)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} limits must be finite and non-negative") from exc
    if not math.isfinite(lower) or not math.isfinite(upper) or lower < 0.0 or upper < 0.0:
        raise ValueError(f"{name} limits must be finite and non-negative")
    if lower > upper:
        raise ValueError(f"{name} minimum must be less than or equal to its maximum")
    return lower, upper


def _number(value: str) -> int | None:
    """Return an integer source number, preserving unknown numbering as absent."""

    try:
        return int(value.strip())
    except (AttributeError, TypeError, ValueError):
        return None


def _insertion_after(first: str | None, second: str | None) -> bool:
    first_code = "" if first is None else first.strip().upper()
    second_code = "" if second is None else second.strip().upper()
    if not second_code:
        return False
    if not first_code:
        return True
    return second_code > first_code


def _source_adjacent(first: ResidueRecord, second: ResidueRecord) -> bool:
    """Check adjacency supported by author numbering and insertion codes."""

    first_number = _number(first.residue_id.auth_seq_id)
    second_number = _number(second.residue_id.auth_seq_id)
    if first_number is None or second_number is None:
        return False
    if second_number == first_number + 1:
        return True
    return second_number == first_number and _insertion_after(
        first.residue_id.insertion_code,
        second.residue_id.insertion_code,
    )


def _compatible_pair(chain: ProteinChain, first: ResidueRecord, second: ResidueRecord) -> bool:
    first_id = first.residue_id
    second_id = second.residue_id
    return (
        first_id.structure_id == chain.structure_id
        and second_id.structure_id == chain.structure_id
        and first_id.model_id == chain.model_id
        and second_id.model_id == chain.model_id
        and first_id.chain_id == chain.chain_id
        and second_id.chain_id == chain.chain_id
    )


def _pairs(chain: ProteinChain) -> tuple[tuple[ResidueRecord, ResidueRecord], ...]:
    records = chain.residue_records
    return tuple(
        (first, second)
        for first, second in zip(records, records[1:], strict=False)
        if _compatible_pair(chain, first, second) and _source_adjacent(first, second)
    )


def _atom(record: ResidueRecord, name: str) -> AtomRecord | None:
    target = name.upper()
    for atom in record.atoms:
        if atom.name.strip().upper() == target:
            return atom
    return None


def _distance(first: AtomRecord, second: AtomRecord) -> float:
    return math.sqrt(
        sum((left - right) ** 2 for left, right in zip(first.coordinate, second.coordinate, strict=True))
    )


def _residue_label(record: ResidueRecord) -> str:
    residue_id = record.residue_id
    insertion = "" if residue_id.insertion_code is None else residue_id.insertion_code.strip()
    return f"{residue_id.structure_id}:{residue_id.model_id}:{residue_id.chain_id}:{residue_id.auth_seq_id}{insertion}:{residue_id.residue_name}"


def _source_label(chain: ProteinChain) -> str:
    return f"{chain.structure_id}:{chain.model_id}:{chain.chain_id}"


def _distance_diagnostic(
    chain: ProteinChain,
    first: ResidueRecord,
    second: ResidueRecord,
    *,
    first_atom_name: str,
    second_atom_name: str,
    minimum: float,
    maximum: float,
    unavailable_code: str,
    outlier_code: str,
    geometry_name: str,
) -> Diagnostic | None:
    first_atom = _atom(first, first_atom_name)
    second_atom = _atom(second, second_atom_name)
    residue_label = _residue_label(first)
    source_label = _source_label(chain)
    if first_atom is None or second_atom is None:
        missing = first_atom_name if first_atom is None else second_atom_name
        position = "i" if first_atom is None else "i+1"
        return Diagnostic(
            code=unavailable_code,
            severity=DiagnosticSeverity.WARNING,
            message=(
                f"{geometry_name} unavailable for {_residue_label(first)} to {_residue_label(second)}: "
                f"missing {missing}({position}) atom."
            ),
            source_id=source_label,
            residue_id=residue_label,
        )
    distance = _distance(first_atom, second_atom)
    if minimum <= distance <= maximum:
        return None
    return Diagnostic(
        code=outlier_code,
        severity=DiagnosticSeverity.WARNING,
        message=(
            f"{geometry_name} outside screening interval for {_residue_label(first)} to "
            f"{_residue_label(second)}: {distance:.3f} Å not in [{minimum:.3f}, {maximum:.3f}] Å."
        ),
        source_id=source_label,
        residue_id=residue_label,
    )


def _screen_chain(
    chain: ProteinChain,
    *,
    first_atom_name: str,
    second_atom_name: str,
    minimum: float,
    maximum: float,
    unavailable_code: str,
    outlier_code: str,
    geometry_name: str,
) -> tuple[Diagnostic, ...]:
    diagnostics: list[Diagnostic] = []
    for first, second in _pairs(chain):
        diagnostic = _distance_diagnostic(
            chain,
            first,
            second,
            first_atom_name=first_atom_name,
            second_atom_name=second_atom_name,
            minimum=minimum,
            maximum=maximum,
            unavailable_code=unavailable_code,
            outlier_code=outlier_code,
            geometry_name=geometry_name,
        )
        if diagnostic is not None:
            diagnostics.append(diagnostic)
    return tuple(diagnostics)


def screen_peptide_continuity(
    chains: ChainInput,
    *,
    min_cn_distance_angstrom: float = DEFAULT_MIN_CN_DISTANCE_ANGSTROM,
    max_cn_distance_angstrom: float = DEFAULT_MAX_CN_DISTANCE_ANGSTROM,
) -> tuple[Diagnostic, ...]:
    """Screen adjacent peptide C(i)-N(i+1) distances in Å.

    The bounds are inclusive and are intended for coordinate sanity screening.
    A source-numbering gap, unknown number, reversed order, or incompatible
    structure/model/chain identity is not treated as a peptide pair.
    """

    minimum, maximum = _validate_limits(
        min_cn_distance_angstrom,
        max_cn_distance_angstrom,
        name="C-N distance",
    )
    diagnostics: list[Diagnostic] = []
    for chain in _chains(chains):
        diagnostics.extend(
            _screen_chain(
                chain,
                first_atom_name="C",
                second_atom_name="N",
                minimum=minimum,
                maximum=maximum,
                unavailable_code="quality.geometry.peptide_cn_unavailable",
                outlier_code="quality.geometry.peptide_cn_gap",
                geometry_name="Peptide observed-coordinate discontinuity C(i)-N(i+1) distance",
            )
        )
    return tuple(diagnostics)


def screen_ca_pseudo_geometry(
    chains: ChainInput,
    *,
    min_ca_distance_angstrom: float = DEFAULT_MIN_CA_DISTANCE_ANGSTROM,
    max_ca_distance_angstrom: float = DEFAULT_MAX_CA_DISTANCE_ANGSTROM,
) -> tuple[Diagnostic, ...]:
    """Screen adjacent Cα(i)-Cα(i+1) pseudo-geometry distances in Å.

    This is a translation/rotation-invariant distance screen, not a
    crystallographic validation or a claim about backbone conformation.
    """

    minimum, maximum = _validate_limits(
        min_ca_distance_angstrom,
        max_ca_distance_angstrom,
        name="C-alpha distance",
    )
    diagnostics: list[Diagnostic] = []
    for chain in _chains(chains):
        diagnostics.extend(
            _screen_chain(
                chain,
                first_atom_name="CA",
                second_atom_name="CA",
                minimum=minimum,
                maximum=maximum,
                unavailable_code="quality.geometry.ca_distance_unavailable",
                outlier_code="quality.geometry.ca_distance_outlier",
                geometry_name="C-alpha(i)-C-alpha(i+1) distance",
            )
        )
    return tuple(diagnostics)


def screen_chain_geometry(
    chains: ChainInput,
    *,
    min_cn_distance_angstrom: float = DEFAULT_MIN_CN_DISTANCE_ANGSTROM,
    max_cn_distance_angstrom: float = DEFAULT_MAX_CN_DISTANCE_ANGSTROM,
    min_ca_distance_angstrom: float = DEFAULT_MIN_CA_DISTANCE_ANGSTROM,
    max_ca_distance_angstrom: float = DEFAULT_MAX_CA_DISTANCE_ANGSTROM,
) -> tuple[Diagnostic, ...]:
    """Run both deterministic geometry screens, in stable screen order."""

    return screen_peptide_continuity(
        chains,
        min_cn_distance_angstrom=min_cn_distance_angstrom,
        max_cn_distance_angstrom=max_cn_distance_angstrom,
    ) + screen_ca_pseudo_geometry(
        chains,
        min_ca_distance_angstrom=min_ca_distance_angstrom,
        max_ca_distance_angstrom=max_ca_distance_angstrom,
    )


__all__ = [
    "ChainInput",
    "DEFAULT_MAX_CA_DISTANCE_ANGSTROM",
    "DEFAULT_MAX_CN_DISTANCE_ANGSTROM",
    "DEFAULT_MIN_CA_DISTANCE_ANGSTROM",
    "DEFAULT_MIN_CN_DISTANCE_ANGSTROM",
    "screen_ca_pseudo_geometry",
    "screen_chain_geometry",
    "screen_peptide_continuity",
]
