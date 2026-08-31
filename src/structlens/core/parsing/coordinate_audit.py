"""Raw coordinate quality gates for PDB and mmCIF sources.

This module runs before Biopython's strict structure constructors.  It keeps a
lenient, immutable copy of every ``ATOM``/``HETATM`` row so malformed numeric
evidence can be reported with its source context instead of being coerced or
silently discarded.  The normalizer may proceed only when the report has no
error diagnostics.
"""

from __future__ import annotations

from collections import OrderedDict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from structlens.core.models import ResidueId
from structlens.core.quality.coordinates import check_coordinate_quality
from structlens.core.quality.models import (
    CoordinateAtom,
    CoordinateQCSettings,
    CoordinateResidue,
    StructureQualityReport,
)

from .limits import ParseLimits, SnapshotLimitError, StructureParseError
from .source_records import column, formal_charge, preflight_pdb, value


class CoordinateQualityError(StructureParseError):
    """Raised when raw coordinate evidence cannot safely reach a parser.

    ``report`` is deliberately stored by identity, not copied, allowing a
    caller to inspect exactly the report that caused parser short-circuiting.
    """

    def __init__(self, report: StructureQualityReport) -> None:
        if not isinstance(report, StructureQualityReport):
            raise TypeError("report must be a StructureQualityReport")
        self.report = report
        self.quality_report = report
        errors = report.error_count
        super().__init__(f"coordinate quality gate rejected source ({errors} error diagnostic(s))")


@dataclass(frozen=True, slots=True)
class CoordinateAudit:
    """Lossless raw rows and the typed report produced for them."""

    residues: tuple[CoordinateResidue, ...]
    report: StructureQualityReport

    def __post_init__(self) -> None:
        residues = tuple(self.residues)
        if any(not isinstance(residue, CoordinateResidue) for residue in residues):
            raise TypeError("residues must contain CoordinateResidue values")
        if not isinstance(self.report, StructureQualityReport):
            raise TypeError("report must be a StructureQualityReport")
        object.__setattr__(self, "residues", residues)

    @property
    def can_normalize(self) -> bool:
        """Whether strict normalization may proceed without losing evidence."""

        return self.report.error_count == 0


def audit_pdb_source(
    source: bytes | bytearray | memoryview | Path,
    *,
    source_id: str | None = None,
    structure_id: str = "source",
    limits: ParseLimits | None = None,
    settings: CoordinateQCSettings | None = None,
    raise_on_error: bool = True,
) -> CoordinateAudit:
    """Audit fixed-column PDB rows before invoking ``PDBParser``."""

    active_limits = limits or ParseLimits()
    raw = _read_source(source, active_limits)
    preflight_pdb(raw, active_limits)
    residues = _pdb_residues(raw, structure_id)
    report = check_coordinate_quality(residues, source_id=source_id, settings=settings)
    _validate_audit_limits(residues, raw, active_limits)
    audit = CoordinateAudit(residues, report)
    _raise_if_needed(audit, raise_on_error)
    return audit


def audit_mmcif_source(
    cif: Mapping[str, object],
    *,
    source_id: str | None = None,
    structure_id: str = "source",
    limits: ParseLimits | None = None,
    settings: CoordinateQCSettings | None = None,
    raise_on_error: bool = True,
) -> CoordinateAudit:
    """Audit an already tokenized ``_atom_site`` table before ``MMCIFParser``."""

    residues = _mmcif_residues(cif, structure_id)
    report = check_coordinate_quality(residues, source_id=source_id, settings=settings)
    _validate_audit_limits(residues, b"", limits or ParseLimits())
    audit = CoordinateAudit(residues, report)
    _raise_if_needed(audit, raise_on_error)
    return audit


# Explicit aliases make the parser boundary readable at call sites and keep
# one implementation for PDB/mmCIF direct-loader integrations.
audit_pdb_coordinates = audit_pdb_source
audit_mmcif_coordinates = audit_mmcif_source


def _raise_if_needed(audit: CoordinateAudit, raise_on_error: bool) -> None:
    if not isinstance(raise_on_error, bool):
        raise TypeError("raise_on_error must be a boolean")
    if raise_on_error and audit.report.error_count:
        raise CoordinateQualityError(audit.report)


def _read_source(
    source: bytes | bytearray | memoryview | Path,
    limits: ParseLimits,
) -> bytes:
    if isinstance(source, Path):
        size = source.stat().st_size
        if size > limits.max_raw_bytes:
            raise SnapshotLimitError("raw_bytes", limits.max_raw_bytes, size)
        with source.open("rb") as stream:
            raw = stream.read(limits.max_raw_bytes + 1)
        if len(raw) > limits.max_raw_bytes:
            raise SnapshotLimitError("raw_bytes", limits.max_raw_bytes, len(raw))
        return raw
    if isinstance(source, bytes):
        raw = source
    elif isinstance(source, (bytearray, memoryview)):
        raw = bytes(source)
    else:
        raise TypeError("source must be bytes or pathlib.Path")
    if len(raw) > limits.max_decompressed_bytes:
        raise SnapshotLimitError("decompressed_bytes", limits.max_decompressed_bytes, len(raw))
    return raw


def _validate_audit_limits(
    residues: Sequence[CoordinateResidue], source: bytes, limits: ParseLimits
) -> None:
    if not isinstance(limits, ParseLimits):
        raise TypeError("limits must be ParseLimits")
    models = {residue.residue_id.model_id for residue in residues}
    chains = {(residue.residue_id.model_id, residue.residue_id.chain_id) for residue in residues}
    atom_count = sum(len(residue.atoms) for residue in residues)
    limits.validate_counts(
        models=len(models),
        chains=len(chains),
        residues=len(residues),
        atoms=atom_count,
    )
    if source:
        limits.validate_counts(atoms=atom_count)


def _pdb_residues(source: bytes, structure_id: str) -> tuple[CoordinateResidue, ...]:
    grouped: OrderedDict[tuple[str, str, str, str, str], list[CoordinateAtom]] = OrderedDict()
    polymer: dict[tuple[str, str, str, str, str], bool] = {}
    try:
        text = source.decode("utf-8", errors="replace")
    except AttributeError as error:
        raise StructureParseError("unable to inspect PDB source") from error
    current_model = "0"
    model_counter = 0
    for line in text.splitlines():
        record = line[:6].strip().upper()
        if record == "MODEL":
            model_counter += 1
            current_model = _canonical_integer(line[10:14].strip() or str(model_counter), "MODEL serial")
            continue
        if record not in {"ATOM", "HETATM"}:
            continue
        serial = line[6:11].strip() or None
        chain = line[21:22].strip()
        auth_seq = line[22:26].strip()
        insertion = line[26:27].strip()
        residue_name = line[17:20].strip().upper()
        key = (current_model, chain, auth_seq, insertion, residue_name)
        grouped.setdefault(key, []).append(
            CoordinateAtom(
                name=line[12:16].strip(),
                element=line[76:78].strip(),
                coordinate=(_number(line[30:38]), _number(line[38:46]), _number(line[46:54])),
                altloc=_optional(line[16:17]),
                occupancy=_optional_number(line[54:60]),
                b_factor=_optional_number(line[60:66]),
                source_atom_id=serial,
                source_serial=serial,
                formal_charge=formal_charge(line[78:80]) if len(line) >= 80 else None,
            )
        )
        polymer[key] = polymer.get(key, True) and record == "ATOM"
    residues = [
        CoordinateResidue(
            residue_id=_residue_id(structure_id, key),
            atoms=tuple(atoms),
            is_polymer=polymer[key],
        )
        for key, atoms in grouped.items()
    ]
    return tuple(residues)


def _mmcif_residues(cif: Mapping[str, object], structure_id: str) -> tuple[CoordinateResidue, ...]:
    ids = column(cif, "_atom_site.id")
    count = len(ids)
    fields = {
        key: column(cif, key, count)
        for key in (
            "_atom_site.group_PDB",
            "_atom_site.type_symbol",
            "_atom_site.label_atom_id",
            "_atom_site.auth_atom_id",
            "_atom_site.label_alt_id",
            "_atom_site.auth_comp_id",
            "_atom_site.label_comp_id",
            "_atom_site.auth_asym_id",
            "_atom_site.label_asym_id",
            "_atom_site.auth_seq_id",
            "_atom_site.label_seq_id",
            "_atom_site.pdbx_PDB_ins_code",
            "_atom_site.Cartn_x",
            "_atom_site.Cartn_y",
            "_atom_site.Cartn_z",
            "_atom_site.occupancy",
            "_atom_site.B_iso_or_equiv",
            "_atom_site.pdbx_PDB_model_num",
            "_atom_site.pdbx_formal_charge",
        )
    }
    grouped: OrderedDict[tuple[str, str, str, str, str], list[CoordinateAtom]] = OrderedDict()
    polymer: dict[tuple[str, str, str, str, str], bool] = {}
    for index, source_id in enumerate(ids):
        model = _canonical_integer(value(fields["_atom_site.pdbx_PDB_model_num"], index, "1"), "model number")
        author = value(fields["_atom_site.auth_asym_id"], index, "")
        label = value(fields["_atom_site.label_asym_id"], index, "")
        chain = author or label
        auth_seq = value(fields["_atom_site.auth_seq_id"], index, "")
        if not auth_seq:
            auth_seq = value(fields["_atom_site.label_seq_id"], index, "")
        insertion = value(fields["_atom_site.pdbx_PDB_ins_code"], index, "")
        residue_name = value(fields["_atom_site.auth_comp_id"], index, "")
        if not residue_name:
            residue_name = value(fields["_atom_site.label_comp_id"], index, "")
        key = (model, chain, auth_seq, insertion, residue_name.upper())
        name = value(fields["_atom_site.auth_atom_id"], index, "")
        if not name:
            name = value(fields["_atom_site.label_atom_id"], index, "")
        atom = CoordinateAtom(
            name=name,
            element=value(fields["_atom_site.type_symbol"], index, ""),
            coordinate=(
                _number(value(fields["_atom_site.Cartn_x"], index, "")),
                _number(value(fields["_atom_site.Cartn_y"], index, "")),
                _number(value(fields["_atom_site.Cartn_z"], index, "")),
            ),
            altloc=_optional(value(fields["_atom_site.label_alt_id"], index, "")),
            occupancy=_optional_number(value(fields["_atom_site.occupancy"], index, "")),
            b_factor=_optional_number(value(fields["_atom_site.B_iso_or_equiv"], index, "")),
            source_atom_id=source_id,
            source_serial=_source_serial(source_id),
            formal_charge=formal_charge(value(fields["_atom_site.pdbx_formal_charge"], index, "")),
        )
        grouped.setdefault(key, []).append(atom)
        group = value(fields["_atom_site.group_PDB"], index, "ATOM").upper()
        polymer[key] = polymer.get(key, True) and group == "ATOM"
    return tuple(
        CoordinateResidue(_residue_id(structure_id, key), tuple(atoms), is_polymer=polymer[key])
        for key, atoms in grouped.items()
    )


def _residue_id(structure_id: str, key: tuple[str, str, str, str, str]) -> ResidueId:
    model, chain, auth_seq, insertion, residue_name = key
    return ResidueId(
        structure_id,
        model,
        chain,
        auth_seq,
        insertion or None,
        residue_name,
    )


def _number(value_text: str) -> object:
    text = str(value_text).strip()
    if not text:
        return ""
    try:
        return float(text)
    except (TypeError, ValueError, OverflowError):
        return text


def _optional_number(value_text: str) -> object | None:
    text = str(value_text).strip()
    if not text or text in {".", "?"}:
        return None
    return _number(text)


def _optional(value_text: str) -> str | None:
    text = str(value_text).strip()
    return None if not text or text in {".", "?"} else text


def _source_serial(value_text: str) -> int | str:
    try:
        return int(value_text)
    except (TypeError, ValueError):
        return value_text


def _canonical_integer(value_text: str, field_name: str) -> str:
    try:
        return str(int(value_text))
    except (TypeError, ValueError) as error:
        raise StructureParseError(f"malformed coordinate source: {field_name} {value_text!r} is not an integer") from error


__all__ = [
    "CoordinateAudit",
    "CoordinateQualityError",
    "audit_mmcif_coordinates",
    "audit_mmcif_source",
    "audit_pdb_coordinates",
    "audit_pdb_source",
]
