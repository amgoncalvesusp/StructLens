"""Source-row metadata and bounded preflight counters for coordinate inputs."""

from __future__ import annotations

import io
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from .limits import ParseLimits, SnapshotLimitError, StructureParseError


@dataclass(frozen=True, slots=True)
class AtomSiteMetadata:
    """Fields unavailable on a generic Biopython ``Atom`` object."""

    source_id: str
    author_chain_id: str
    label_chain_id: str | None
    entity_id: str | None
    label_seq_id: str | None
    auth_seq_id: str
    insertion_code: str | None
    residue_name: str
    formal_charge: int | None
    group: str


def preflight_pdb(source: bytes, limits: ParseLimits) -> None:
    """Count PDB records from bounded bytes before constructing Bio objects."""

    models: set[str] = set()
    chains: set[tuple[str, str]] = set()
    residues: set[tuple[str, str, str, str, str, str]] = set()
    atom_ids: set[tuple[str, str]] = set()
    atom_count = 0
    current_model = "0"
    saw_model = False
    model_number = 0
    stream = io.BytesIO(source)
    while raw_line := stream.readline():
        line = raw_line.decode("utf-8", errors="replace").rstrip("\r\n")
        record = line[:6].strip().upper()
        if record == "MODEL":
            saw_model = True
            model_number += 1
            current_model = _canonical_pdb_integer(
                line[10:14].strip() or str(model_number), "MODEL serial"
            )
            models.add(current_model)
            limits.validate_counts(models=len(models), chains=len(chains), residues=len(residues), atoms=atom_count)
            continue
        if record not in {"ATOM", "HETATM"}:
            continue
        if not saw_model:
            models.add(current_model)
            limits.validate_counts(models=len(models))
        serial_text = line[6:11].strip()
        if not serial_text:
            raise StructureParseError("malformed PDB atom record: missing atom serial")
        try:
            serial = str(int(serial_text))
        except ValueError as error:
            raise StructureParseError(
                f"malformed PDB atom record: atom serial {serial_text!r} is not preservable"
            ) from error
        atom_key = (current_model, serial)
        if atom_key in atom_ids:
            raise StructureParseError(f"duplicate PDB atom serial {serial!r} in model {current_model!r}")
        atom_ids.add(atom_key)
        chain = line[21:22].strip()
        residue = line[17:20].strip().upper()
        auth_seq = line[22:26].strip()
        insertion = line[26:27].strip()
        chains.add((current_model, chain))
        residues.add((current_model, chain, auth_seq, insertion, residue, record))
        atom_count += 1
        limits.validate_counts(
            models=len(models),
            chains=len(chains),
            residues=len(residues),
            atoms=atom_count,
        )


def preflight_mmcif(cif: Mapping[str, object], limits: ParseLimits) -> None:
    """Count atom-site rows from a parsed mmCIF table before Bio objects."""

    ids = column(cif, "_atom_site.id")
    models = column(cif, "_atom_site.pdbx_PDB_model_num", len(ids))
    authors = column(cif, "_atom_site.auth_asym_id", len(ids))
    labels = column(cif, "_atom_site.label_asym_id", len(ids))
    auth_seq = column(cif, "_atom_site.auth_seq_id", len(ids))
    label_seq = column(cif, "_atom_site.label_seq_id", len(ids))
    insertion = column(cif, "_atom_site.pdbx_PDB_ins_code", len(ids))
    auth_names = column(cif, "_atom_site.auth_comp_id", len(ids))
    label_names = column(cif, "_atom_site.label_comp_id", len(ids))
    names = auth_names if any(item.strip() for item in auth_names) else label_names
    model_values: set[str] = set()
    chain_values: set[tuple[str, str]] = set()
    residue_values: set[tuple[str, str, str, str, str]] = set()
    for index in range(len(ids)):
        model = _canonical_mmcif_model(value(models, index, "1"))
        chain = value(authors, index, value(labels, index, ""))
        sequence = value(auth_seq, index, value(label_seq, index, ""))
        code = value(insertion, index, "")
        name = value(names, index, "").upper()
        model_values.add(model)
        chain_values.add((model, chain))
        residue_values.add((model, chain, sequence, code, name))
    limits.validate_counts(
        models=len(model_values), chains=len(chain_values), residues=len(residue_values), atoms=len(ids)
    )


def preflight_mmcif_bytes(source: bytes, limits: ParseLimits) -> None:
    """Tokenize an mmCIF stream and enforce row limits before materialization.

    ``MMCIF2Dict`` is intentionally not used here: its normal constructor
    stores every loop value. The actual parser may still use it after this
    bounded pass for metadata extraction.
    """

    try:
        text = source.decode("utf-8")
        tokenizer = MMCIF2Dict.__new__(MMCIF2Dict)
        tokenizer.quote_chars = ["'", '"']
        tokenizer.whitespace_chars = [" ", "\t"]
        stream = _TokenStream(iter(tokenizer._tokenize(io.StringIO(text))))  # type: ignore[no-untyped-call]
        state: dict[str, set[Any]] = {}
        while True:
            token = stream.next()
            if token is None:
                return
            if token != "loop_":
                continue
            keys: list[str] = []
            first_value: str | None = None
            while True:
                token = stream.next()
                if token is None:
                    break
                if token.startswith("_"):
                    keys.append(token)
                    continue
                first_value = token
                break
            if not keys:
                raise StructureParseError("malformed mmCIF loop: no column names")
            if any(key.startswith("_atom_site.") for key in keys):
                if "_atom_site.id" not in keys:
                    raise StructureParseError("malformed mmCIF atom_site loop: missing _atom_site.id")
                _consume_atom_rows(stream, keys, first_value, limits, state)
            else:
                _consume_ignored_loop(stream, keys, first_value)
    except (SnapshotLimitError, StructureParseError):
        raise
    except (UnicodeDecodeError, ValueError, TypeError, AttributeError) as error:
        raise StructureParseError(f"malformed mmCIF token stream: {error}") from error


class _TokenStream:
    def __init__(self, tokens: Any) -> None:
        self._tokens = iter(tokens)
        self._pending: str | None = None

    def next(self) -> str | None:
        if self._pending is not None:
            token = self._pending
            self._pending = None
            return token
        try:
            return str(next(self._tokens))
        except StopIteration:
            return None

    def push(self, token: str) -> None:
        if self._pending is not None:
            raise StructureParseError("malformed mmCIF token stream: duplicate pushback")
        self._pending = token


def _consume_atom_rows(
    stream: _TokenStream,
    keys: Sequence[str],
    first_value: str | None,
    limits: ParseLimits,
    state: dict[str, set[Any]],
) -> None:
    index = {key: position for position, key in enumerate(keys)}
    row: list[str] = []
    if first_value is not None:
        row.append(first_value)
    while True:
        if len(row) == len(keys):
            _record_atom_row(row, index, limits, state)
            row = []
            continue
        token = stream.next()
        if token is None:
            if row:
                if len(row) != len(keys):
                    raise StructureParseError("malformed mmCIF atom_site loop: incomplete row")
                _record_atom_row(row, index, limits, state)
            return
        if _is_loop_boundary(token):
            if row:
                raise StructureParseError("malformed mmCIF atom_site loop: incomplete row")
            stream.push(token)
            return
        row.append(token)


def _record_atom_row(
    row: Sequence[str],
    index: Mapping[str, int],
    limits: ParseLimits,
    state: dict[str, set[Any]],
) -> None:
    source_id = _row_value(row, index, "_atom_site.id")
    if source_id in {"", ".", "?"}:
        raise StructureParseError("malformed mmCIF atom_site row: missing atom-site id")
    model = _canonical_mmcif_model(
        _row_value(row, index, "_atom_site.pdbx_PDB_model_num") or "1"
    )
    author = _row_value(row, index, "_atom_site.auth_asym_id")
    label = _row_value(row, index, "_atom_site.label_asym_id")
    chain = author or label
    sequence = _row_value(row, index, "_atom_site.auth_seq_id") or _row_value(row, index, "_atom_site.label_seq_id")
    residue_name = _row_value(row, index, "_atom_site.auth_comp_id") or _row_value(
        row, index, "_atom_site.label_comp_id"
    )
    insertion = _row_value(row, index, "_atom_site.pdbx_PDB_ins_code")
    ids = state.setdefault("atom_ids", set())
    key = (model, _canonical_mmcif_atom_id(source_id))
    if key in ids:
        raise StructureParseError(f"duplicate mmCIF atom_site id {source_id!r} in model {model!r}")
    ids.add(key)
    models = state.setdefault("models", set())
    chains = state.setdefault("chains", set())
    residues = state.setdefault("residues", set())
    models.add(model)
    chains.add((model, chain))
    residues.add((model, chain, sequence, insertion, residue_name))
    limits.validate_counts(models=len(models), chains=len(chains), residues=len(residues), atoms=len(ids))


def _consume_ignored_loop(
    stream: _TokenStream,
    keys: Sequence[str],
    first_value: str | None,
) -> None:
    row_length = len(keys)
    row_size = 1 if first_value is not None else 0
    while True:
        if row_size == row_length:
            row_size = 0
            continue
        token = stream.next()
        if token is None:
            if row_size:
                raise StructureParseError("malformed mmCIF loop: incomplete row")
            return
        if _is_loop_boundary(token):
            if row_size:
                raise StructureParseError("malformed mmCIF loop: incomplete row")
            stream.push(token)
            return
        row_size += 1


def _row_value(row: Sequence[str], index: Mapping[str, int], key: str) -> str:
    position = index.get(key)
    if position is None or position >= len(row):
        return ""
    value_text = row[position].strip()
    return "" if value_text in {".", "?"} else value_text


def _is_loop_boundary(token: str) -> bool:
    return token == "loop_" or token.startswith("data_") or token.startswith("_")


def pdb_atom_metadata(source: bytes) -> dict[tuple[str, str], AtomSiteMetadata]:
    """Extract PDB source IDs and formal charges without losing identity."""

    result: dict[tuple[str, str], AtomSiteMetadata] = {}
    current_model = "0"
    saw_model = False
    for line in source.decode("utf-8", errors="replace").splitlines():
        record = line[:6].strip().upper()
        if record == "MODEL":
            saw_model = True
            current_model = _canonical_pdb_integer(
                line[10:14].strip() or current_model, "MODEL serial"
            )
            continue
        if record not in {"ATOM", "HETATM"}:
            continue
        serial_text = line[6:11].strip()
        if not serial_text:
            raise StructureParseError("malformed PDB atom record: missing atom serial")
        serial = _canonical_pdb_integer(serial_text, "atom serial")
        author = line[21:22].strip()
        auth_seq = line[22:26].strip()
        result[(current_model, serial)] = AtomSiteMetadata(
            source_id=serial,
            author_chain_id=author,
            label_chain_id=None,
            entity_id=None,
            label_seq_id=None,
            auth_seq_id=auth_seq,
            insertion_code=clean_optional(line[26:27]),
            residue_name=line[17:20].strip().upper(),
            formal_charge=formal_charge(line[78:80]) if len(line) >= 80 else None,
            group=record,
        )
        if not saw_model:
            current_model = "0"
    return result


def mmcif_atom_metadata(cif: Mapping[str, object]) -> dict[tuple[str, str], AtomSiteMetadata]:
    ids = column(cif, "_atom_site.id")
    count = len(ids)
    values = {
        key: column(cif, key, count)
        for key in (
            "_atom_site.pdbx_PDB_model_num",
            "_atom_site.auth_asym_id",
            "_atom_site.label_asym_id",
            "_atom_site.label_entity_id",
            "_atom_site.label_seq_id",
            "_atom_site.auth_seq_id",
            "_atom_site.pdbx_PDB_ins_code",
            "_atom_site.label_comp_id",
            "_atom_site.auth_comp_id",
            "_atom_site.pdbx_formal_charge",
            "_atom_site.group_PDB",
        )
    }
    result: dict[tuple[str, str], AtomSiteMetadata] = {}
    for index, source_id in enumerate(ids):
        model = _canonical_mmcif_model(
            value(values["_atom_site.pdbx_PDB_model_num"], index, "1")
        )
        author = value(values["_atom_site.auth_asym_id"], index, "")
        label = clean_optional(value(values["_atom_site.label_asym_id"], index, ""))
        auth_seq = value(values["_atom_site.auth_seq_id"], index, value(values["_atom_site.label_seq_id"], index, ""))
        residue_name = value(
            values["_atom_site.auth_comp_id"], index, value(values["_atom_site.label_comp_id"], index, "")
        ).upper()
        result[(model, _canonical_mmcif_atom_id(source_id))] = AtomSiteMetadata(
            source_id=source_id,
            author_chain_id=author,
            label_chain_id=label,
            entity_id=clean_optional(value(values["_atom_site.label_entity_id"], index, "")),
            label_seq_id=clean_optional(value(values["_atom_site.label_seq_id"], index, "")),
            auth_seq_id=auth_seq,
            insertion_code=clean_optional(value(values["_atom_site.pdbx_PDB_ins_code"], index, "")),
            residue_name=residue_name,
            formal_charge=formal_charge(value(values["_atom_site.pdbx_formal_charge"], index, "")),
            group=value(values["_atom_site.group_PDB"], index, "ATOM").upper(),
        )
    return result


def column(cif: Mapping[str, object], key: str, count: int | None = None) -> list[str]:
    raw = cif.get(key, [])
    if isinstance(raw, str):
        values = [raw]
    elif isinstance(raw, Sequence):
        values = [str(item) for item in raw]
    else:
        values = []
    if count is not None and len(values) < count:
        values.extend([""] * (count - len(values)))
    return values


def value(values: Sequence[str], index: int, default: str) -> str:
    if index >= len(values):
        return default
    item = values[index].strip()
    return default if item in {"", ".", "?"} else item


def first_value(cif: Mapping[str, object], key: str) -> str | None:
    values = column(cif, key)
    return clean_optional(values[0]) if values else None


def formal_charge(value_text: str) -> int | None:
    value_text = value_text.strip()
    if value_text in {"", ".", "?"}:
        return None
    try:
        return int(value_text)
    except ValueError:
        if value_text[-1:] in {"+", "-"} and value_text[:-1].isdigit():
            amount = int(value_text[:-1])
            return amount if value_text[-1] == "+" else -amount
        return None


def clean_optional(value_text: str) -> str | None:
    value_text = value_text.strip().strip("'\"")
    return None if value_text in {"", ".", "?"} else value_text


def clean_method(value: object) -> str | None:
    if value is None:
        return None
    text = str(value).strip().strip("'\"")
    return text.upper() if text and text not in {".", "?"} else None


def _canonical_pdb_integer(value: str, field_name: str) -> str:
    try:
        return str(int(value))
    except ValueError as error:
        raise StructureParseError(
            f"malformed PDB record: {field_name} {value!r} is not preservable"
        ) from error


def _canonical_mmcif_model(value: str) -> str:
    try:
        return str(int(value))
    except ValueError as error:
        raise StructureParseError(
            f"malformed mmCIF atom_site row: model number {value!r} is not an integer"
        ) from error


def _canonical_mmcif_atom_id(value: str) -> str:
    try:
        return str(int(value))
    except ValueError:
        # Biopython also retains non-numeric mmCIF atom-site IDs as strings.
        return value


__all__ = [
    "AtomSiteMetadata",
    "clean_optional",
    "clean_method",
    "column",
    "first_value",
    "formal_charge",
    "mmcif_atom_metadata",
    "pdb_atom_metadata",
    "preflight_mmcif",
    "preflight_mmcif_bytes",
    "preflight_pdb",
    "value",
]
