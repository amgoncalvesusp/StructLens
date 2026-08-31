"""Content-addressed coordinate parsing and lossless normalization.

The parser boundary is deliberately the only place that knows about
Biopython. A source is captured into :class:`SourceSnapshot` first and the
Biopython adapters receive a ``StringIO`` over those bytes; they never reopen
the caller's path.
"""

from __future__ import annotations

import io
import warnings
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser

from structlens.core.evidence.status import Diagnostic, DiagnosticSeverity
from structlens.core.models import (
    AtomRecord,
    ComponentKind,
    ProteinChain,
    ProteinStructure,
    ResidueId,
    ResidueNumbering,
    ResidueRecord,
    StructureComponent,
)

from .classification import (
    ION_NAMES as _ION_NAMES,
)
from .classification import (
    KNOWN_LIGANDS as _KNOWN_LIGANDS,
)
from .classification import (
    LIGAND_TABLE_VERSION as _LIGAND_TABLE_VERSION,
)
from .classification import (
    MODIFIED_POLYMER_NAMES as _MODIFIED_POLYMER_NAMES,
)
from .classification import (
    WATER_NAMES as _WATER_NAMES,
)
from .classification import (
    altloc_inventory as _altloc_inventory,
)
from .classification import (
    chain_selected as _chain_selected,
)
from .classification import (
    classify_residue as _classify_residue,
)
from .classification import (
    component_id as _component_id,
)
from .classification import (
    selected_atoms as _selected_atoms,
)
from .limits import ParseLimits, SnapshotLimitError, StructureParseError
from .models import (
    AltlocPolicy,
    AssemblyScope,
    InputSelection,
    ParsedStructure,
    StructureFormat,
    StructureMetadata,
)
from .snapshot import SourceSnapshot, capture_snapshot
from .source_records import (
    AtomSiteMetadata as _AtomSiteMetadata,
)
from .source_records import (
    clean_method as _clean_method,
)
from .source_records import (
    clean_optional as _clean_optional,
)
from .source_records import (
    first_value as _first_value,
)
from .source_records import (
    mmcif_atom_metadata as _mmcif_atom_metadata,
)
from .source_records import (
    pdb_atom_metadata as _pdb_atom_metadata,
)
from .source_records import (
    preflight_mmcif_bytes as _preflight_mmcif_bytes,
)
from .source_records import (
    preflight_pdb as _preflight_pdb,
)

_ONE_LETTER_BY_THREE_LETTER = {name.upper(): letter for name, letter in protein_letters_3to1.items()}
KNOWN_LIGANDS = _KNOWN_LIGANDS
LIGAND_TABLE_VERSION = _LIGAND_TABLE_VERSION
ION_NAMES = _ION_NAMES


def load_structure(path: Path) -> ProteinStructure:
    """Load a PDB/mmCIF path using the compatibility return type."""

    return load_structure_legacy(path)


def load_structure_legacy(path: Path) -> ProteinStructure:
    """Load v0.3-compatible all-model and non-water HETATM records."""

    snapshot = capture_snapshot(path)
    structure_id = _structure_id(snapshot.display_name, snapshot.logical_format)
    if snapshot.logical_format == "pdb":
        _preflight_pdb(snapshot.decompressed_bytes, ParseLimits())
        structure, atom_metadata, _, _, _ = _parse_pdb(snapshot, structure_id)
    else:
        structure, atom_metadata, _, _, _ = _parse_mmcif(snapshot, structure_id, ParseLimits())
    return _normalize_legacy_structure(structure, structure_id, path, atom_metadata)


def load_structure_evidence(
    source: Path | SourceSnapshot,
    *,
    selection: InputSelection | None = None,
    limits: ParseLimits | None = None,
) -> ParsedStructure:
    """Parse one bounded source and retain normalized structural evidence."""

    active_limits = limits if limits is not None else ParseLimits()
    if isinstance(source, SourceSnapshot):
        snapshot = source
    elif isinstance(source, Path):
        snapshot = capture_snapshot(source, active_limits)
    else:
        raise TypeError("source must be a pathlib.Path or SourceSnapshot")
    source_path = str(source) if isinstance(source, Path) else None

    if selection is not None:
        if selection.assembly_scope is AssemblyScope.BIOLOGICAL_ASSEMBLY:
            raise ValueError("biological assembly parsing is not built; use asymmetric_unit")
        if selection.altloc_policy is AltlocPolicy.ALL:
            raise ValueError("altloc policy 'all' is not supported until conformer-aware analysis is built")
        _validate_selection(selection, snapshot)

    structure_id = _structure_id(snapshot.display_name, snapshot.logical_format)
    if snapshot.logical_format == "pdb":
        _preflight_pdb(snapshot.decompressed_bytes, active_limits)
        structure, atom_metadata, diagnostics, method, resolution = _parse_pdb(snapshot, structure_id)
        fmt = StructureFormat.PDB
    else:
        structure, atom_metadata, diagnostics, method, resolution = _parse_mmcif(snapshot, structure_id, active_limits)
        fmt = StructureFormat.MMCIF

    available_models = tuple(_model_id(model) for model in structure)
    active_limits.validate_counts(
        models=len(available_models),
        chains=_chain_count(structure),
        residues=_residue_count(structure),
        atoms=_atom_count(structure),
    )
    if not available_models:
        raise ValueError("structure contains no models")

    requested_model = selection.model_id if selection is not None else None
    if requested_model is None:
        if len(available_models) != 1:
            raise ValueError(
                f"model_id is required when the source contains multiple models ({', '.join(available_models)})"
            )
        requested_model = available_models[0]
    if requested_model not in available_models:
        raise ValueError(f"unknown model_id {requested_model!r}; available models: {available_models}")

    active_selection = _effective_selection(
        selection,
        snapshot,
        fmt,
        requested_model,
        structure,
        atom_metadata,
        source_path,
    )
    protein_structure, components, metadata = _normalize_structure(
        structure,
        structure_id,
        snapshot,
        active_selection,
        available_models,
        atom_metadata,
        method,
        resolution,
    )
    return ParsedStructure(
        protein_structure=protein_structure,
        components=components,
        selection=active_selection,
        metadata=metadata,
        raw_source_hash=snapshot.raw_sha256,
        diagnostics=diagnostics,
    )


def normalize_biopython_structure(
    structure: Any,
    path: Path,
    label_numbering: Mapping[tuple[str, str, str, str | None, str], str | None] | None = None,
) -> ProteinStructure:
    """Retain the v0.3 normalizer for callers using Biopython directly."""

    return _normalize_legacy_structure(
        structure,
        str(structure.id),
        path,
        {},
        label_numbering,
    )


def _parse_pdb(
    snapshot: SourceSnapshot,
    structure_id: str,
) -> tuple[Any, dict[tuple[str, str], _AtomSiteMetadata], tuple[Diagnostic, ...], str | None, float | None]:
    text = _decode_source(snapshot)
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            parser = PDBParser(QUIET=False, PERMISSIVE=True)  # type: ignore[no-untyped-call]
            structure = parser.get_structure(structure_id, io.StringIO(text))  # type: ignore[no-untyped-call]
    except (SnapshotLimitError, StructureParseError):
        raise
    except Exception as error:
        raise StructureParseError(f"unable to parse PDB source {snapshot.display_name}: {error}") from error
    header = getattr(structure, "header", {}) or {}
    method = _clean_method(header.get("structure_method"))
    resolution = _positive_float(header.get("resolution"))
    diagnostics = _warning_diagnostics(caught, snapshot.display_name)
    return structure, _pdb_atom_metadata(snapshot.decompressed_bytes), diagnostics, method, resolution


def _parse_mmcif(
    snapshot: SourceSnapshot,
    structure_id: str,
    limits: ParseLimits,
) -> tuple[Any, dict[tuple[str, str], _AtomSiteMetadata], tuple[Diagnostic, ...], str | None, float | None]:
    text = _decode_source(snapshot)
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            _preflight_mmcif_bytes(snapshot.decompressed_bytes, limits)
            cif = MMCIF2Dict(io.StringIO(text))  # type: ignore[no-untyped-call]
            parser = MMCIFParser(QUIET=False)  # type: ignore[no-untyped-call]
            structure = parser.get_structure(structure_id, io.StringIO(text))  # type: ignore[no-untyped-call]
    except (SnapshotLimitError, StructureParseError):
        raise
    except Exception as error:
        raise StructureParseError(f"unable to parse mmCIF source {snapshot.display_name}: {error}") from error
    atom_metadata = _mmcif_atom_metadata(cif)
    method = _clean_method(_first_value(cif, "_exptl.method"))
    resolution = _positive_float(_first_value(cif, "_refine.ls_d_res_high"))
    diagnostics = _warning_diagnostics(caught, snapshot.display_name)
    return structure, atom_metadata, diagnostics, method, resolution


def _normalize_structure(
    structure: Any,
    structure_id: str,
    snapshot: SourceSnapshot,
    selection: InputSelection,
    available_models: tuple[str, ...],
    atom_metadata: Mapping[tuple[str, str], _AtomSiteMetadata],
    method: str | None,
    resolution: float | None,
    label_numbering: Mapping[tuple[str, str, str, str | None, str], str | None] | None = None,
) -> tuple[ProteinStructure, tuple[StructureComponent, ...], StructureMetadata]:
    grouped_records: dict[tuple[str, str | None, str | None], list[ResidueRecord]] = {}
    components: list[StructureComponent] = []
    selected_model_object = next(model for model in structure if _model_id(model) == selection.model_id)
    for source_chain in selected_model_object:
        author_chain_id = str(source_chain.id).strip()
        for residue in source_chain:
            residue_name = str(residue.get_resname()).strip().upper()
            residue_atoms = _all_atoms(residue)
            source_meta = _first_atom_metadata(residue_atoms, selection.model_id, atom_metadata)
            auth_seq_id, insertion_code, label_seq_id, entity_id, group = _residue_identifiers(
                residue,
                source_meta,
                label_numbering,
                selection.model_id,
                author_chain_id,
            )
            selected_for_analysis = _chain_selected(
                selection, author_chain_id, source_meta.label_chain_id if source_meta else None, entity_id
            )
            kind = _classify_residue(residue_name, group, residue_atoms)
            residue_id = ResidueId(
                structure_id,
                selection.model_id,
                author_chain_id,
                auth_seq_id,
                insertion_code,
                residue_name,
            )
            primary_atoms = tuple(
                _atom_record(atom, selection.model_id, atom_metadata)
                for atom in _selected_atoms(residue, selection.altloc_policy)
            )
            all_atoms = tuple(_atom_record(atom, selection.model_id, atom_metadata) for atom in residue_atoms)
            component = StructureComponent(
                component_id=_component_id(
                    selection.model_id,
                    author_chain_id,
                    source_meta.label_chain_id if source_meta else None,
                    entity_id,
                    auth_seq_id,
                    insertion_code,
                    residue_name,
                ),
                kind=kind,
                atoms=all_atoms,
                metadata={
                    "altloc_inventory": _altloc_inventory(residue),
                    "classification_table_version": _LIGAND_TABLE_VERSION,
                    "source_group": group,
                    "selected_for_analysis": selected_for_analysis,
                    "retention_scope": "selected_model",
                },
                residue_id=residue_id,
                model_id=selection.model_id,
                author_chain_id=author_chain_id,
                label_chain_id=source_meta.label_chain_id if source_meta else None,
                entity_id=entity_id,
                residue_name=residue_name,
                auth_seq_id=auth_seq_id,
                insertion_code=insertion_code,
            )
            components.append(component)
            if not selected_for_analysis or kind is not ComponentKind.POLYMER_RESIDUE:
                continue
            one_letter = _ONE_LETTER_BY_THREE_LETTER.get(residue_name)
            record = ResidueRecord(
                residue_id=residue_id,
                numbering=ResidueNumbering(auth_seq_id, label_seq_id, insertion_code),
                residue_name=residue_name,
                one_letter=one_letter,
                atoms=primary_atoms,
                is_standard=one_letter is not None and residue_name not in _MODIFIED_POLYMER_NAMES,
            )
            key = (author_chain_id, source_meta.label_chain_id if source_meta else None, entity_id)
            grouped_records.setdefault(key, []).append(record)
    chains: list[ProteinChain] = []
    for (author_chain_id, label_chain_id, entity_id), records in grouped_records.items():
        chains.append(
            ProteinChain(
                structure_id=structure_id,
                model_id=selection.model_id,
                chain_id=author_chain_id,
                residues=tuple(record.residue_id for record in records),
                sequence="".join(record.one_letter or "X" for record in records),
                residue_records=tuple(records),
                source_path=selection.path,
                author_chain_id=author_chain_id,
                label_chain_id=label_chain_id,
                entity_id=entity_id,
            )
        )
    analyzed_auth = tuple(dict.fromkeys(chain.author_chain_id or chain.chain_id for chain in chains))
    analyzed_label = tuple(
        label for label in dict.fromkeys(chain.label_chain_id for chain in chains) if label is not None
    )
    analyzed_entities = tuple(entity for entity in dict.fromkeys(chain.entity_id for chain in chains) if entity is not None)
    metadata = StructureMetadata(
        format=selection.format,
        selected_model=selection.model_id,
        available_models=available_models,
        analyzed_chain_ids=analyzed_auth,
        experimental_method=method,
        resolution_angstrom=resolution,
        assembly_scope=selection.assembly_scope,
        author_chain_ids=analyzed_auth,
        label_chain_ids=analyzed_label,
        entity_ids=analyzed_entities,
    )
    protein_structure = ProteinStructure(
        structure_id=structure_id,
        chains=tuple(chains),
        source_path=selection.path,
        metadata={
            "format": selection.format.value,
            "selected_model": selection.model_id,
            "source_content_id": snapshot.content_id,
        },
    )
    return protein_structure, tuple(components), metadata


def _normalize_legacy_structure(
    structure: Any,
    structure_id: str,
    path: Path,
    atom_metadata: Mapping[tuple[str, str], _AtomSiteMetadata],
    label_numbering: Mapping[tuple[str, str, str, str | None, str], str | None] | None = None,
) -> ProteinStructure:
    """Normalize the deliberately lossy v0.3 chain/sequence view."""

    chains: list[ProteinChain] = []
    for model in structure:
        model_id = _model_id(model)
        for source_chain in model:
            chain_id = str(source_chain.id).strip()
            records: list[ResidueRecord] = []
            for residue in source_chain:
                residue_name = str(residue.get_resname()).strip().upper()
                if residue_name in _WATER_NAMES:
                    continue
                residue_atoms = _all_atoms(residue)
                source_meta = _first_atom_metadata(residue_atoms, model_id, atom_metadata)
                auth_seq_id, insertion_code, label_seq_id, _, _ = _residue_identifiers(
                    residue,
                    source_meta,
                    label_numbering,
                    model_id,
                    chain_id,
                )
                residue_id = ResidueId(
                    structure_id,
                    model_id,
                    chain_id,
                    auth_seq_id,
                    insertion_code,
                    residue_name,
                )
                atoms = tuple(
                    _atom_record(atom, model_id, atom_metadata)
                    for atom in _selected_atoms(residue, AltlocPolicy.HIGHEST_OCCUPANCY)
                )
                one_letter = _ONE_LETTER_BY_THREE_LETTER.get(residue_name)
                records.append(
                    ResidueRecord(
                        residue_id=residue_id,
                        numbering=ResidueNumbering(auth_seq_id, label_seq_id, insertion_code),
                        residue_name=residue_name,
                        one_letter=one_letter,
                        atoms=atoms,
                        is_standard=one_letter is not None and residue_name not in _MODIFIED_POLYMER_NAMES,
                    )
                )
            if records:
                chains.append(
                    ProteinChain(
                        structure_id=structure_id,
                        model_id=model_id,
                        chain_id=chain_id,
                        residues=tuple(record.residue_id for record in records),
                        sequence="".join(record.one_letter or "X" for record in records),
                        residue_records=tuple(records),
                        source_path=str(path),
                    )
                )
    return ProteinStructure(structure_id=structure_id, chains=tuple(chains), source_path=str(path))


def _effective_selection(
    selection: InputSelection | None,
    snapshot: SourceSnapshot,
    fmt: StructureFormat,
    model_id: str,
    structure: Any,
    atom_metadata: Mapping[tuple[str, str], _AtomSiteMetadata],
    source_path: str | None,
) -> InputSelection:
    if selection is not None:
        if selection.model_id is None or (selection.path is None and source_path is not None):
            return InputSelection(
                content_id=selection.content_id,
                display_name=selection.display_name,
                format=selection.format,
                model_id=model_id,
                author_chain_ids=selection.author_chain_ids,
                label_chain_ids=selection.label_chain_ids,
                altloc_policy=selection.altloc_policy,
                assembly_scope=selection.assembly_scope,
                path=selection.path or source_path,
                chain_locators=selection.chain_locators,
            )
        return selection
    model = next(source_model for source_model in structure if _model_id(source_model) == model_id)
    authors: list[str] = []
    for chain in model:
        author = str(chain.id)
        # A default selection describes analyzed polymer chains only. A chain
        # containing solely solvent/ions must not make ParsedStructure's chain
        # coherence contract fail, since those components are intentionally not
        # represented as ProteinChain values.
        if not any(
            _classify_residue(
                str(residue.get_resname()).strip().upper(),
                _residue_group(residue, model_id, atom_metadata),
                _all_atoms(residue),
            )
            is ComponentKind.POLYMER_RESIDUE
            for residue in chain
        ):
            continue
        authors.append(author)
    return InputSelection(
        content_id=snapshot.content_id,
        display_name=snapshot.display_name,
        format=fmt,
        model_id=model_id,
        author_chain_ids=tuple(authors),
        label_chain_ids=(),
        altloc_policy=AltlocPolicy.HIGHEST_OCCUPANCY,
        assembly_scope=AssemblyScope.ASYMMETRIC_UNIT,
        path=source_path,
    )


def _validate_selection(selection: InputSelection, snapshot: SourceSnapshot) -> None:
    if selection.content_id != snapshot.content_id:
        raise ValueError("selection content_id does not match captured source")
    if selection.format.value != snapshot.logical_format:
        raise ValueError("selection format does not match captured source")


def _decode_source(snapshot: SourceSnapshot) -> str:
    try:
        return snapshot.decompressed_bytes.decode("utf-8")
    except UnicodeDecodeError as error:
        raise StructureParseError(
            f"unable to decode structure source {snapshot.display_name}: source is not valid UTF-8 text"
        ) from error


def _structure_id(display_name: str, logical_format: str) -> str:
    name = display_name
    suffixes = (".gz", ".pdb", ".ent") if logical_format == "pdb" else (".gz", ".cif", ".mmcif")
    for suffix in suffixes:
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
    return name


def _model_id(model: Any) -> str:
    serial = getattr(model, "serial_num", None)
    if serial is None:
        serial = getattr(model, "id", "0")
    return str(serial)


def _all_atoms(residue: Any) -> list[Any]:
    atoms: list[Any] = []
    for atom in residue.child_dict.values():
        if atom.is_disordered() == 2:
            atoms.extend(atom.disordered_get_list())
        else:
            atoms.append(atom)
    return atoms


def _chain_count(structure: Any) -> int:
    return sum(1 for model in structure for _ in model)


def _residue_count(structure: Any) -> int:
    return sum(1 for model in structure for chain in model for _ in chain)


def _atom_count(structure: Any) -> int:
    return sum(len(_all_atoms(residue)) for model in structure for chain in model for residue in chain)


def _first_atom_metadata(
    atoms: Sequence[Any],
    model_id: str,
    atom_metadata: Mapping[tuple[str, str], _AtomSiteMetadata],
) -> _AtomSiteMetadata | None:
    for atom in atoms:
        serial = getattr(atom, "serial_number", None)
        if serial is not None:
            found = atom_metadata.get((model_id, str(serial)))
            if found is not None:
                return found
    return None


def _residue_identifiers(
    residue: Any,
    source_meta: _AtomSiteMetadata | None,
    label_numbering: Mapping[tuple[str, str, str, str | None, str], str | None] | None,
    model_id: str,
    chain_id: str,
) -> tuple[str, str | None, str | None, str | None, str]:
    _, raw_auth_seq, raw_icode = residue.id
    auth_seq_id = str(raw_auth_seq)
    insertion_code = _none_if_unknown(str(raw_icode))
    if source_meta is not None:
        auth_seq_id = source_meta.auth_seq_id
        insertion_code = source_meta.insertion_code
        label_seq_id = source_meta.label_seq_id
        entity_id = source_meta.entity_id
        group = source_meta.group
    else:
        label_seq_id = None
        entity_id = None
        group = "ATOM" if residue.id[0] == " " else "HETATM"
    if label_numbering is not None:
        label_seq_id = label_numbering.get(
            (model_id, chain_id, auth_seq_id, insertion_code, str(residue.get_resname()).upper())
        )
    return auth_seq_id, insertion_code, label_seq_id, entity_id, group


def _residue_group(
    residue: Any,
    model_id: str,
    atom_metadata: Mapping[tuple[str, str], _AtomSiteMetadata],
) -> str:
    metadata = _first_atom_metadata(_all_atoms(residue), model_id, atom_metadata)
    if metadata is not None:
        return metadata.group
    return "ATOM" if residue.id[0] == " " else "HETATM"


def _atom_record(
    atom: Any,
    model_id: str,
    atom_metadata: Mapping[tuple[str, str], _AtomSiteMetadata],
) -> AtomRecord:
    serial = getattr(atom, "serial_number", None)
    source_meta = atom_metadata.get((model_id, str(serial))) if serial is not None else None
    name = str(atom.get_name()).strip()
    element = str(getattr(atom, "element", "")).strip().upper() or name[0]
    source_serial: int | str | None
    if isinstance(serial, int):
        source_serial = serial
    elif serial is None:
        source_serial = None
    else:
        source_serial = str(serial)
    b_factor = atom.get_bfactor()
    occupancy = atom.get_occupancy()
    return AtomRecord(
        name=name,
        element=element,
        coordinate=tuple(float(value) for value in atom.get_coord()),
        altloc=_none_if_unknown(str(atom.get_altloc())),
        occupancy=float(occupancy) if occupancy is not None else None,
        b_factor=float(b_factor) if b_factor is not None else None,
        source_atom_id=source_meta.source_id
        if source_meta is not None
        else (str(serial) if serial is not None else None),
        source_serial=source_serial,
        formal_charge=source_meta.formal_charge if source_meta is not None else None,
    )


def _positive_float(value: object) -> float | None:
    try:
        parsed = float(str(value)) if value is not None else None
    except (TypeError, ValueError):
        return None
    return parsed if parsed is not None and parsed > 0.0 else None


def _none_if_unknown(value: str) -> str | None:
    return _clean_optional(value)


def _warning_diagnostics(warnings_seen: Sequence[Any], source_name: str) -> tuple[Diagnostic, ...]:
    return tuple(
        Diagnostic(
            code="biopython_warning",
            severity=DiagnosticSeverity.WARNING,
            message=str(item.message),
            source_id=source_name,
        )
        for item in warnings_seen
    )


__all__ = [
    "KNOWN_LIGANDS",
    "LIGAND_TABLE_VERSION",
    "ION_NAMES",
    "load_structure",
    "load_structure_legacy",
    "load_structure_evidence",
    "normalize_biopython_structure",
]
