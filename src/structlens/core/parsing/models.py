"""Immutable parser-boundary contracts for selections and retained evidence."""

from __future__ import annotations

import hashlib
import json
import math
import re
from dataclasses import dataclass
from enum import Enum
from typing import TYPE_CHECKING

from structlens.core.evidence.status import Diagnostic
from structlens.core.models.components import StructureComponent
from structlens.core.models.structure import ProteinStructure

if TYPE_CHECKING:
    from collections.abc import Sequence


class StructureFormat(str, Enum):
    """Coordinate formats supported by the parser boundary."""

    PDB = "pdb"
    MMCIF = "mmcif"
    CIF = "mmcif"


class AltlocPolicy(str, Enum):
    """Policy for the primary representation of alternate locations."""

    HIGHEST_OCCUPANCY = "highest_occupancy"
    FIRST = "first"
    ALL = "all"


class AssemblyScope(str, Enum):
    """Coordinate assembly represented by a parsed input."""

    ASYMMETRIC_UNIT = "asymmetric_unit"
    BIOLOGICAL_ASSEMBLY = "biological_assembly"


_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_MISSING_HASH = object()


def _enum(value: object, enum_type: type[Enum], name: str) -> Enum:
    if isinstance(value, enum_type):
        return value
    try:
        return enum_type(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"unknown {name}: {value!r}") from exc


def _text(value: object, name: str) -> str:
    if value is None:
        raise ValueError(f"{name} must not be None")
    if not isinstance(value, str):
        value = str(value)
    value = value.strip()
    if not value:
        raise ValueError(f"{name} must not be empty")
    return value


def _ids(values: Sequence[object], name: str) -> tuple[str, ...]:
    normalized = tuple(_text(value, name) for value in values)
    if len(set(normalized)) != len(normalized):
        raise ValueError(f"{name} must not contain duplicates")
    return tuple(sorted(normalized))


def _chain_ids(values: Sequence[object], name: str) -> tuple[str, ...]:
    normalized: list[str] = []
    for value in values:
        if value is None:
            raise ValueError(f"{name} must not be None")
        normalized.append(str(value).strip())
    if len(set(normalized)) != len(normalized):
        raise ValueError(f"{name} must not contain duplicates")
    return tuple(sorted(normalized))


def _optional_text(value: object, name: str) -> str | None:
    if value is None:
        return None
    return _text(value, name)


def _optional_chain_id(value: object, name: str) -> str | None:
    if value is None:
        return None
    return str(value).strip()


@dataclass(frozen=True, slots=True)
class ChainLocator:
    """A frozen author/label/entity chain correspondence."""

    author_chain_id: str | None = None
    label_chain_id: str | None = None
    entity_id: str | None = None

    def __post_init__(self) -> None:
        author = _optional_chain_id(self.author_chain_id, "author_chain_id")
        label = _optional_chain_id(self.label_chain_id, "label_chain_id")
        entity = _optional_text(self.entity_id, "entity_id")
        if author is None and label is None and entity is None:
            raise ValueError("chain locator must identify an author, label, or entity chain")
        object.__setattr__(self, "author_chain_id", author)
        object.__setattr__(self, "label_chain_id", label)
        object.__setattr__(self, "entity_id", entity)


@dataclass(frozen=True, slots=True, init=False)
class InputSelection:
    """User-selected, content-addressed parser input.

    ``path`` and ``display_name`` are presentation metadata.  They are not
    included in :attr:`selection_id`, so replacing or renaming a path cannot
    change scientific identity.
    """

    content_id: str
    display_name: str
    format: StructureFormat
    model_id: str
    author_chain_ids: tuple[str, ...]
    label_chain_ids: tuple[str, ...]
    altloc_policy: AltlocPolicy
    assembly_scope: AssemblyScope
    path: str | None
    chain_locators: tuple[ChainLocator, ...]

    def __init__(
        self,
        content_id: str,
        display_name: str,
        format: StructureFormat | str | None = None,
        model_id: str | int | None = None,
        author_chain_ids: Sequence[object] = (),
        label_chain_ids: Sequence[object] = (),
        altloc_policy: AltlocPolicy | str = AltlocPolicy.HIGHEST_OCCUPANCY,
        assembly_scope: AssemblyScope | str = AssemblyScope.ASYMMETRIC_UNIT,
        path: str | None = None,
        chain_locators: Sequence[ChainLocator] = (),
        *,
        source_format: StructureFormat | str | None = None,
        selected_model: str | int | None = None,
        author_chain_id: str | None = None,
        label_chain_id: str | None = None,
    ) -> None:
        normalized_content_id = _text(content_id, "content_id").lower()
        if not _SHA256_RE.fullmatch(normalized_content_id):
            raise ValueError("content_id must be a 64-character SHA-256 hexadecimal digest")
        normalized_display_name = _text(display_name, "display_name")
        if format is None:
            format = source_format
        elif source_format is not None and _enum(format, StructureFormat, "structure format") != _enum(
            source_format, StructureFormat, "source format"
        ):
            raise ValueError("format and source_format must match")
        if format is None:
            raise ValueError("format is required")
        normalized_format = _enum(format, StructureFormat, "structure format")
        if model_id is None:
            model_id = selected_model
        elif selected_model is not None and str(model_id).strip() != str(selected_model).strip():
            raise ValueError("model_id and selected_model must match")
        normalized_model = _text(model_id, "model_id")
        author_values = tuple(author_chain_ids)
        if author_chain_id is not None:
            if author_values and author_values != (author_chain_id,):
                raise ValueError("author_chain_id and author_chain_ids must match")
            author_values = (author_chain_id,)
        label_values = tuple(label_chain_ids)
        if label_chain_id is not None:
            if label_values and label_values != (label_chain_id,):
                raise ValueError("label_chain_id and label_chain_ids must match")
            label_values = (label_chain_id,)
        normalized_author_ids = _chain_ids(author_values, "author_chain_id")
        normalized_label_ids = _chain_ids(label_values, "label_chain_id")
        locators = tuple(chain_locators)
        if any(not isinstance(locator, ChainLocator) for locator in locators):
            raise TypeError("chain_locators must contain ChainLocator values")
        if len(normalized_author_ids) > 1 and len(normalized_label_ids) > 1 and not locators:
            raise ValueError("multi-chain author/label selection is ambiguous; provide chain_locators")
        if locators:
            if len(set(locators)) != len(locators):
                raise ValueError("chain_locators must not contain duplicates")
            locators = tuple(
                sorted(
                    locators,
                    key=lambda locator: (
                        locator.author_chain_id or "",
                        locator.label_chain_id or "",
                        locator.entity_id or "",
                    ),
                )
            )
            locator_authors = tuple(sorted({item.author_chain_id for item in locators if item.author_chain_id is not None}))
            locator_labels = tuple(sorted({item.label_chain_id for item in locators if item.label_chain_id is not None}))
            if normalized_author_ids and normalized_author_ids != locator_authors:
                raise ValueError("author_chain_ids must match chain_locators")
            if normalized_label_ids and normalized_label_ids != locator_labels:
                raise ValueError("label_chain_ids must match chain_locators")
            normalized_author_ids = locator_authors
            normalized_label_ids = locator_labels
        else:
            locator_count = max(len(normalized_author_ids), len(normalized_label_ids))
            locators = tuple(
                ChainLocator(
                    normalized_author_ids[index] if index < len(normalized_author_ids) else None,
                    normalized_label_ids[index] if index < len(normalized_label_ids) else None,
                )
                for index in range(locator_count)
            )
        object.__setattr__(self, "content_id", normalized_content_id)
        object.__setattr__(self, "display_name", normalized_display_name)
        object.__setattr__(self, "format", normalized_format)
        object.__setattr__(self, "model_id", normalized_model)
        object.__setattr__(self, "author_chain_ids", normalized_author_ids)
        object.__setattr__(self, "label_chain_ids", normalized_label_ids)
        object.__setattr__(self, "altloc_policy", _enum(altloc_policy, AltlocPolicy, "altloc policy"))
        object.__setattr__(self, "assembly_scope", _enum(assembly_scope, AssemblyScope, "assembly scope"))
        if path is not None and not isinstance(path, str):
            path = str(path)
        object.__setattr__(self, "path", path.strip() if path is not None else None)
        object.__setattr__(self, "chain_locators", locators)

    @property
    def selected_model(self) -> str:
        """Alias used by metadata and parser adapters."""

        return self.model_id

    @property
    def source_format(self) -> StructureFormat:
        return self.format

    @property
    def author_chain_id(self) -> str | None:
        return self.author_chain_ids[0] if len(self.author_chain_ids) == 1 else None

    @property
    def label_chain_id(self) -> str | None:
        return self.label_chain_ids[0] if len(self.label_chain_ids) == 1 else None

    @property
    def selection_id(self) -> str:
        """Deterministic identifier for content and scientific settings."""

        canonical = {
            "assembly_scope": self.assembly_scope.value,
            "author_chain_ids": self.author_chain_ids,
            "content_id": self.content_id,
            "format": self.format.value,
            "label_chain_ids": self.label_chain_ids,
            "model_id": self.model_id,
            "altloc_policy": self.altloc_policy.value,
            "chain_locators": tuple(
                {
                    "author_chain_id": locator.author_chain_id,
                    "entity_id": locator.entity_id,
                    "label_chain_id": locator.label_chain_id,
                }
                for locator in self.chain_locators
            ),
        }
        serialized = json.dumps(canonical, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(serialized.encode("utf-8")).hexdigest()

    @property
    def identifier(self) -> str:
        return self.selection_id


@dataclass(frozen=True, slots=True, init=False)
class StructureMetadata:
    """Validated format/model/chain metadata retained from a source."""

    format: StructureFormat
    selected_model: str
    available_models: tuple[str, ...]
    analyzed_chain_ids: tuple[str, ...]
    experimental_method: str | None
    resolution_angstrom: float | None
    assembly_scope: AssemblyScope
    author_chain_ids: tuple[str, ...]
    label_chain_ids: tuple[str, ...]
    entity_ids: tuple[str, ...]

    def __init__(
        self,
        format: StructureFormat | str | None = None,
        selected_model: str | int | None = None,
        available_models: Sequence[object] = (),
        analyzed_chain_ids: Sequence[object] = (),
        experimental_method: str | None = None,
        resolution_angstrom: float | None = None,
        assembly_scope: AssemblyScope | str = AssemblyScope.ASYMMETRIC_UNIT,
        author_chain_ids: Sequence[object] = (),
        label_chain_ids: Sequence[object] = (),
        entity_ids: Sequence[object] = (),
        *,
        source_format: StructureFormat | str | None = None,
        model_id: str | int | None = None,
        analyzed_chains: Sequence[object] | None = None,
    ) -> None:
        if format is None:
            format = source_format
        elif source_format is not None and _enum(format, StructureFormat, "structure format") != _enum(
            source_format, StructureFormat, "source format"
        ):
            raise ValueError("format and source_format must match")
        if format is None:
            raise ValueError("format is required")
        if selected_model is None:
            selected_model = model_id
        elif model_id is not None and str(selected_model).strip() != str(model_id).strip():
            raise ValueError("selected_model and model_id must match")
        model = _text(selected_model, "selected_model")
        models = _ids(available_models or (model,), "available_model")
        if model not in models:
            raise ValueError("selected_model must be included in available_models")
        chains = analyzed_chain_ids if analyzed_chains is None else analyzed_chains
        resolution = None if resolution_angstrom is None else float(resolution_angstrom)
        if resolution is not None and (not math.isfinite(resolution) or resolution <= 0.0):
            raise ValueError("resolution_angstrom must be finite and positive")
        method = None if experimental_method is None else _text(experimental_method, "experimental_method")
        object.__setattr__(self, "format", _enum(format, StructureFormat, "structure format"))
        object.__setattr__(self, "selected_model", model)
        object.__setattr__(self, "available_models", models)
        object.__setattr__(self, "analyzed_chain_ids", _chain_ids(chains, "analyzed_chain_id"))
        object.__setattr__(self, "experimental_method", method)
        object.__setattr__(self, "resolution_angstrom", resolution)
        object.__setattr__(self, "assembly_scope", _enum(assembly_scope, AssemblyScope, "assembly scope"))
        object.__setattr__(self, "author_chain_ids", _chain_ids(author_chain_ids, "author_chain_id"))
        object.__setattr__(self, "label_chain_ids", _chain_ids(label_chain_ids, "label_chain_id"))
        object.__setattr__(self, "entity_ids", _ids(entity_ids, "entity_id"))

    @property
    def source_format(self) -> StructureFormat:
        return self.format

    @property
    def analyzed_chains(self) -> tuple[str, ...]:
        return self.analyzed_chain_ids


@dataclass(frozen=True, slots=True, init=False)
class ParsedStructure:
    """Normalized protein structure plus lossless parser evidence."""

    protein_structure: ProteinStructure
    components: tuple[StructureComponent, ...]
    selection: InputSelection
    metadata: StructureMetadata
    raw_source_hash: str
    diagnostics: tuple[Diagnostic, ...]

    def __init__(
        self,
        protein_structure: ProteinStructure | None = None,
        components: Sequence[StructureComponent] = (),
        selection: InputSelection | None = None,
        metadata: StructureMetadata | None = None,
        raw_source_hash: str | object = _MISSING_HASH,
        diagnostics: Sequence[Diagnostic] = (),
        *,
        structure: ProteinStructure | None = None,
        source_hash: str | None = None,
        raw_hash: str | None = None,
    ) -> None:
        if protein_structure is None:
            protein_structure = structure
        elif structure is not None and structure is not protein_structure:
            raise ValueError("structure and protein_structure must refer to one value")
        if protein_structure is None or not isinstance(protein_structure, ProteinStructure):
            raise TypeError("protein_structure must be a ProteinStructure")
        if selection is None:
            raise TypeError("selection is required")
        if metadata is None:
            raise TypeError("metadata is required")
        if not isinstance(selection, InputSelection):
            raise TypeError("selection must be an InputSelection")
        if not isinstance(metadata, StructureMetadata):
            raise TypeError("metadata must be a StructureMetadata")
        if metadata.format is not selection.format:
            raise ValueError("metadata format must match selection format")
        if metadata.selected_model != selection.model_id:
            raise ValueError("metadata selected model must match selection model")
        if metadata.assembly_scope is not selection.assembly_scope:
            raise ValueError("metadata assembly scope must match selection assembly scope")
        structure_chains = tuple(protein_structure.chains)
        if not structure_chains:
            if selection.author_chain_ids or selection.label_chain_ids or selection.chain_locators:
                raise ValueError("selected chain is not present in an empty structure")
            if metadata.analyzed_chain_ids:
                raise ValueError("metadata analyzed chains are not present in the structure")
        else:
            structure_models = {str(chain.model_id) for chain in structure_chains}
            if selection.model_id not in structure_models:
                raise ValueError("selection model is not present in the structure")
            selected_chains = tuple(chain for chain in structure_chains if str(chain.model_id) == selection.model_id)

            def identifiers(namespace: str) -> set[str]:
                values: set[str] = set()
                for chain in selected_chains:
                    value = getattr(chain, namespace)
                    if value is not None:
                        values.add(str(value))
                    if namespace == "author_chain_id" and value is None:
                        values.add(str(chain.chain_id))
                return values

            author_identifiers = identifiers("author_chain_id")
            label_identifiers = identifiers("label_chain_id")
            for chain_id in selection.author_chain_ids:
                if chain_id not in author_identifiers:
                    raise ValueError(f"selected author chain is not present: {chain_id}")
            for chain_id in selection.label_chain_ids:
                if chain_id not in label_identifiers:
                    raise ValueError(f"selected label chain is not present: {chain_id}")
            for locator in selection.chain_locators:
                matches = selected_chains
                if locator.author_chain_id is not None:
                    matches = tuple(
                        chain
                        for chain in matches
                        if (chain.author_chain_id if chain.author_chain_id is not None else chain.chain_id)
                        == locator.author_chain_id
                    )
                if locator.label_chain_id is not None:
                    matches = tuple(
                        chain
                        for chain in matches
                        if chain.label_chain_id == locator.label_chain_id
                    )
                if locator.entity_id is not None:
                    matches = tuple(chain for chain in matches if chain.entity_id == locator.entity_id)
                if not matches:
                    raise ValueError("selected chain locator is not present in the structure")
            expected_analyzed = {str(chain.chain_id) for chain in selected_chains}
            expected_author = set(author_identifiers)
            expected_label = set(label_identifiers)
            structure_entities = {
                str(chain.entity_id) for chain in selected_chains if chain.entity_id is not None
            }
            expected_entities = structure_entities
            if set(metadata.analyzed_chain_ids) != expected_analyzed:
                raise ValueError("metadata analyzed chains are not coherent with the structure")
            if set(metadata.author_chain_ids) != expected_author:
                raise ValueError("metadata author chains are not coherent with the structure")
            if set(metadata.label_chain_ids) != expected_label:
                raise ValueError("metadata label chains are not coherent with the structure")
            if set(metadata.entity_ids) != expected_entities:
                raise ValueError("metadata entities are not coherent with the structure")
            if metadata.author_chain_ids and not set(selection.author_chain_ids).issubset(
                metadata.author_chain_ids
            ):
                raise ValueError("selection author chains are not coherent with metadata")
            if metadata.label_chain_ids and not set(selection.label_chain_ids).issubset(
                metadata.label_chain_ids
            ):
                raise ValueError("selection label chains are not coherent with metadata")
        normalized_components = tuple(components)
        if any(not isinstance(item, StructureComponent) for item in normalized_components):
            raise TypeError("components must contain StructureComponent values")
        component_ids = tuple(item.component_id for item in normalized_components)
        if len(set(component_ids)) != len(component_ids):
            raise ValueError("components must have unique identities")
        normalized_hash = raw_source_hash
        if normalized_hash is _MISSING_HASH:
            if source_hash is None and raw_hash is None:
                raise TypeError("raw_source_hash is required")
            if source_hash is not None and raw_hash is not None and source_hash != raw_hash:
                raise ValueError("source_hash and raw_hash aliases must match")
            normalized_hash = source_hash if source_hash is not None else raw_hash
        elif not isinstance(normalized_hash, str):
            raise TypeError("raw_source_hash must be a string")
        else:
            if source_hash is not None and source_hash != normalized_hash:
                raise ValueError("source_hash and raw hash aliases must match")
            if raw_hash is not None and raw_hash != normalized_hash:
                raise ValueError("raw_hash and raw hash aliases must match")
            if source_hash is not None and raw_hash is not None and source_hash != raw_hash:
                raise ValueError("source_hash and raw hash aliases must match")
        if normalized_hash is None:
            raise TypeError("raw_source_hash is required")
        normalized_hash = _text(normalized_hash, "raw_source_hash").lower()
        if not _SHA256_RE.fullmatch(normalized_hash):
            raise ValueError("raw_source_hash must be a 64-character SHA-256 hexadecimal digest")
        normalized_diagnostics = tuple(diagnostics)
        if any(not isinstance(item, Diagnostic) for item in normalized_diagnostics):
            raise TypeError("diagnostics must contain Diagnostic values")
        object.__setattr__(self, "protein_structure", protein_structure)
        object.__setattr__(self, "components", normalized_components)
        object.__setattr__(self, "selection", selection)
        object.__setattr__(self, "metadata", metadata)
        object.__setattr__(self, "raw_source_hash", normalized_hash)
        object.__setattr__(self, "diagnostics", normalized_diagnostics)

    @property
    def structure(self) -> ProteinStructure:
        return self.protein_structure

    @property
    def source_hash(self) -> str:
        return self.raw_source_hash

    @property
    def content_id(self) -> str:
        """SHA-256 identity of the logical content selected for parsing."""

        return self.selection.content_id


__all__ = [
    "AltlocPolicy",
    "AssemblyScope",
    "ChainLocator",
    "InputSelection",
    "ParsedStructure",
    "StructureFormat",
    "StructureMetadata",
]
