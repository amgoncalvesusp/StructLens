"""Persistence boundary for canonical :class:`AnalysisReport` values.

Only report JSON is persisted as scientific identity.  Audit events are
workflow metadata and are deliberately ignored by ``serialize_report``.  A
coordinate path is never read unless the caller explicitly asks for its hash
to be verified against the report's content-addressed selection.
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any, cast

from structlens.application.report_pocket_serialization import pocket_report_from_dict
from structlens.application.report_snapshot_io import (
    MAX_SNAPSHOT_BYTES,
    load_snapshot,
    snapshot_path,
    store_snapshot,
    verify_report_inputs,
    verify_snapshot,
    verify_source_path,
)
from structlens.core.difference_maps import ResidueDisplacementVector
from structlens.core.evidence import (
    Availability,
    Diagnostic,
    DiagnosticSeverity,
    EvidenceCard,
    EvidenceQuality,
    InteractionEvidence,
    PocketEvidenceChannel,
    PocketEvidenceSections,
    SequenceEvidence,
    SiteEvidence,
    StructureEvidence,
)
from structlens.core.interactions import (
    InteractionChange,
    InteractionDifference,
    InteractionRecord,
    InteractionType,
    ReferenceInteractionKey,
)
from structlens.core.models import (
    ComparisonMode,
    CorrespondenceStatus,
    MutationEvent,
    MutationKind,
    ResidueId,
    StructuralTransform,
)
from structlens.core.msa import (
    AnalysisSequence,
    MSAColumn,
    MSAResidueCell,
    MultipleSequenceAlignment,
    SequenceResidueRef,
)
from structlens.core.parsing import (
    AltlocPolicy,
    AssemblyScope,
    ChainLocator,
    InputSelection,
    SourceSnapshot,
    StructureFormat,
)
from structlens.core.provenance import AuditEvent, MethodProvenance
from structlens.core.quality import CoordinateQCSettings, StructureQualityReport
from structlens.core.reports import (
    AnalysisReport,
    AnalysisSnapshot,
    CorrespondenceSnapshot,
    DistanceMapSnapshot,
    InputQualityBundle,
)
from structlens.core.reports.safe_io import atomic_write_bytes, read_bounded_bytes
from structlens.core.reports.schema import validate_analysis_report_payload
from structlens.core.sites import SiteDefinition, SiteDefinitionMode, SiteMetrics

_MAX_SNAPSHOT_BYTES = MAX_SNAPSHOT_BYTES


def report_to_dict(report: AnalysisReport) -> dict[str, Any]:
    if not isinstance(report, AnalysisReport):
        raise TypeError("report must be an AnalysisReport")
    payload = cast(dict[str, Any], report.to_json())
    validate_analysis_report_payload(payload)
    return payload


def serialize_report(
    report: AnalysisReport,
    *,
    snapshots: Sequence[SourceSnapshot] = (),
    source_paths: Sequence[str | Path] = (),
    snapshot_dir: str | Path | None = None,
    audit_event: AuditEvent | None = None,
) -> bytes:
    """Serialize one canonical report, optionally verifying source evidence.

    ``audit_event`` is accepted for API symmetry but never enters canonical
    bytes.  Verification is opt-in so a report can be exported from a project
    that has persisted snapshots without requiring the original paths.
    """

    if audit_event is not None and not isinstance(audit_event, AuditEvent):
        raise TypeError("audit_event must be an AuditEvent or None")
    payload = report_to_dict(report)
    if snapshot_dir is not None:
        # Validate identities and mutable paths before creating any persisted
        # artifacts; the second pass verifies the stored bytes and metadata.
        if snapshots or source_paths:
            verify_report_inputs(report, tuple(snapshots), tuple(source_paths), None)
        for snapshot in snapshots:
            store_snapshot(snapshot, snapshot_dir)
    verify_report_inputs(report, tuple(snapshots), tuple(source_paths), snapshot_dir)
    encoded = json.dumps(payload, ensure_ascii=False, allow_nan=False, sort_keys=True, separators=(",", ":"))
    return encoded.encode("utf-8")


def deserialize_report(
    payload: Mapping[str, Any] | str | bytes,
    *,
    snapshots: Sequence[SourceSnapshot] = (),
    source_paths: Sequence[str | Path] = (),
    snapshot_dir: str | Path | None = None,
) -> AnalysisReport:
    """Restore and verify a canonical report from bytes, text, or a mapping."""

    values = validate_analysis_report_payload(payload)
    report = _report_from_dict(values)
    verify_report_inputs(report, tuple(snapshots), tuple(source_paths), snapshot_dir)
    if values["report_id"] != report.report_id:
        raise ValueError("report_id does not match canonical report content")
    if _canonical_payload_bytes(values) != _canonical_payload_bytes(report_to_dict(report)):
        raise ValueError("analysis-report payload is not a canonical typed round-trip")
    return report


def _canonical_payload_bytes(payload: Mapping[str, Any]) -> bytes:
    return json.dumps(payload, ensure_ascii=False, allow_nan=False, sort_keys=True, separators=(",", ":")).encode(
        "utf-8"
    )


def save_report(
    report: AnalysisReport,
    path: str | Path,
    *,
    snapshots: Sequence[SourceSnapshot] = (),
    source_paths: Sequence[str | Path] = (),
    snapshot_dir: str | Path | None = None,
    audit_event: AuditEvent | None = None,
) -> None:
    output = Path(path)
    atomic_write_bytes(
        output,
        serialize_report(
            report,
            snapshots=snapshots,
            source_paths=source_paths,
            snapshot_dir=snapshot_dir,
            audit_event=audit_event,
        ),
        max_bytes=_MAX_SNAPSHOT_BYTES,
        label="report file",
        idempotent=True,
    )


def load_report(
    path: str | Path,
    *,
    snapshots: Sequence[SourceSnapshot] = (),
    source_paths: Sequence[str | Path] = (),
    snapshot_dir: str | Path | None = None,
) -> AnalysisReport:
    source = Path(path)
    data = read_bounded_bytes(source, max_bytes=_MAX_SNAPSHOT_BYTES, label="report file")
    return deserialize_report(data, snapshots=snapshots, source_paths=source_paths, snapshot_dir=snapshot_dir)


def _report_from_dict(payload: Mapping[str, Any]) -> AnalysisReport:
    return AnalysisReport(
        reference_selection=_selection_from_dict(_object(payload["reference_selection"], "reference_selection")),
        target_selection=_selection_from_dict(_object(payload["target_selection"], "target_selection")),
        input_quality=_input_quality_from_dict(_object(payload["input_quality"], "input_quality")),
        analysis=None if payload["analysis"] is None else _analysis_from_dict(_object(payload["analysis"], "analysis")),
        msa=None if payload["msa"] is None else _msa_from_dict(_object(payload["msa"], "msa")),
        interactions=_interactions_from_json(payload["interactions"]),
        sites=_sites_from_json(payload["sites"]),
        distance_map=None
        if payload["distance_map"] is None
        else _distance_map_from_dict(_object(payload["distance_map"], "distance_map")),
        displacement_vectors=tuple(
            _vector_from_dict(_object(item, "vector")) for item in _array(payload, "displacement_vectors")
        ),
        evidence_cards=tuple(
            _evidence_card_from_dict(_object(item, "evidence card")) for item in _array(payload, "evidence_cards")
        ),
        site_definitions=tuple(
            _site_definition_from_dict(_object(item, "site definition")) for item in _array(payload, "site_definitions")
        ),
        schema_version=str(payload["schema_version"]),
        comparison_mode=ComparisonMode(str(payload["comparison_mode"])),
        availability=_availability_from_dict(_object(payload["availability"], "availability")),
        diagnostics=tuple(
            _diagnostic_from_dict(_object(item, "diagnostic")) for item in _array(payload, "diagnostics")
        ),
        provenance=None
        if payload["provenance"] is None
        else MethodProvenance.from_json(_object(payload["provenance"], "provenance")),
        pockets=None
        if payload.get("pockets") is None
        else pocket_report_from_dict(_object(payload["pockets"], "pockets")),
    )


def _selection_from_dict(payload: Mapping[str, Any]) -> InputSelection:
    values = dict(payload)
    locators = tuple(ChainLocator(**dict(item)) for item in _array(values, "chain_locators"))
    selection = InputSelection(
        content_id=str(values["content_id"]),
        display_name=str(values.get("display_name", values["content_id"])),
        format=StructureFormat(str(values["format"])),
        model_id=str(values["model_id"]),
        author_chain_ids=tuple(str(item) for item in _array(values, "author_chain_ids")),
        label_chain_ids=tuple(str(item) for item in _array(values, "label_chain_ids")),
        altloc_policy=AltlocPolicy(str(values["altloc_policy"])),
        assembly_scope=AssemblyScope(str(values["assembly_scope"])),
        path=None,
        chain_locators=locators,
    )
    if values.get("selection_id") != selection.selection_id:
        raise ValueError("selection_id does not match canonical selection content")
    return selection


def _input_quality_from_dict(payload: Mapping[str, Any]) -> InputQualityBundle:
    return InputQualityBundle(
        _quality_from_dict(_object(payload["reference"], "reference")),
        _quality_from_dict(_object(payload["target"], "target")),
    )


def _quality_from_dict(payload: Mapping[str, Any]) -> StructureQualityReport:
    settings = dict(_object(payload["settings"], "settings"))
    return StructureQualityReport(
        Availability(str(payload["availability"])),
        diagnostics=tuple(
            _diagnostic_from_dict(_object(item, "diagnostic")) for item in _array(payload, "diagnostics")
        ),
        counts={str(key): int(value) for key, value in _object(payload["counts"], "counts").items()},
        settings=CoordinateQCSettings(**settings),
        provenance=None
        if payload["provenance"] is None
        else MethodProvenance.from_json(_object(payload["provenance"], "provenance")),
    )


def _analysis_from_dict(payload: Mapping[str, Any]) -> AnalysisSnapshot:
    transform = payload.get("transform")
    return AnalysisSnapshot(
        reference_id=str(payload["reference_id"]),
        target_id=str(payload["target_id"]),
        correspondences=tuple(
            _correspondence_from_dict(_object(item, "correspondence")) for item in _array(payload, "correspondences")
        ),
        mutations=tuple(_mutation_from_dict(_object(item, "mutation")) for item in _array(payload, "mutations")),
        sequence_identity=float(payload["sequence_identity"]),
        sequence_coverage=float(payload["sequence_coverage"]),
        alignment_decision=str(payload["alignment_decision"]),
        sequence_similarity=_optional_float(payload.get("sequence_similarity")),
        strict_rmsd_angstrom=_optional_float(payload.get("strict_rmsd_angstrom")),
        refined_rmsd_angstrom=_optional_float(payload.get("refined_rmsd_angstrom")),
        mapped_residue_count=int(payload.get("mapped_residue_count", 0)),
        refined_residue_count=_optional_int(payload.get("refined_residue_count")),
        excluded_alignment_indices=tuple(int(item) for item in _array(payload, "excluded_alignment_indices")),
        tm_score=_optional_float(payload.get("tm_score")),
        legacy_provenance={
            str(key): str(value) for key, value in _object(payload["legacy_provenance"], "legacy_provenance").items()
        },
        transform=None if transform is None else _transform_from_dict(_object(transform, "transform")),
        method_provenance=None
        if payload["method_provenance"] is None
        else MethodProvenance.from_json(_object(payload["method_provenance"], "method_provenance")),
    )


def _correspondence_from_dict(payload: Mapping[str, Any]) -> CorrespondenceSnapshot:
    return CorrespondenceSnapshot(
        alignment_index=int(payload["alignment_index"]),
        reference=_residue_from_json(payload.get("reference")),
        target=_residue_from_json(payload.get("target")),
        reference_one_letter=payload.get("reference_one_letter"),
        target_one_letter=payload.get("target_one_letter"),
        status=CorrespondenceStatus(str(payload["status"])),
        sequence_score=_optional_float(payload.get("sequence_score")),
        ca_displacement_angstrom=_optional_float(payload.get("ca_displacement_angstrom")),
        backbone_rmsd_angstrom=_optional_float(payload.get("backbone_rmsd_angstrom")),
        sidechain_rmsd_angstrom=_optional_float(payload.get("sidechain_rmsd_angstrom")),
        all_heavy_atom_rmsd_angstrom=_optional_float(payload.get("all_heavy_atom_rmsd_angstrom")),
        is_outlier=bool(payload.get("is_outlier", False)),
        is_key_residue=bool(payload.get("is_key_residue", False)),
        mapping_source=str(payload.get("mapping_source", "unknown")),
        mapping_locked=bool(payload.get("mapping_locked", False)),
    )


def _mutation_from_dict(payload: Mapping[str, Any]) -> MutationEvent:
    return MutationEvent(
        int(payload["alignment_index"]),
        MutationKind(str(payload["kind"])),
        _residue_from_json(payload.get("reference")),
        _residue_from_json(payload.get("target")),
        None if payload.get("reference_aa") is None else str(payload["reference_aa"]),
        None if payload.get("target_aa") is None else str(payload["target_aa"]),
        str(payload["reference_label"]),
        str(payload["target_label"]),
        str(payload["canonical_notation"]),
        _optional_int(payload.get("blosum62_score")),
        _optional_int(payload.get("grantham_distance")),
        None if payload.get("physicochemical_class") is None else str(payload["physicochemical_class"]),
    )


def _msa_from_dict(payload: Mapping[str, Any]) -> MultipleSequenceAlignment:
    def ref(value: object) -> SequenceResidueRef:
        item = _object(value, "sequence residue")
        return SequenceResidueRef(
            int(item["sequence_index"]), str(item["one_letter"]), _residue_from_json(item.get("residue_id"))
        )

    sequences = tuple(
        AnalysisSequence(
            str(row["structure_id"]),
            row.get("chain_id"),
            str(row["sequence"]),
            tuple(ref(item) for item in _array(row, "residues")),
            cast(Any, str(row["source"])),
        )
        for row in (_object(item, "MSA sequence") for item in _array(payload, "sequences"))
    )
    columns = []
    for raw in _array(payload, "columns"):
        column = _object(raw, "MSA column")
        cells = []
        for cell_raw in _array(column, "cells"):
            cell = _object(cell_raw, "MSA cell")
            cells.append(
                MSAResidueCell(
                    str(cell["structure_id"]),
                    int(cell["alignment_column"]),
                    None if cell["residue"] is None else ref(cell["residue"]),
                    str(cell["character"]),
                )
            )
        columns.append(
            MSAColumn(
                int(column["index"]),
                str(column["reference_label"]),
                _residue_from_json(column.get("reference_residue")),
                tuple(cells),
                int(column["non_gap_count"]),
                float(column["gap_fraction"]),
                float(column["ambiguous_fraction"]),
                _optional_float(column.get("conservation_score")),
                _optional_float(column.get("entropy_bits")),
            )
        )
    return MultipleSequenceAlignment(
        sequences,
        tuple(
            (str(item[0]), str(item[1]))
            for item in (_array_value(row, "aligned row") for row in _array(payload, "aligned_rows"))
        ),
        tuple(columns),
        payload.get("reference_structure_id"),
        str(payload["algorithm"]),
        tuple(str(item) for item in _array(payload, "provenance")),
    )


def _interactions_from_json(
    value: object,
) -> InteractionEvidence | tuple[InteractionDifference | InteractionRecord, ...] | None:
    if value is None:
        return None
    if isinstance(value, Mapping) and "differences" in value:
        return InteractionEvidence(
            tuple(
                _interaction_difference_from_dict(_object(item, "interaction difference"))
                for item in _array(value, "differences")
            ),
            tuple(
                _interaction_record_from_dict(_object(item, "interaction record"))
                for item in _array(value, "reference_interactions")
            ),
            tuple(
                _interaction_record_from_dict(_object(item, "interaction record"))
                for item in _array(value, "target_interactions")
            ),
        )
    return tuple(
        _interaction_difference_from_dict(_object(item, "interaction difference"))
        if isinstance(item, Mapping) and "key" in item
        else _interaction_record_from_dict(_object(item, "interaction record"))
        for item in _array_value(value, "interactions")
    )


def _interaction_record_from_dict(payload: Mapping[str, Any]) -> InteractionRecord:
    return InteractionRecord(
        str(payload["structure_id"]),
        InteractionType(str(payload["interaction_type"])),
        _required_residue(payload.get("residue_a"), "residue_a"),
        _residue_from_json(payload.get("residue_b")),
        payload.get("atom_a"),
        payload.get("atom_b"),
        float(payload["distance_angstrom"]),
        _optional_float(payload.get("angle_degrees")),
        payload.get("ligand_or_metal_id"),
        str(payload.get("evidence_mode", "heavy_atom_geometry")),
    )


def _interaction_difference_from_dict(payload: Mapping[str, Any]) -> InteractionDifference:
    key = _object(payload["key"], "interaction key")
    return InteractionDifference(
        ReferenceInteractionKey(
            InteractionType(str(key["interaction_type"])),
            str(key["reference_position_a"]),
            key.get("reference_position_b"),
            key.get("external_partner_id"),
        ),
        InteractionChange(str(payload["change"])),
        None
        if payload.get("reference_record") is None
        else _interaction_record_from_dict(_object(payload["reference_record"], "reference record")),
        None
        if payload.get("target_record") is None
        else _interaction_record_from_dict(_object(payload["target_record"], "target record")),
    )


def _sites_from_json(value: object) -> SiteEvidence | tuple[SiteMetrics, ...] | None:
    if value is None:
        return None
    values = _array(value, "metrics") if isinstance(value, Mapping) else _array_value(value, "sites")
    metrics = tuple(SiteMetrics(**_object(item, "site metric")) for item in values)
    return SiteEvidence(metrics) if isinstance(value, Mapping) else metrics


def _distance_map_from_dict(payload: Mapping[str, Any]) -> DistanceMapSnapshot:
    return DistanceMapSnapshot(
        tuple(str(item) for item in _array(payload, "reference_positions")),
        tuple(
            tuple(float(item) for item in _array_value(row, "distance row"))
            for row in _array(payload, "reference_distances_angstrom")
        ),
        tuple(
            tuple(float(item) for item in _array_value(row, "distance row"))
            for row in _array(payload, "target_distances_angstrom")
        ),
        tuple(
            tuple(float(item) for item in _array_value(row, "distance row"))
            for row in _array(payload, "delta_angstrom")
        ),
        tuple(tuple(bool(item) for item in _array_value(row, "mask row")) for row in _array(payload, "valid_mask")),
    )


def _vector_from_dict(payload: Mapping[str, Any]) -> ResidueDisplacementVector:
    return ResidueDisplacementVector(
        str(payload["reference_position"]),
        _required_residue(payload["reference_residue"], "reference_residue"),
        _required_residue(payload["target_residue"], "target_residue"),
        _coordinate(payload, "start_xyz"),
        _coordinate(payload, "end_xyz"),
        _coordinate(payload, "vector_xyz"),
        float(payload["magnitude_angstrom"]),
    )


def _site_definition_from_dict(payload: Mapping[str, Any]) -> SiteDefinition:
    return SiteDefinition(
        str(payload["site_id"]),
        str(payload["name"]),
        SiteDefinitionMode(str(payload["mode"])),
        tuple(_required_residue(item, "reference_residues item") for item in _array(payload, "reference_residues")),
        _residue_from_json(payload.get("center_residue")),
        payload.get("ligand_id"),
        _optional_float(payload.get("radius_angstrom")),
    )


def _evidence_card_from_dict(payload: Mapping[str, Any]) -> EvidenceCard:
    residue_ref_payload = _object(payload["residue_ref"], "residue_ref")
    residue_ref = SequenceResidueRef(
        int(residue_ref_payload["sequence_index"]),
        str(residue_ref_payload["one_letter"]),
        _residue_from_json(residue_ref_payload.get("residue_id")),
    )
    sequence = _object(payload["sequence"], "sequence")
    structure = _object(payload["structure"], "structure")
    quality = _object(payload["quality"], "quality")
    return EvidenceCard(
        residue_ref,
        payload.get("target_id"),
        SequenceEvidence(
            str(sequence["reference_one_letter"]) if sequence.get("reference_one_letter") is not None else None,
            str(sequence["target_one_letter"]) if sequence.get("target_one_letter") is not None else None,
            _optional_int(sequence.get("alignment_index")),
            _optional_float(sequence.get("sequence_identity")),
            _optional_float(sequence.get("conservation_fraction")),
            _optional_float(sequence.get("entropy_bits")),
            _optional_float(sequence.get("gap_fraction")),
            _optional_float(sequence.get("ambiguous_fraction")),
            tuple(
                SequenceResidueRef(
                    int(item["sequence_index"]), str(item["one_letter"]), _residue_from_json(item.get("residue_id"))
                )
                for item in _array(sequence, "source_refs")
            ),
        ),
        StructureEvidence(**{str(key): value for key, value in structure.items()}),
        _required_interactions(payload["interactions"]),
        _required_sites(payload["site"]),
        EvidenceQuality(
            str(quality["overall_status"]),
            tuple(str(item) for item in _array(quality, "available_sections")),
            tuple(str(item) for item in _array(quality, "unavailable_sections")),
            tuple(str(item) for item in _array(quality, "warnings")),
            _optional_float(quality.get("coverage_fraction")),
            _optional_int(quality.get("source_count")),
        ),
        str(payload["schema_version"]),
        tuple(str(item) for item in _array(payload, "provenance")),
        _pocket_sections_from_dict(payload.get("pocket_concordance")),
    )


def _pocket_sections_from_dict(value: object) -> PocketEvidenceSections | None:
    if value is None:
        return None
    payload = _object(value, "pocket_concordance")
    expected = {
        "geometry",
        "ligand_support",
        "parameter_persistence",
        "volume_sensitivity",
        "match_ambiguity",
        "qc",
        "interactions",
    }
    unknown = set(payload).difference(expected)
    if unknown:
        raise ValueError(f"pocket_concordance contains unknown channels: {', '.join(sorted(unknown))}")
    if set(payload) != expected:
        raise ValueError("pocket_concordance must contain all seven channels")
    channels: dict[str, PocketEvidenceChannel] = {}
    for name in sorted(expected):
        channel = _object(payload[name], f"pocket_concordance.{name}")
        channels[name] = PocketEvidenceChannel(
            Availability(str(channel["availability"])),
            measure=channel.get("measure"),
            units=dict(_object(channel.get("units", {}), f"{name}.units")),
            diagnostics=tuple(
                _diagnostic_from_dict(_object(item, f"{name} diagnostic")) for item in _array(channel, "diagnostics")
            ),
            provenance=channel.get("provenance"),
        )
    return PocketEvidenceSections(**channels)


def _availability_from_dict(payload: Mapping[str, Any]) -> Any:
    from structlens.core.reports import SectionAvailability

    return SectionAvailability(**{str(key): Availability(str(value)) for key, value in payload.items()})


def _diagnostic_from_dict(payload: Mapping[str, Any]) -> Diagnostic:
    return Diagnostic(
        str(payload["code"]),
        DiagnosticSeverity(str(payload["severity"])),
        str(payload["message"]),
        payload.get("source_id"),
        payload.get("atom_id"),
        payload.get("residue_id"),
        payload.get("remediation"),
    )


def _transform_from_dict(payload: Mapping[str, Any]) -> StructuralTransform:
    return StructuralTransform(
        tuple(_coordinate(row, "rotation row") for row in _array(payload, "rotation")),
        _coordinate(payload, "translation"),
    )


def _residue_from_json(value: object) -> ResidueId | None:
    if value is None:
        return None
    return ResidueId(**_object(value, "residue"))


def _required_residue(value: object, name: str) -> ResidueId:
    residue = _residue_from_json(value)
    if residue is None:
        raise ValueError(f"{name} must contain a residue")
    return residue


def _coordinate(payload: object, name: str) -> tuple[float, float, float]:
    raw: object = payload.get(name) if isinstance(payload, Mapping) else payload
    values = tuple(float(item) for item in _array_value(raw, name))
    if len(values) != 3:
        raise ValueError(f"{name} must contain exactly three values")
    return values[0], values[1], values[2]


def _required_interactions(value: object) -> InteractionEvidence:
    interactions = _interactions_from_json(value)
    if not isinstance(interactions, InteractionEvidence):
        raise ValueError("evidence card interactions must be an InteractionEvidence object")
    return interactions


def _required_sites(value: object) -> SiteEvidence:
    sites = _sites_from_json(value)
    if not isinstance(sites, SiteEvidence):
        raise ValueError("evidence card site must be a SiteEvidence object")
    return sites


def _object(value: object, name: str) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be an object")
    return {str(key): item for key, item in value.items()}


def _array(value: Mapping[str, Any], name: str) -> list[Any]:
    return _array_value(value.get(name), name)


def _array_value(value: object, name: str) -> list[Any]:
    if not isinstance(value, (list, tuple)):
        raise ValueError(f"{name} must be an array")
    if len(value) > 100_000:
        raise ValueError(f"{name} exceeds the maximum item limit")
    return list(value)


def _optional_float(value: object) -> float | None:
    return None if value is None else float(cast(float, value))


def _optional_int(value: object) -> int | None:
    return None if value is None else int(cast(int, value))


# Descriptive aliases used by application and CLI callers during migration.
serialize_analysis_report = serialize_report
deserialize_analysis_report = deserialize_report
save_analysis_report = save_report
load_analysis_report = load_report
persist_snapshot = store_snapshot
canonical_report_bytes = serialize_report
serialize_report_json = serialize_report
deserialize_report_json = deserialize_report


__all__ = [
    "deserialize_analysis_report",
    "deserialize_report",
    "load_analysis_report",
    "load_report",
    "load_snapshot",
    "persist_snapshot",
    "report_to_dict",
    "save_analysis_report",
    "save_report",
    "serialize_analysis_report",
    "serialize_report",
    "serialize_report_json",
    "deserialize_report_json",
    "canonical_report_bytes",
    "snapshot_path",
    "store_snapshot",
    "verify_snapshot",
    "verify_source_path",
]
