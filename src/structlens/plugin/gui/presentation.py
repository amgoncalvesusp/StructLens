"""Immutable, report-first projections for the desktop GUI.

The Qt layer is deliberately not a scientific layer.  This module turns the
typed :class:`~structlens.core.reports.AnalysisReport` into small display
records, preserving native availability, diagnostics, units and provenance.
No metric is recalculated here; every number is copied from the report.
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from dataclasses import dataclass, field
from enum import Enum
from types import MappingProxyType
from typing import Any

from structlens.core.evidence import Availability, InteractionEvidence, SiteEvidence
from structlens.core.reports import AnalysisReport


class FrozenMap(Mapping[str, Any]):
    """Small immutable mapping that remains friendly to ``dataclasses.asdict``."""

    __slots__ = ("_values",)

    def __init__(self, values: Mapping[str, Any] | None = None) -> None:
        self._values = MappingProxyType({str(key): freeze_json(value) for key, value in (values or {}).items()})

    def __getitem__(self, key: str) -> Any:
        return self._values[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self._values)

    def __len__(self) -> int:
        return len(self._values)

    def __repr__(self) -> str:
        return repr(dict(self._values))

    def __deepcopy__(self, memo: dict[int, Any]) -> FrozenMap:
        del memo
        return self


def freeze_json(value: Any) -> Any:
    """Return a recursively immutable copy of JSON-compatible data."""

    if isinstance(value, FrozenMap):
        return value
    if isinstance(value, Mapping):
        return FrozenMap(value)
    if isinstance(value, (list, tuple)):
        return tuple(freeze_json(item) for item in value)
    return value


def thaw_json(value: Any) -> Any:
    """Return a fresh mutable JSON copy for compatibility adapter boundaries."""

    if isinstance(value, Mapping):
        return {str(key): thaw_json(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [thaw_json(item) for item in value]
    return value


_freeze = freeze_json


def _text(value: Any) -> str:
    if value is None:
        return "Unavailable — value not reported"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, Enum):
        return str(value.value)
    return str(value)


@dataclass(frozen=True, slots=True)
class PresentedSection:
    """One section of the GUI projection."""

    status: Availability
    summary: str
    rows: tuple[FrozenMap, ...] = ()
    data: FrozenMap = field(default_factory=FrozenMap)
    reason: str | None = None


@dataclass(frozen=True, slots=True)
class PresentationSections:
    summary: PresentedSection
    quality: PresentedSection
    msa: PresentedSection
    mutations: PresentedSection
    correspondences: PresentedSection
    interactions: PresentedSection
    sites: PresentedSection
    distance_map: PresentedSection
    vectors: PresentedSection
    evidence_cards: PresentedSection
    diagnostics: PresentedSection
    export: PresentedSection
    pymol: PresentedSection


@dataclass(frozen=True, slots=True)
class ReportPresentation:
    """Deterministic GUI data for one canonical report."""

    report_id: str
    sections: PresentationSections
    units: FrozenMap
    explanations: FrozenMap
    bundle_payloads: FrozenMap
    chart_datasets: FrozenMap


def _diagnostic_reason(report: AnalysisReport, name: str) -> str | None:
    state = getattr(report.availability, name, Availability.NOT_APPLICABLE)
    native_reason = {
        Availability.NOT_APPLICABLE: "not applicable for this report",
        Availability.NOT_DETECTED: "no evidence was detected",
        Availability.INVALID_INPUT: "the selected input is invalid",
        Availability.DEPENDENCY_UNAVAILABLE: "the required dependency is unavailable",
        Availability.NUMERICAL_FAILURE: "the calculation could not be completed",
    }.get(state)
    prefixes = (f"report.{name}",) if name else ()
    for diagnostic in report.diagnostics:
        if diagnostic.code.startswith(prefixes):
            suffix = f"; {native_reason}" if native_reason is not None else ""
            return f"{diagnostic.code}: {diagnostic.message}{suffix}"
    return native_reason


def _section(
    report: AnalysisReport,
    name: str,
    data: Mapping[str, Any] | None = None,
    rows: tuple[Mapping[str, Any], ...] = (),
    *,
    summary: str | None = None,
) -> PresentedSection:
    status = getattr(report.availability, name)
    reason = None if status is Availability.AVAILABLE else _diagnostic_reason(report, name)
    rendered = summary or (f"{len(rows)} record(s)" if rows else "Evidence available")
    if reason is not None:
        rendered = f"Unavailable — {reason}"
    return PresentedSection(
        status=status,
        summary=rendered,
        rows=tuple(FrozenMap({key: _text(value) for key, value in row.items()}) for row in rows),
        data=FrozenMap(_mapping(data)),
        reason=reason,
    )


def _mapping(value: Any) -> Mapping[str, Any]:
    if isinstance(value, Mapping):
        return value
    return {"items": value}


def _report_rows(report: AnalysisReport) -> dict[str, Any]:
    """Return fresh canonical JSON once; it is a projection, not a calculation."""

    payload = report.to_json()
    return payload


def present_report(report: AnalysisReport) -> ReportPresentation:
    """Project one report into immutable, deterministic GUI values."""

    if not isinstance(report, AnalysisReport):
        raise TypeError("report must be an AnalysisReport")
    payload = _report_rows(report)
    analysis = report.analysis
    analysis_payload = payload.get("analysis") or {}
    reference = report.reference_selection
    target = report.target_selection
    summary_rows: tuple[Mapping[str, Any], ...] = (
        {
            "reference": reference.display_name,
            "target": target.display_name,
            "decision": analysis.alignment_decision if analysis else None,
            "sequence_identity": analysis.sequence_identity if analysis else None,
            "sequence_similarity": analysis.sequence_similarity if analysis else None,
            "sequence_coverage": analysis.sequence_coverage if analysis else None,
            "strict_rmsd_angstrom": analysis.strict_rmsd_angstrom if analysis else None,
            "refined_rmsd_angstrom": analysis.refined_rmsd_angstrom if analysis else None,
            "tm_score": analysis.tm_score if analysis else None,
            "mapped_residue_count": analysis.mapped_residue_count if analysis else None,
        },
    )
    summary = _section(
        report,
        "analysis",
        {"reference": reference.display_name, "target": target.display_name, **analysis_payload},
        summary_rows,
        summary=(
            f"{reference.display_name} → {target.display_name} · "
            f"{analysis.mapped_residue_count if analysis else 'Unavailable'} mapped residue pair(s)"
        ),
    )
    quality = _section(
        report, "input_quality", payload["input_quality"], summary="Input quality is authoritative and descriptive."
    )

    msa_rows: tuple[Mapping[str, Any], ...] = ()
    if report.msa is not None:
        msa_rows = tuple(
            {
                "structure": structure_id,
                "aligned_sequence": row,
                "source": next(
                    (item.source for item in report.msa.sequences if item.structure_id == structure_id), "unknown"
                ),
            }
            for structure_id, row in report.msa.aligned_rows
        )
    msa = _section(report, "msa", payload.get("msa") or {}, msa_rows)

    mutation_payloads = tuple(analysis_payload.get("mutations", ()))
    mutation_rows = tuple(
        {
            "index": item.get("alignment_index"),
            "kind": item.get("kind"),
            "reference": item.get("reference_label"),
            "target": item.get("target_label"),
            "notation": item.get("canonical_notation"),
            "blosum62": item.get("blosum62_score"),
            "grantham": item.get("grantham_distance"),
            "class": item.get("physicochemical_class"),
        }
        for item in mutation_payloads
    )
    mutations = _section(
        report,
        "analysis",
        {"mutations": mutation_payloads},
        mutation_rows,
        summary=f"{len(mutation_rows)} mutation event(s)",
    )

    correspondence_payloads = tuple(analysis_payload.get("correspondences", ()))
    correspondence_rows = tuple(
        {
            "index": item.get("alignment_index"),
            "reference": item.get("reference"),
            "target": item.get("target"),
            "status": item.get("status"),
            "ca_displacement_angstrom": item.get("ca_displacement_angstrom"),
            "backbone_rmsd_angstrom": item.get("backbone_rmsd_angstrom"),
            "sidechain_rmsd_angstrom": item.get("sidechain_rmsd_angstrom"),
            "all_heavy_atom_rmsd_angstrom": item.get("all_heavy_atom_rmsd_angstrom"),
            "outlier": item.get("is_outlier"),
            "key": item.get("is_key_residue"),
        }
        for item in correspondence_payloads
    )
    correspondences = _section(
        report,
        "analysis",
        {"correspondences": correspondence_payloads},
        correspondence_rows,
        summary=f"{len(correspondence_rows)} aligned position(s)",
    )

    interactions = _section(
        report, "interactions", _mapping(payload.get("interactions") or {}), tuple(_interaction_rows(report))
    )
    sites = _section(report, "sites", _mapping(payload.get("sites") or {}), tuple(_site_rows(report)))
    distance = _section(report, "distance_map", payload.get("distance_map") or {}, _distance_rows(report))
    vectors = _section(
        report, "displacement_vectors", {"vectors": payload.get("displacement_vectors", ())}, _vector_rows(report)
    )
    cards = _section(report, "evidence_cards", {"cards": payload.get("evidence_cards", ())}, _card_rows(report))
    diagnostics = _section(
        report,
        "analysis",
        {"diagnostics": payload.get("diagnostics", ())},
        tuple(
            {
                "code": item.code,
                "severity": item.severity.value,
                "message": item.message,
                "remediation": item.remediation,
            }
            for item in report.diagnostics
        ),
        summary=f"{len(report.diagnostics)} diagnostic(s)",
    )
    export = _section(
        report,
        "analysis",
        {
            "report_id": report.report_id,
            "schema_version": report.schema_version,
            "input_hashes": payload.get("provenance", {}).get("input_hashes", {}) if payload.get("provenance") else {},
        },
        summary="Canonical JSON, CSV, TSV, and XLSX exports use this report.",
    )
    pymol = _section(
        report,
        "analysis",
        {
            "reference_selection": payload["reference_selection"],
            "target_selection": payload["target_selection"],
            "analysis": analysis_payload,
        },
        summary="PyMOL data is reversible and derived from the canonical report.",
    )

    bundles = {
        "msa_summary": {"status": report.availability.msa.value, "data": payload.get("msa")},
        "conservation": {
            "status": report.availability.msa.value,
            "columns": (payload.get("msa") or {}).get("columns", []),
        },
        "interactions": {"status": report.availability.interactions.value, "data": payload.get("interactions")},
        "sites": {"status": report.availability.sites.value, "data": payload.get("sites")},
        "evidence": {"status": report.availability.evidence_cards.value, "cards": payload.get("evidence_cards", [])},
        "vectors": {
            "status": report.availability.displacement_vectors.value,
            "vectors": payload.get("displacement_vectors", []),
        },
    }
    charts = {
        "msa_conservation_profile": {
            "unit": "fraction",
            "columns": tuple(
                {"index": column.index, "value": column.conservation_score} for column in report.msa.columns
            )
            if report.msa
            else (),
        },
        "structural_deviation_profile": {
            "unit": "Å",
            "points": tuple(
                {
                    "position": _numeric_residue_position(item.reference.auth_seq_id) if item.reference else None,
                    "label": (
                        f"{item.reference.chain_id}:{item.reference.auth_seq_id}{item.reference.insertion_code or ''}"
                        if item.reference
                        else "Unavailable"
                    ),
                    "value": item.ca_displacement_angstrom,
                    "metadata": {
                        "alignment_index": item.alignment_index,
                        "reference_residue": _residue_metadata(item.reference),
                        "target_residue": _residue_metadata(item.target),
                    },
                }
                for item in analysis.correspondences
            )
            if analysis
            else (),
        },
        "distance_difference_heatmap": {"unit": "Å", "data": payload.get("distance_map")},
    }
    return ReportPresentation(
        report.report_id,
        PresentationSections(
            summary,
            quality,
            msa,
            mutations,
            correspondences,
            interactions,
            sites,
            distance,
            vectors,
            cards,
            diagnostics,
            export,
            pymol,
        ),
        FrozenMap({"length": "Å", "area": "Å²", "volume": "Å³", "fraction": "fraction", "angle": "degrees"}),
        FrozenMap(
            {
                "method": "Values are projected from the canonical report; no GUI metric is recalculated.",
                "cα": "Cα distances are aligned residue displacements in Å.",
                "atomic_envelope": "Atomic envelope volume is the reported coordinate-envelope measure.",
                "pocket_free_volume": "Pocket free volume is reported by the validated pocket service in Å³.",
                "limits": "Unavailable values remain unavailable; method limits and diagnostics are retained.",
            }
        ),
        freeze_json(bundles),
        freeze_json(charts),
    )


def _interaction_rows(report: AnalysisReport) -> list[Mapping[str, Any]]:
    value = report.interactions
    records = value.differences if isinstance(value, InteractionEvidence) else value or ()
    rows: list[Mapping[str, Any]] = []
    for item in records:
        key = getattr(item, "key", None)
        row: dict[str, Any] = {
            "key": key,
            "change": getattr(item, "change", None),
            "interaction_type": key.interaction_type if key is not None else getattr(item, "interaction_type", None),
            "reference_position_a": getattr(key, "reference_position_a", None),
            "reference_position_b": getattr(key, "reference_position_b", None),
            "external_partner_id": getattr(key, "external_partner_id", None),
            "distance_unit": "Å",
            "angle_unit": "degrees",
        }
        for prefix, record in (
            ("reference", getattr(item, "reference_record", None)),
            ("target", getattr(item, "target_record", None)),
        ):
            row.update(
                {
                    f"{prefix}_distance_angstrom": getattr(record, "distance_angstrom", None),
                    f"{prefix}_angle_degrees": getattr(record, "angle_degrees", None),
                    f"{prefix}_residue_a": getattr(record, "residue_a", None),
                    f"{prefix}_residue_b": getattr(record, "residue_b", None),
                    f"{prefix}_atom_a": getattr(record, "atom_a", None),
                    f"{prefix}_atom_b": getattr(record, "atom_b", None),
                    f"{prefix}_evidence_mode": getattr(record, "evidence_mode", None),
                }
            )
        rows.append(row)
    return rows


def _site_rows(report: AnalysisReport) -> list[Mapping[str, Any]]:
    value = report.sites
    records = value.metrics if isinstance(value, SiteEvidence) else value or ()
    return [
        {
            "site_id": item.site_id,
            "structure_id": item.structure_id,
            "mapped_residue_count": item.mapped_residue_count,
            "coverage_fraction": item.coverage_fraction,
            "global_frame_backbone_rmsd_angstrom": item.global_frame_backbone_rmsd_angstrom,
            "site_fitted_backbone_rmsd_angstrom": item.site_fitted_backbone_rmsd_angstrom,
            "centroid_displacement_angstrom": item.centroid_displacement_angstrom,
            "radius_of_gyration_angstrom": item.radius_of_gyration_angstrom,
            "atomic_envelope_volume_angstrom3": item.atomic_envelope_volume_angstrom3,
            "sasa_angstrom2": item.sasa_angstrom2,
            "polar_residue_fraction": item.polar_residue_fraction,
            "charged_residue_fraction": item.charged_residue_fraction,
            "length_unit": "Å",
            "area_unit": "Å²",
            "volume_unit": "Å³",
            "composition_unit": "fraction",
        }
        for item in records
    ]


def _distance_rows(report: AnalysisReport) -> tuple[Mapping[str, Any], ...]:
    value = report.distance_map
    if value is None:
        return ()
    return tuple(
        {
            "position": position,
            "reference_distances_angstrom": reference,
            "target_distances_angstrom": target,
            "delta_angstrom": delta,
        }
        for position, reference, target, delta in zip(
            value.reference_positions,
            value.reference_distances_angstrom,
            value.target_distances_angstrom,
            value.delta_angstrom,
            strict=True,
        )
    )


def _vector_rows(report: AnalysisReport) -> tuple[Mapping[str, Any], ...]:
    return tuple(
        {
            "reference_position": item.reference_position,
            "reference_residue": item.reference_residue,
            "target_residue": item.target_residue,
            "vector_xyz": item.vector_xyz,
            "magnitude_angstrom": item.magnitude_angstrom,
        }
        for item in report.displacement_vectors
    )


def _card_rows(report: AnalysisReport) -> tuple[Mapping[str, Any], ...]:
    return tuple(
        {
            "reference_residue": item.reference_residue,
            "target_id": item.target_id,
            "alignment_index": item.sequence.alignment_index,
            "quality": item.quality.overall_status,
        }
        for item in report.evidence_cards
    )


def _numeric_residue_position(value: str) -> float | None:
    try:
        return float(value)
    except ValueError:
        return None


def _residue_metadata(residue: Any) -> Mapping[str, Any] | None:
    if residue is None:
        return None
    return {
        "structure_id": residue.structure_id,
        "model_id": residue.model_id,
        "chain_id": residue.chain_id,
        "auth_seq_id": residue.auth_seq_id,
        "insertion_code": residue.insertion_code,
        "residue_name": residue.residue_name,
    }


__all__ = [
    "FrozenMap",
    "PresentedSection",
    "PresentationSections",
    "ReportPresentation",
    "freeze_json",
    "present_report",
    "thaw_json",
]
