"""Reproducible tabular exports for scientific results."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict
from pathlib import Path
from typing import Any

from openpyxl import Workbook
from openpyxl.styles import Font
from PIL import Image

from structlens.core.models import AnalysisResult


def export_analysis_xlsx(result: AnalysisResult, path: str | Path) -> None:
    workbook = Workbook()
    summary = workbook.active
    if summary is None:
        raise RuntimeError("Workbook did not create a summary sheet")
    summary.title = "Summary"
    summary.append(["Metric", "Value", "Units"])
    summary.append(["Sequence identity", result.sequence_identity, "fraction"])
    summary.append(["Sequence coverage", result.sequence_coverage, "fraction"])
    summary.append(["Strict Cα RMSD", result.strict_rmsd_angstrom, "Å"])
    summary.append(["Refined Cα RMSD", result.refined_rmsd_angstrom, "Å"])
    summary.append(["Mapped residues", result.mapped_residue_count, "residues"])
    summary.append(["Mutation count", result.mutation_count, "events"])
    summary.append(["Alignment decision", result.alignment_decision, ""])
    for cell in summary[1]:
        cell.font = Font(bold=True)

    residues = workbook.create_sheet("Residues")
    residues.append(
        [
            "Reference residue",
            "Target residue",
            "Status",
            "Mutation notation",
            "Cα displacement (Å)",
            "Backbone RMSD (Å)",
            "Side-chain RMSD (Å)",
            "All-heavy-atom RMSD (Å)",
            "Outlier",
            "Key residue",
        ]
    )
    by_index = {
        event.alignment_index: event.canonical_notation for event in result.mutations
    }
    for item in result.correspondences:
        residues.append(
            [
                _residue_label(item.reference),
                _residue_label(item.target),
                item.status.value,
                by_index.get(item.alignment_index, ""),
                item.ca_displacement_angstrom,
                item.backbone_rmsd_angstrom,
                item.sidechain_rmsd_angstrom,
                item.all_heavy_atom_rmsd_angstrom,
                item.is_outlier,
                item.is_key_residue,
            ]
        )

    mutations = workbook.create_sheet("Mutations")
    mutations.append(
        [
            "Alignment index",
            "Kind",
            "Reference",
            "Target",
            "Notation",
            "BLOSUM62",
            "Grantham distance",
            "Physicochemical class",
        ]
    )
    for event in result.mutations:
        mutations.append(
            [
                event.alignment_index,
                event.kind.value,
                _residue_label(event.reference),
                _residue_label(event.target),
                event.canonical_notation,
                event.blosum62_score,
                event.grantham_distance,
                event.physicochemical_class,
            ]
        )
    for sheet in workbook.worksheets:
        sheet.freeze_panes = "A2"
        sheet.auto_filter.ref = sheet.dimensions
    workbook.save(path)


def export_analysis_csv(result: AnalysisResult, path: str | Path) -> None:
    rows = [
        {
            "reference_residue": _residue_label(item.reference),
            "target_residue": _residue_label(item.target),
            "status": item.status.value,
            "ca_displacement_angstrom": item.ca_displacement_angstrom,
            "backbone_rmsd_angstrom": item.backbone_rmsd_angstrom,
            "outlier": item.is_outlier,
        }
        for item in result.correspondences
    ]
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]) if rows else ["status"]
        )
        writer.writeheader()
        writer.writerows(rows)


def export_analysis_json(result: AnalysisResult, path: str | Path) -> None:
    payload: dict[str, Any] = {
        "reference_id": result.reference_id,
        "target_id": result.target_id,
        "sequence_identity": result.sequence_identity,
        "sequence_coverage": result.sequence_coverage,
        "alignment_decision": result.alignment_decision,
        "strict_rmsd_angstrom": result.strict_rmsd_angstrom,
        "refined_rmsd_angstrom": result.refined_rmsd_angstrom,
        "mutations": [asdict(event) for event in result.mutations],
    }
    Path(path).write_text(
        json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8"
    )


def export_publication_image(
    image: Image.Image,
    path: str | Path,
    *,
    width_mm: float = 100.0,
    dpi: int = 300,
) -> None:
    """Render a supplied molecular scene at publication dimensions and DPI."""

    if dpi not in {300, 600}:
        raise ValueError("Publication image dpi must be 300 or 600")
    if width_mm <= 0:
        raise ValueError("width_mm must be positive")
    output_path = Path(path)
    suffix = output_path.suffix.lower()
    if suffix not in {".jpg", ".jpeg", ".tif", ".tiff"}:
        raise ValueError("Publication image format must be JPEG or TIFF")
    if suffix in {".jpg", ".jpeg"} and (
        "A" in image.getbands() or image.mode in {"RGBA", "LA"}
    ):
        raise ValueError(
            "JPEG export does not support transparency; provide an opaque image"
        )
    width_px = round(width_mm / 25.4 * dpi)
    height_px = max(1, round(width_px * image.height / image.width))
    rendered = image.resize((width_px, height_px), Image.Resampling.LANCZOS)
    if suffix in {".jpg", ".jpeg"}:
        rendered = rendered.convert("RGB")
        rendered.save(
            output_path, format="JPEG", quality=95, dpi=(dpi, dpi), optimize=True
        )
    else:
        rendered.save(
            output_path, format="TIFF", compression="tiff_deflate", dpi=(dpi, dpi)
        )


def _residue_label(residue: Any) -> str:
    if residue is None:
        return ""
    insertion = residue.insertion_code or ""
    return f"{residue.chain_id}:{residue.auth_seq_id}{insertion} {residue.residue_name}"


__all__ = [
    "export_analysis_csv",
    "export_analysis_json",
    "export_analysis_xlsx",
    "export_publication_image",
]
