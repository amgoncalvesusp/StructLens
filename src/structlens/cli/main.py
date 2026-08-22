"""English-only command-line interface for StructLens."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from structlens.application.analysis_service import AnalysisService
from structlens.application.export_service import (
    export_analysis_csv,
    export_analysis_json,
    export_analysis_xlsx,
)
from structlens.core.errors import StructLensError
from structlens.core.models import AlignmentMode, AnalysisSettings


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="structlens",
        description="Reproducible sequence and structure comparison for proteins.",
    )
    subparsers = parser.add_subparsers(dest="command")
    compare = subparsers.add_parser(
        "compare", help="compare one reference and one target"
    )
    compare.add_argument("reference", type=Path)
    compare.add_argument("target", type=Path)
    compare.add_argument(
        "--mode",
        choices=[mode.value for mode in AlignmentMode],
        default=AlignmentMode.AUTO.value,
        help="residue mapping mode",
    )
    compare.add_argument("--output", type=Path, help="write an XLSX result workbook")
    compare.add_argument("--csv", type=Path, help="write a residue CSV")
    compare.add_argument("--json", type=Path, help="write a JSON summary")
    return parser


def main(argv: list[str] | None = None) -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    parser = build_parser()
    try:
        try:
            args = parser.parse_args(argv)
        except SystemExit as exc:
            return exc.code if isinstance(exc.code, int) else 2
        if args.command != "compare":
            parser.print_help()
            return 0
        settings = AnalysisSettings(alignment_mode=AlignmentMode(args.mode))
        result = AnalysisService().analyze_paths(args.reference, args.target, settings)
        print(f"Reference: {result.reference_id}")
        print(f"Target: {result.target_id}")
        print(f"Sequence identity: {result.sequence_identity:.3f}")
        print(f"Sequence coverage: {result.sequence_coverage:.3f}")
        print(f"Strict Cα RMSD: {_format_metric(result.strict_rmsd_angstrom)} Å")
        print(f"Refined Cα RMSD: {_format_metric(result.refined_rmsd_angstrom)} Å")
        print(f"Mapped residues: {result.mapped_residue_count}")
        print(f"Mutations: {result.mutation_count}")
        print(f"Alignment: {result.alignment_decision}")
        if args.output:
            export_analysis_xlsx(result, args.output)
        if args.csv:
            export_analysis_csv(result, args.csv)
        if args.json:
            export_analysis_json(result, args.json)
        return 0
    except (StructLensError, OSError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2


def _format_metric(value: float | None) -> str:
    return "unavailable" if value is None else f"{value:.3f}"


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = ["build_parser", "main"]
