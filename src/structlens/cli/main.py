"""English-only StructLens command-line entry point."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

from structlens.cli.commands import WorkflowResult, run_compare, run_pocket_compare, run_pockets, run_qc
from structlens.core.models import AlignmentMode
from structlens.core.parsing import AltlocPolicy, AssemblyScope


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="structlens",
        description="Reproducible sequence, structure, quality, and pocket analysis for proteins.",
    )
    subparsers = parser.add_subparsers(dest="command")
    _add_compare(subparsers)
    _add_single_structure(subparsers, "qc", "assess coordinate quality", include_json=True)
    _add_single_structure(subparsers, "pockets", "detect and measure pockets", include_json=True)
    _add_pocket_compare(subparsers)
    return parser


def main(argv: list[str] | None = None) -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    parser = build_parser()
    try:
        args = parser.parse_args(argv)
    except SystemExit as error:
        code = error.code
        return code if isinstance(code, int) else 2
    if args.command is None:
        parser.print_help()
        return 0
    result = _dispatch(args)
    if result.error is not None:
        print(f"Error: {result.error}", file=sys.stderr)
        return result.exit_code
    if result.summary:
        stream = sys.stderr if result.exit_code != 0 else sys.stdout
        print(result.summary, file=stream)
    return result.exit_code


def _dispatch(args: argparse.Namespace) -> WorkflowResult:
    if args.command == "compare":
        return run_compare(
            args.reference,
            args.target,
            mode=args.mode,
            reference_model=args.reference_model,
            target_model=args.target_model,
            reference_author_chains=tuple(args.reference_author_chain),
            target_author_chains=tuple(args.target_author_chain),
            reference_label_chains=tuple(args.reference_label_chain),
            target_label_chains=tuple(args.target_label_chain),
            reference_altloc_policy=args.reference_altloc_policy,
            target_altloc_policy=args.target_altloc_policy,
            reference_assembly_scope=args.reference_assembly_scope,
            target_assembly_scope=args.target_assembly_scope,
            json_path=args.json,
            report_path=args.report,
            csv_path=args.csv,
            tsv_path=args.tsv,
            xlsx_path=args.output,
        )
    if args.command == "qc":
        return run_qc(
            args.structure,
            model=args.model,
            author_chains=tuple(args.author_chain),
            label_chains=tuple(args.label_chain),
            altloc_policy=args.altloc_policy,
            assembly_scope=args.assembly_scope,
            json_path=args.json,
        )
    if args.command == "pockets":
        return run_pockets(
            args.structure,
            model=args.model,
            author_chains=tuple(args.author_chain),
            label_chains=tuple(args.label_chain),
            altloc_policy=args.altloc_policy,
            assembly_scope=args.assembly_scope,
            json_path=args.json,
        )
    if args.command == "pocket-compare":
        return run_pocket_compare(
            args.reference,
            args.target,
            mode=args.mode,
            reference_model=args.reference_model,
            target_model=args.target_model,
            reference_author_chains=tuple(args.reference_author_chain),
            target_author_chains=tuple(args.target_author_chain),
            reference_label_chains=tuple(args.reference_label_chain),
            target_label_chains=tuple(args.target_label_chain),
            reference_altloc_policy=args.reference_altloc_policy,
            target_altloc_policy=args.target_altloc_policy,
            reference_assembly_scope=args.reference_assembly_scope,
            target_assembly_scope=args.target_assembly_scope,
            json_path=args.json,
        )
    raise ValueError(f"unknown command: {args.command}")


def _add_compare(subparsers: argparse._SubParsersAction[Any]) -> None:
    compare = subparsers.add_parser("compare", help="compare one reference and one target")
    compare.add_argument("reference", type=Path)
    compare.add_argument("target", type=Path)
    compare.add_argument("--mode", choices=[mode.value for mode in AlignmentMode], default=AlignmentMode.AUTO.value)
    compare.add_argument("--output", type=Path, help="write a canonical XLSX report")
    compare.add_argument("--csv", type=Path, help="write a canonical CSV report")
    compare.add_argument("--tsv", type=Path, help="write a canonical TSV report")
    report_group = compare.add_mutually_exclusive_group()
    report_group.add_argument("--json", type=Path, help="write a canonical JSON report")
    report_group.add_argument("--report", type=Path, help="alias for canonical JSON report output")
    _add_compare_selection(compare, "reference")
    _add_compare_selection(compare, "target")


def _add_pocket_compare(subparsers: argparse._SubParsersAction[Any]) -> None:
    command = subparsers.add_parser("pocket-compare", help="compare detected pockets in two structures")
    command.add_argument("reference", type=Path)
    command.add_argument("target", type=Path)
    command.add_argument("--mode", choices=[mode.value for mode in AlignmentMode], default=AlignmentMode.AUTO.value)
    command.add_argument("--json", type=Path, help="write the typed pocket-comparison JSON envelope")
    _add_compare_selection(command, "reference")
    _add_compare_selection(command, "target")


def _add_compare_selection(parser: argparse.ArgumentParser, role: str) -> None:
    parser.add_argument(f"--{role}-model", default="1", type=_nonempty_argument("model"))
    parser.add_argument(f"--{role}-author-chain", action="append", default=[], type=_nonempty_argument("author chain"))
    parser.add_argument(f"--{role}-label-chain", action="append", default=[], type=_nonempty_argument("label chain"))
    parser.add_argument(
        f"--{role}-altloc-policy",
        choices=[policy.value for policy in AltlocPolicy],
        default=AltlocPolicy.HIGHEST_OCCUPANCY.value,
    )
    parser.add_argument(
        f"--{role}-assembly-scope",
        choices=[scope.value for scope in AssemblyScope],
        default=AssemblyScope.ASYMMETRIC_UNIT.value,
    )


def _add_single_structure(
    subparsers: argparse._SubParsersAction[Any],
    name: str,
    help_text: str,
    *,
    include_json: bool,
) -> None:
    command = subparsers.add_parser(name, help=help_text)
    command.add_argument("structure", type=Path)
    command.add_argument("--model", default="1", type=_nonempty_argument("model"))
    command.add_argument("--author-chain", action="append", default=[], type=_nonempty_argument("author chain"))
    command.add_argument("--label-chain", action="append", default=[], type=_nonempty_argument("label chain"))
    command.add_argument(
        "--altloc-policy",
        choices=[policy.value for policy in AltlocPolicy],
        default=AltlocPolicy.HIGHEST_OCCUPANCY.value,
    )
    command.add_argument(
        "--assembly-scope",
        choices=[scope.value for scope in AssemblyScope],
        default=AssemblyScope.ASYMMETRIC_UNIT.value,
    )
    if include_json:
        command.add_argument("--json", type=Path, help="write a deterministic JSON result")


def _nonempty_argument(label: str) -> Callable[[str], str]:
    def parse(value: str) -> str:
        if not value.strip():
            raise argparse.ArgumentTypeError(f"{label} must not be empty")
        return value.strip()

    return parse


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = ["build_parser", "main"]
