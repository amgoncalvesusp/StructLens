"""Deterministic tabular rows for aggregate and concordance pocket evidence."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

from structlens.core.reports import PocketReportSnapshot


@dataclass(frozen=True, slots=True)
class PocketExportRow:
    """One lossless row specification consumed by every tabular writer."""

    section: str
    role: str
    metric: str
    value: Any
    units: str
    status: str
    provenance: Any = ""


def _display_units(units: dict[str, str], serialize: Callable[[Any], str]) -> str:
    if not units:
        return ""
    if len(units) == 1:
        return next(iter(units.values()))
    return serialize(units)


def pocket_export_rows(
    pockets: PocketReportSnapshot,
    *,
    serialize: Callable[[Any], str],
) -> tuple[PocketExportRow, ...]:
    """Return aggregate units, native concordance values, and diagnostics."""

    status = pockets.availability.value if pockets.availability is not None else "not_applicable"
    rows = [
        PocketExportRow(
            "pockets",
            "report",
            "pocket_report_units",
            dict(pockets.units),
            "",
            status,
        )
    ]
    channels = () if pockets.concordance is None else sorted(pockets.concordance.to_mapping().items())
    for name, channel in channels:
        channel_status = channel.availability.value
        channel_payload = channel.to_json()
        provenance = "" if channel_payload["provenance"] is None else channel_payload["provenance"]
        rows.append(
            PocketExportRow(
                "pockets",
                "concordance",
                name,
                channel_payload["measure"],
                _display_units(dict(channel.units), serialize),
                channel_status,
                provenance,
            )
        )
        rows.extend(
            PocketExportRow(
                "diagnostics",
                f"concordance:{name}",
                diagnostic.code,
                diagnostic.message,
                "",
                channel_status,
                provenance,
            )
            for diagnostic in channel.diagnostics
        )
    return tuple(rows)


__all__ = ["PocketExportRow", "pocket_export_rows"]
