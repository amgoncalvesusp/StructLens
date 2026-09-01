"""Bounded XLSX helpers shared by the report export facade."""

from __future__ import annotations

import io
import re
import zipfile
from collections.abc import Callable, Iterable
from datetime import datetime
from typing import Any

from openpyxl import Workbook


def append_safe(sheet: Any, values: Iterable[Any], safe: Callable[[Any], Any]) -> None:
    sheet.append([safe(value) for value in values])


def normalise_workbook(workbook: Workbook) -> None:
    stamp = datetime(1980, 1, 1)
    workbook.properties.created = stamp
    workbook.properties.modified = stamp
    workbook.properties.creator = "StructLens"
    workbook.properties.lastModifiedBy = "StructLens"
    workbook.calculation.fullCalcOnLoad = False
    workbook.calculation.forceFullCalc = False


def deterministic_xlsx(data: bytes) -> bytes:
    """Rewrite ZIP member timestamps/ordering so identical reports hash alike."""

    output = io.BytesIO()
    with (
        zipfile.ZipFile(io.BytesIO(data), "r") as source,
        zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as target,
    ):
        for info in sorted(source.infolist(), key=lambda item: item.filename):
            member = zipfile.ZipInfo(info.filename, date_time=(1980, 1, 1, 0, 0, 0))
            member.compress_type = zipfile.ZIP_DEFLATED
            member.external_attr = 0o600 << 16
            content = source.read(info.filename)
            if info.filename == "docProps/core.xml":
                content = re.sub(
                    rb"(<dcterms:modified\b[^>]*>)[^<]*(</dcterms:modified>)",
                    rb"\g<1>1980-01-01T00:00:00Z\g<2>",
                    content,
                )
            target.writestr(member, content)
    return output.getvalue()


__all__ = ["append_safe", "deterministic_xlsx", "normalise_workbook"]
