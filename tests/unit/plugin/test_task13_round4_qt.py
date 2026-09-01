"""Round-four RED tests for the real Qt presentation boundary."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest
from PySide6.QtWidgets import QApplication, QGroupBox, QLabel, QPushButton, QTableWidget

from structlens.application.report_service import ReportService
from structlens.core.evidence import Availability
from structlens.core.reports import DistanceMapSnapshot
from structlens.plugin.gui.main_panel import build_qt_panel

from ._task13_round4_fixtures import request, rich_report


@pytest.fixture(scope="module")
def application() -> QApplication:
    app = QApplication.instance() or QApplication([])
    yield app
    app.quit()


def _populate(controller: object, request_value: object, report: object) -> None:
    reference_path = Path(request_value.reference_selection.path)  # type: ignore[attr-defined]
    target_path = Path(request_value.target_selection.path)  # type: ignore[attr-defined]
    assert controller._load_source("reference", reference_path)  # type: ignore[attr-defined]
    assert controller._load_source("target", target_path)  # type: ignore[attr-defined]
    controller._report_request = request_value  # type: ignore[attr-defined]
    controller._analysis_finished(report)  # type: ignore[attr-defined]


def _widget_text(controller: object) -> str:
    chunks: list[str] = []
    widget = controller.widget  # type: ignore[attr-defined]
    for child in widget.findChildren(QLabel):
        chunks.append(child.text())
    for child in widget.findChildren(QGroupBox):
        chunks.append(child.title())
    for child in widget.findChildren(QPushButton):
        chunks.append(child.text())
    for table in widget.findChildren(QTableWidget):
        for row in range(table.rowCount()):
            for column in range(table.columnCount()):
                item = table.item(row, column)
                if item is not None:
                    chunks.append(item.text())
        for column in range(table.columnCount()):
            chunks.append(table.horizontalHeaderItem(column).text())
    return "\n".join(chunks)


def test_required_pyside6_binding_is_a_direct_release_gate() -> None:
    """This module must fail loudly if the official Qt binding is unavailable."""

    assert QApplication is not None


def test_qt_renders_every_available_report_section_from_one_presentation(application, tmp_path: Path) -> None:
    """Qt must expose report-owned interactions, maps, vectors, cards, and diagnostics."""

    request_value, report = rich_report(tmp_path)
    assert report.analysis is not None and report.analysis.correspondences
    distance = DistanceMapSnapshot(
        ("Z:900", "Z:901"),
        ((0.0, 2.5), (2.5, 0.0)),
        ((0.0, 3.75), (3.75, 0.0)),
        ((0.0, 1.25), (1.25, 0.0)),
        ((True, True), (True, True)),
    )
    report = replace(report, distance_map=distance)
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    _populate(controller, request_value, report)
    rendered = _widget_text(controller)

    assert "round4-site" in rendered
    assert "hbond_geometric" in rendered
    assert "Z:900" in rendered
    assert "5.0" in rendered or "5.000" in rendered
    assert "round4.integration.diagnostic" in rendered
    assert "12.5" in rendered or "12.500" in rendered
    assert "Evidence Card" in rendered
    assert "Distance" in rendered
    assert "Vector" in rendered
    controller.close()
    panel.deleteLater()


def test_qt_preserves_native_optional_failure_reason_while_rendering_independent_sections(
    application, tmp_path: Path
) -> None:
    """Unavailable MSA text must retain its diagnostic and not hide site data."""

    class FailingMSA:
        def align(self, *_args: object, **_kwargs: object) -> object:
            raise RuntimeError("MSA backend absent")

    request_value = request(tmp_path)
    report = ReportService(msa_engine=FailingMSA()).analyze(request_value)
    assert report.availability.msa is Availability.NUMERICAL_FAILURE
    assert report.msa is None
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    _populate(controller, request_value, report)
    rendered = _widget_text(controller)

    assert "Unavailable" in rendered
    assert "report.msa.failed" in rendered
    assert "could not be calculated" in rendered
    assert "active" in rendered
    controller.close()
    panel.deleteLater()


def test_qt_report_population_does_not_create_a_competing_legacy_analysis(application, tmp_path: Path) -> None:
    request_value, report = rich_report(tmp_path)
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    _populate(controller, request_value, report)

    assert controller.model.report is report
    assert controller.model.analysis is None
    assert controller._report_presentation is not None
    controller.close()
    panel.deleteLater()
