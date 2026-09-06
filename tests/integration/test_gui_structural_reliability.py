"""Off-screen Task 14 workflow contracts for one canonical report."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from PySide6 import QtWidgets  # noqa: E402

from structlens.application.dto import AnalysisReportRequest  # noqa: E402
from structlens.application.report_service import ReportService  # noqa: E402
from structlens.core.models import AnalysisSettings  # noqa: E402
from structlens.core.parsing import InputSelection, StructureFormat, capture_snapshot  # noqa: E402
from structlens.core.reports.pockets import PocketReportSnapshot  # noqa: E402
from structlens.plugin.gui.main_panel import build_qt_panel  # noqa: E402

_FIXTURE = Path("tests/fixtures/parsing/enriched.pdb")


@pytest.fixture(scope="module")
def application():
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    yield app
    app.quit()


def _request(tmp_path: Path) -> AnalysisReportRequest:
    source = _FIXTURE.read_bytes()
    reference_path = tmp_path / "reference.pdb"
    target_path = tmp_path / "target.pdb"
    reference_path.write_bytes(source)
    target_path.write_bytes(source)
    reference = capture_snapshot(reference_path)
    target = capture_snapshot(target_path)
    return AnalysisReportRequest(
        reference,
        target,
        InputSelection(
            reference.content_id,
            reference.display_name,
            StructureFormat.PDB,
            "1",
            author_chain_ids=("A",),
            path=str(reference_path),
        ),
        InputSelection(
            target.content_id,
            target.display_name,
            StructureFormat.PDB,
            "1",
            author_chain_ids=("A",),
            path=str(target_path),
        ),
        analysis_settings=AnalysisSettings(),
    )


def test_canonical_compare_populates_typed_pockets_and_gui_uses_same_report_id(application, tmp_path: Path) -> None:
    request = _request(tmp_path)
    report = ReportService().analyze(request)

    # Pocket evidence is part of the canonical compare artifact.  A separate
    # GUI-local detector would create a second, untraceable report identity.
    assert isinstance(report.pockets, PocketReportSnapshot)

    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    controller._populate_report(report, request=request)
    application.processEvents()

    assert controller.model.report is report
    assert controller._report_presentation is not None
    assert controller._report_presentation.report_id == report.report_id
    labels = {widget.text() for widget in panel.findChildren(QtWidgets.QLabel) if widget.text()}
    assert "Quality" in labels
    assert any("Pocket" in label or "Pockets" in label for label in labels)

    actions = [button.text().casefold() for button in panel.findChildren(QtWidgets.QPushButton)]
    assert any("detect" in action and "pocket" in action for action in actions)
    assert any("measure" in action and "volume" in action for action in actions)
    assert any("export" in action for action in actions)

    # Navigation and off-screen event processing must remain responsive after
    # the report is installed; this catches GUI work accidentally doing a
    # second scientific calculation synchronously.
    for index in range(controller.nav.count()):
        controller.nav.setCurrentRow(index)
        application.processEvents()
    assert controller.model.report.report_id == report.report_id
    controller.close()
    panel.deleteLater()
