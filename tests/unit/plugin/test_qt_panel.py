import os
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
qt_widgets = pytest.importorskip("PySide6.QtWidgets")

from structlens.plugin.gui.main_panel import build_qt_panel  # noqa: E402


@pytest.fixture(scope="module")
def application():
    app = qt_widgets.QApplication.instance() or qt_widgets.QApplication([])
    yield app
    app.quit()


def test_qt_panel_builds_operate_workflow(application) -> None:
    panel = build_qt_panel()
    controller = panel._structlens_controller

    assert panel.objectName() == "structlensPanel"
    assert controller.nav.count() == 6
    assert controller.pages.count() == 6
    assert controller.compare_button.text() == "Compare"
    assert controller.mode_combo.count() == 4
    assert controller.mutation_table.columnCount() == 8
    assert controller.residue_table.columnCount() == 10

    controller.close()
    panel.deleteLater()


def test_failed_source_reload_clears_previous_structure(application) -> None:
    panel = build_qt_panel()
    controller = panel._structlens_controller
    fixture = Path("tests/fixtures/parsing/numbering_altloc.pdb")
    assert controller._load_source("reference", fixture) is True
    assert controller.reference_structure is not None
    controller._show_error = lambda _: None

    assert controller._load_source("reference", Path("missing-input.pdb")) is False
    assert controller.reference_structure is None
    assert controller.reference_chain_combo.count() == 0

    controller.close()
    panel.deleteLater()
