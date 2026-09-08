"""Visible workflow checks for the existing scientific pages."""

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest
from PySide6 import QtWidgets

from structlens.plugin.gui.main_panel import build_qt_panel


@pytest.fixture
def controller():
    application = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    panel = build_qt_panel()
    yield panel._structlens_controller
    panel._structlens_controller.close()
    panel.deleteLater()
    application.processEvents()


def test_specialist_controls_collapse_without_losing_values(controller):
    assert controller.advanced_group.isHidden()
    assert not controller.advanced_toggle.isChecked()
    controller.advanced_toggle.click()
    assert not controller.advanced_group.isHidden()
    controller.identity_spin.setValue(0.42)
    controller.advanced_toggle.click()
    controller.advanced_toggle.click()
    assert controller.identity_spin.value() == 0.42
    assert controller.comparison_combo.isHidden()
    assert controller.comparison_combo.currentData() == "pairwise"
    assert controller.manual_group.isHidden()
    controller.mode_combo.setCurrentIndex(3)
    assert not controller.manual_group.isHidden()


def test_secondary_details_remain_accessible(controller):
    assert controller.backend_group.isHidden()
    controller.backend_toggle.click()
    assert not controller.backend_group.isHidden()
    assert controller.report_details_group.isHidden()
    assert controller.report_section_labels
    controller.report_details_toggle.click()
    assert not controller.report_details_group.isHidden()


def test_site_and_pocket_controls_describe_actual_actions(controller):
    assert controller.site_define_button.text() == "Add site to analysis"
    assert "next comparison" in controller.site_status_label.text().lower()
    assert controller.pocket_detect_button.isHidden()
    assert controller.pocket_measure_button.isHidden()
    assert controller.load_sources_button.text() == "Load sources"


def test_exports_fit_two_columns(controller):
    group = next(
        group for group in controller.widget.findChildren(QtWidgets.QGroupBox) if group.title() == "Evidence exports"
    )
    layout = group.layout()
    assert isinstance(layout, QtWidgets.QGridLayout)
    assert layout.columnCount() == 2
    assert layout.count() == 7
    assert len(controller.chart_export_buttons) == 3
