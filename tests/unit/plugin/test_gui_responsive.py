"""Visible controls and native scrolling on small desktop screens."""

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from PySide6 import QtCore, QtWidgets  # noqa: E402

from structlens import __version__  # noqa: E402
from structlens.plugin.gui.main_panel import build_qt_panel  # noqa: E402


@pytest.fixture(scope="module")
def application():
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


@pytest.mark.parametrize("width,height", [(800, 600), (1024, 600), (640, 480), (1440, 900)])
def test_every_page_and_primary_action_remain_reachable(application, width, height):
    panel = build_qt_panel(command=None)
    panel.setAttribute(QtCore.Qt.WidgetAttribute.WA_DontShowOnScreen, True)
    controller = panel._structlens_controller
    try:
        panel.resize(width, height)
        panel.show()
        application.processEvents()
        assert panel.width() <= width
        assert panel.height() <= height
        assert controller.compare_button.text() == "Compare structures"
        assert controller.workflow_context.wordWrap()
        assert __version__ in controller.footer_meta.text()
        for index in range(controller.nav.count()):
            controller.nav.setCurrentRow(index)
            controller.nav.scrollToItem(controller.nav.item(index))
            application.processEvents()
            assert controller.pages.currentIndex() == index
            assert controller.nav.viewport().rect().intersects(
                controller.nav.visualItemRect(controller.nav.item(index))
            )
            button_position = controller.compare_button.mapTo(panel, QtCore.QPoint())
            assert panel.rect().contains(QtCore.QRect(button_position, controller.compare_button.size()))
            scroll = controller.pages.currentWidget().findChild(QtWidgets.QScrollArea)
            assert scroll.viewport().width() > 100
            assert scroll.viewport().height() > 100
    finally:
        controller.close()
        panel.close()
        panel.deleteLater()


def test_overflow_can_scroll_to_bottom_right(application):
    panel = build_qt_panel(command=None)
    panel.setAttribute(QtCore.Qt.WidgetAttribute.WA_DontShowOnScreen, True)
    controller = panel._structlens_controller
    try:
        panel.resize(800, 600)
        scroll = controller.pages.widget(0).findChild(QtWidgets.QScrollArea)
        scroll.widget().setMinimumSize(1400, 1200)
        panel.show()
        application.processEvents()
        for bar in (scroll.horizontalScrollBar(), scroll.verticalScrollBar()):
            assert bar.isVisible()
            assert bar.maximum() > 0
            bar.setValue(bar.maximum())
            assert bar.value() == bar.maximum()
        assert scroll.widget().palette().window().color().lightness() < 64
    finally:
        controller.close()
        panel.close()
        panel.deleteLater()


def test_long_source_details_wrap_without_widening_page(application):
    panel = build_qt_panel(command=None)
    panel.setAttribute(QtCore.Qt.WidgetAttribute.WA_DontShowOnScreen, True)
    controller = panel._structlens_controller
    try:
        panel.resize(800, 600)
        panel.show()
        application.processEvents()
        page = controller.pages.widget(0).findChild(QtWidgets.QScrollArea).widget()
        original_width = page.minimumSizeHint().width()
        controller.reference_meta.setText("Loaded source with multiple models and chains. " * 30)
        application.processEvents()
        assert page.minimumSizeHint().width() <= max(original_width, 800)
    finally:
        controller.close()
        panel.close()
        panel.deleteLater()


@pytest.mark.parametrize("width,height", [(800, 600), (1024, 600), (1920, 1080)])
def test_initial_geometry_uses_available_logical_screen(application, width, height):
    from structlens.gui.main import _fit_window_to_screen

    available = QtCore.QRect(120, 50, width, height)

    class Screen:
        def availableGeometry(self):
            return available

    class Window(QtWidgets.QWidget):
        def screen(self):
            return Screen()

    panel = Window()
    try:
        _fit_window_to_screen(panel)
        assert available.contains(panel.geometry())
        assert panel.width() <= 1280
        assert panel.height() <= 820
        assert panel.height() < height
    finally:
        panel.deleteLater()
