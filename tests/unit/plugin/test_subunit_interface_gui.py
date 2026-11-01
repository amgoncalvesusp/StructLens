"""Exercise interchain contacts from the standalone Project page."""

import csv
import os
import time
from dataclasses import replace
from threading import Event

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from PySide6 import QtCore, QtTest, QtWidgets  # noqa: E402

from structlens.plugin.gui.main_panel import build_qt_panel  # noqa: E402
from tests.unit.plugin.test_variant_batch import make_source  # noqa: E402


@pytest.fixture
def panel():
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    widget = build_qt_panel(command=None)
    widget.show()
    yield widget._structlens_controller
    widget._structlens_controller.close()
    widget.close()
    widget.deleteLater()
    app.processEvents()


def click(button):
    QtTest.QTest.mouseClick(button, QtCore.Qt.MouseButton.LeftButton)
    QtWidgets.QApplication.processEvents()


def wait_result(window):
    end = time.monotonic() + 15
    while window.future and time.monotonic() < end:
        QtWidgets.QApplication.processEvents()
        time.sleep(.015)
    assert window.future is None


def test_interfaces_from_reference_only_selects_distinct_chains_and_exports(panel, tmp_path, monkeypatch):
    source = make_source(tmp_path / "complex.pdb", ("A", "B"))
    assert panel._load_source("reference", source)
    click(panel.interface_button)
    window = panel.interface_window
    assert window.chain_a.currentData() == "A" and window.chain_b.currentData() == "B"
    click(window.run_button)
    wait_result(window)
    assert window.result is not None and window.result.contacts
    assert window.pairs.rowCount() > 0 and window.residues.rowCount() > 0
    output = tmp_path / "interface.csv"
    monkeypatch.setattr(QtWidgets.QFileDialog, "getSaveFileName", lambda *_: (str(output), ""))
    click(window.export_button)
    with output.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["source_sha256"] == panel.reference_loaded_source.snapshot.raw_sha256
    pair_rows = [r for r in rows if r["record"] == "pair"]
    assert len(pair_rows) == len(window.result.contacts)
    assert {r["chain_a"] for r in pair_rows} == {"A"}
    assert {r["chain_b"] for r in pair_rows} == {"B"}
    window.cutoff.setValue(4.0)
    assert window.result is None and not window.export_button.isEnabled()


def test_single_chain_cannot_run_or_export(panel, tmp_path):
    assert panel._load_source("reference", make_source(tmp_path / "single.pdb"))
    click(panel.interface_button)
    window = panel.interface_window
    assert not window.run_button.isEnabled() and not window.export_button.isEnabled()
    assert "two" in window.status.text().lower()


def test_nucleic_acid_chain_is_not_offered(panel, tmp_path):
    path = make_source(tmp_path / "mixed.pdb")
    dna = "ATOM      5    P  DA D   1       3.000   0.000   0.000  1.00 20.00           P\n"
    path.write_text(path.read_text().replace("END\n", dna + "END\n"))
    assert panel._load_source("reference", path)
    click(panel.interface_button)
    interface = panel.interface_window
    assert interface.chain_a.count() == 1 and interface.chain_a.currentData() == "A"
    assert not interface.run_button.isEnabled()


def test_selected_residue_drills_down_by_chain_and_author_id(panel, tmp_path):
    assert panel._load_source("reference", make_source(tmp_path / "dimer.pdb", ("A", "B")))
    click(panel.interface_button)
    window = panel.interface_window
    click(window.run_button)
    wait_result(window)
    window.tabs.setCurrentIndex(1)
    window.residues.setCurrentCell(0, 0)
    label = window.residues.item(0, 0).text()
    click(window.partners_button)
    assert window.tabs.currentIndex() == 0
    assert 0 < window.pairs.rowCount() < len(window.result.contacts)
    assert all(label in (window.pairs.item(i, 0).text(), window.pairs.item(i, 1).text())
               for i in range(window.pairs.rowCount()))
    click(window.clear_filter_button)
    assert window.pairs.rowCount() == len(window.result.contacts)


def test_closed_dialog_does_not_install_an_inflight_result(panel, tmp_path, monkeypatch):
    from structlens.plugin.gui import subunit_interface as gui

    assert panel._load_source("reference", make_source(tmp_path / "complex.pdb", ("A", "B")))
    click(panel.interface_button)
    old = panel.interface_window
    entered, release = Event(), Event()
    original = gui.analyze_subunit_interface

    def delayed(*args, **kwargs):
        entered.set()
        assert release.wait(10)
        return original(*args, **kwargs)

    monkeypatch.setattr(gui, "analyze_subunit_interface", delayed)
    click(old.run_button)
    try:
        assert entered.wait(3)
        assert not old.cutoff.isEnabled()
        old.dialog.close()
        assert old.cancel_event.is_set()
        click(panel.interface_button)
        fresh = panel.interface_window
        assert fresh is not old and fresh.result is None
    finally:
        release.set()
    click(fresh.run_button)
    wait_result(fresh)
    assert old.result is None
    assert fresh.result and fresh.result.contacts


def test_interface_csv_escapes_external_formula_tokens(panel, tmp_path):
    from structlens.plugin.gui.subunit_interface import export_interface_csv

    assert panel._load_source("reference", make_source(tmp_path / "complex.pdb", ("A", "B")))
    click(panel.interface_button)
    window = panel.interface_window
    click(window.run_button)
    wait_result(window)
    result = window.result
    original = result.contacts[0]
    malicious = replace(original, atom_a='=danger', residue_a=replace(original.residue_a, chain_id='@user'))
    path = tmp_path / "escaped.csv"
    export_interface_csv(path, replace(result, contacts=(malicious,)))
    with path.open(encoding="utf-8-sig", newline="") as handle:
        pair = next(row for row in csv.DictReader(handle) if row["record"] == "pair")
    assert pair["atom_a"] == "'=danger"
    assert pair["chain_a"] == "'@user"
