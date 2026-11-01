"""Desktop reference-versus-variants workflow through the real Project action."""

import importlib.util
import os
import subprocess
import sys
import time
from threading import Event

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from PySide6 import QtWidgets  # noqa: E402

from structlens.application.report_serialization import load_report  # noqa: E402
from structlens.plugin.gui.main_panel import build_qt_panel  # noqa: E402


@pytest.fixture
def batch_panel(tmp_path):
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    panel = build_qt_panel(command=None)
    control = panel._structlens_controller
    ref = make_source(tmp_path / "reference.pdb")
    assert control._load_source("reference", ref)
    control.mode_combo.setCurrentText("Sequence")
    control.variant_batch_button.click()
    yield control, control.variant_batch, ref
    control.close()
    panel.deleteLater()
    app.processEvents()


def make_source(path, chains=("A",), third="SER"):
    lines = []
    for chain in chains:
        for i, aa in enumerate(("ALA", "VAL", third, "GLY"), 1):
            lines.append(f"ATOM  {len(lines)+1:5d}  CA  {aa} {chain}{i:4d}    "
                         f"{i*3.8:8.3f}{float(i%2):8.3f}{float(i%3):8.3f}  1.00 20.00           C\n")
    path.write_text("".join(lines) + "END\n", encoding="ascii")
    return path


def wait_idle(batch):
    deadline = time.monotonic() + 30
    while batch.busy and time.monotonic() < deadline:
        QtWidgets.QApplication.processEvents()
        time.sleep(.01)
    assert not batch.busy


@pytest.mark.parametrize("binding", ["PySide6", "PyQt5"])
def test_batch_tooltips_display_untrusted_paths_and_errors_as_text(binding):
    script = '''
from importlib import import_module
from types import SimpleNamespace
import sys

from structlens.plugin.gui.variant_batch import VariantBatch
widgets = import_module(sys.argv[1] + ".QtWidgets")
gui = import_module(sys.argv[1] + ".QtGui")
app = widgets.QApplication([])
table = widgets.QTableWidget(1, 1)
batch = SimpleNamespace(panel=SimpleNamespace(w=widgets), table=table)
payload = """a & <img src="local-test.png"> 'line'"""
for label, tip, expected in [("variant", payload, payload), ("Failed: " + payload, "", "Failed: " + payload)]:
    VariantBatch.cell(batch, 0, 0, label, tip)
    item = table.item(0, 0)
    assert item.text() == label and "<img" not in item.toolTip()
    rendered = gui.QTextDocument()
    rendered.setHtml(item.toolTip())
    assert rendered.toPlainText() == expected
'''
    if importlib.util.find_spec(binding) is None:
        pytest.skip(f"{binding} not installed")
    subprocess.run([sys.executable, "-c", script, binding], check=True, timeout=60)


def test_variant_batch_reports_scope_and_open_exact_export(batch_panel, tmp_path, monkeypatch):
    control, batch, ref = batch_panel
    changed = make_source(tmp_path / "variant.pdb", ("A", "B"), third="THR")
    bad = tmp_path / "missing.pdb"
    batch.add_paths([ref, changed, bad])
    wait_idle(batch)
    assert batch.table.rowCount() == 3
    assert batch.table.cellWidget(1, 1).currentData() is None
    assert "Failed" in batch.table.item(2, 2).text()
    batch.table.cellWidget(1, 1).setCurrentIndex(2)  # model 0, chain B
    batch.run_button.click()
    wait_idle(batch)
    assert batch.table.item(0, 6).text() == "0"
    assert batch.table.item(1, 6).text() == "1"
    assert batch.table.item(1, 3).text() == "100.0%"
    assert batch.table.item(1, 4).text() == "4"
    assert "Failed" in batch.table.item(2, 2).text()
    batch.table.cellClicked.emit(1, 0)
    binding = control._canonical_binding()
    assert binding.request.target_selection.author_chain_ids == ("B",)
    assert binding.request.target_selection.path == str(changed)
    assert control.result_labels["mutation_count"].text() == "1"
    destination = tmp_path / "pair.json"
    monkeypatch.setattr(QtWidgets.QFileDialog, "getSaveFileName", lambda *_: (str(destination), ""))
    control._export_json()
    assert load_report(destination, snapshots=binding.snapshots).report_id == binding.report.report_id
    control.variant_batch_button.click()
    assert control.variant_batch is batch
    batch.table.cellClicked.emit(0, 0)
    assert control._canonical_binding().request.target_selection.path == str(ref)
    assert len(control._report_history) == 2


def test_variant_batch_cancel_preserves_completed_pairs(batch_panel, tmp_path, monkeypatch):
    control, batch, _ = batch_panel
    targets = [make_source(tmp_path / f"v{i}.pdb") for i in range(3)]
    batch.add_paths(targets)
    wait_idle(batch)
    entered, release = Event(), Event()
    from structlens.application.report_service import ReportService
    original = ReportService.analyze

    def controlled(service, request, *args, **kwargs):
        if request.target_selection.path == str(targets[1]):
            entered.set()
            assert release.wait(10)
        return original(service, request, *args, **kwargs)

    monkeypatch.setattr(ReportService, "analyze", controlled)
    batch.run_button.click()
    try:
        end = time.monotonic() + 15
        while not entered.is_set() and time.monotonic() < end:
            QtWidgets.QApplication.processEvents()
            time.sleep(.01)
        assert entered.is_set()
        batch.cancel_button.click()
    finally:
        release.set()
    wait_idle(batch)
    assert batch.table.item(0, 2).text() == "Complete"
    assert batch.table.item(1, 2).text() == "Cancelled"
    assert batch.table.item(2, 2).text() == "Cancelled"
    assert control.model.report is None


def test_queue_keeps_original_settings_and_rejects_file_changed_after_capture(batch_panel, tmp_path):
    control, batch, reference = batch_panel
    changed = make_source(tmp_path / "changed.pdb")
    batch.add_paths([changed, reference])
    wait_idle(batch)
    make_source(changed, third="THR")
    control.mode_combo.setCurrentIndex(control.mode_combo.findData("manual"))
    batch.run_button.click()
    wait_idle(batch)
    assert "source changed after it was loaded" in batch.table.item(0, 2).text()
    assert batch.table.item(1, 2).text() == "Complete"
    assert batch.rows[1].binding.request.analysis_settings.alignment_mode.value == "sequence"


def test_manual_pairs_are_verified_separately_for_each_selected_target(batch_panel, tmp_path):
    control, batch, reference = batch_panel
    control.mode_combo.setCurrentIndex(control.mode_combo.findData("manual"))
    control.manual_edit.setPlainText("A:1 -> A:1\nA:2 -> A:2\nA:3 -> A:3")
    other_chain = make_source(tmp_path / "chain_b.pdb", ("B",))
    batch.add_paths([reference, other_chain])
    wait_idle(batch)
    control.manual_edit.clear()
    batch.run_button.click()
    wait_idle(batch)
    assert batch.table.item(0, 2).text() == "Complete"
    assert "names a residue not in the selected chains" in batch.table.item(1, 2).text()
