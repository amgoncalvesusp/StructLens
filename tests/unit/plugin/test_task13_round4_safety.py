"""Round-four RED tests for deep immutability and safe chart output."""

from __future__ import annotations

from pathlib import Path

import pytest
from matplotlib.figure import Figure
from openpyxl import Workbook
from PySide6.QtWidgets import QApplication

from structlens.application.chart_data import ChartDataset, ChartSeries
from structlens.application.chart_export import export_chart_image, export_chart_xlsx
from structlens.plugin.gui import presentation as presentation_module
from structlens.plugin.gui.presentation import present_report

from ._task13_round4_fixtures import gapped_card, rich_report


@pytest.fixture(scope="module")
def application() -> QApplication:
    app = QApplication.instance() or QApplication([])
    yield app
    app.quit()


def _dataset() -> ChartDataset:
    return ChartDataset(
        "round4",
        "Round-four chart",
        "position",
        "value",
        "Å",
        (ChartSeries("observed", ((1.0, 2.5),), ("Z:900",)),),
        "Descriptive fixture chart.",
    )


def test_report_presentation_freezes_nested_bundle_and_chart_payloads(tmp_path: Path) -> None:
    """A frozen presentation must not leak mutable nested export structures."""

    _, report = rich_report(tmp_path)
    presentation = present_report(report)

    with pytest.raises((TypeError, AttributeError)):
        presentation.bundle_payloads["evidence"]["cards"][0]["target_id"] = "tampered"
    with pytest.raises((TypeError, AttributeError)):
        presentation.chart_datasets["distance_difference_heatmap"]["data"]["delta_angstrom"][0][0] = 99.0

    first = present_report(report)
    second = present_report(report)
    assert first == second
    assert first.bundle_payloads["evidence"] == second.bundle_payloads["evidence"]
    assert first.chart_datasets["distance_difference_heatmap"] == second.chart_datasets["distance_difference_heatmap"]


def test_thaw_json_returns_an_independent_recursive_copy(tmp_path: Path) -> None:
    """Compatibility consumers receive a fresh mutable copy at the boundary."""

    _, report = rich_report(tmp_path)
    presentation = present_report(report)
    thaw = getattr(presentation_module, "thaw_json", None)
    assert callable(thaw), "presentation must expose explicit recursive thaw_json conversion"
    copy = thaw(presentation.bundle_payloads)
    copy["evidence"]["cards"].append({"target_id": "new"})
    copy["evidence"]["cards"][0]["target_id"] = "tampered"

    assert len(presentation.bundle_payloads["evidence"]["cards"]) != len(copy["evidence"]["cards"])
    assert presentation.bundle_payloads["evidence"]["cards"][0]["target_id"] != "tampered"


def test_compatibility_bundle_setter_deep_copies_nested_values(application, tmp_path: Path) -> None:
    """Legacy GUI setters cannot mutate the canonical export payload through aliases."""

    panel = __import__("structlens.plugin.gui.main_panel", fromlist=["build_qt_panel"]).build_qt_panel(command=None)
    controller = panel._structlens_controller
    payload = {"cards": [{"target_id": "target", "nested": [1, 2]}]}
    controller.set_v03_bundle_payloads(evidence=payload)
    payload["cards"][0]["target_id"] = "tampered"
    payload["cards"][0]["nested"].append(3)
    first = controller._v03_bundle_kwargs()["evidence"]
    assert first["cards"][0]["target_id"] == "target"
    assert first["cards"][0]["nested"] == [1, 2]

    first["cards"][0]["target_id"] = "mutated getter"
    second = controller._v03_bundle_kwargs()["evidence"]
    assert second["cards"][0]["target_id"] == "target"
    controller.close()
    panel.deleteLater()


def test_gapped_evidence_card_row_uses_authoritative_alignment_index(tmp_path: Path) -> None:
    """A gap/insertion cannot relabel the card with the residue sequence index."""

    _, base = rich_report(tmp_path)
    report = gapped_card(base)
    presentation = present_report(report)
    rows = presentation.sections.evidence_cards.rows

    assert rows
    assert rows[0]["alignment_index"] == "17"


def test_chart_xlsx_failure_leaves_an_existing_destination_unchanged(tmp_path: Path, monkeypatch) -> None:
    """Chart serialization must render before atomically replacing the destination."""

    output = tmp_path / "chart.xlsx"
    output.write_bytes(b"previous-valid-output")

    def partial_save(_workbook: Workbook, target: object) -> None:
        if hasattr(target, "write"):
            target.write(b"partial-output")  # type: ignore[attr-defined]
        else:
            Path(target).write_bytes(b"partial-output")
        raise OSError("simulated chart write failure")

    monkeypatch.setattr(Workbook, "save", partial_save)
    with pytest.raises((OSError, RuntimeError, ValueError)):
        export_chart_xlsx(_dataset(), output)

    assert output.read_bytes() == b"previous-valid-output"


def test_chart_image_failure_leaves_an_existing_destination_unchanged(tmp_path: Path, monkeypatch) -> None:
    """A failed image render cannot publish a partial file over valid output."""

    output = tmp_path / "chart.png"
    output.write_bytes(b"previous-valid-output")

    def partial_savefig(_figure: Figure, target: object, **_kwargs: object) -> None:
        if hasattr(target, "write"):
            target.write(b"partial-output")  # type: ignore[attr-defined]
        else:
            Path(target).write_bytes(b"partial-output")
        raise OSError("simulated image write failure")

    monkeypatch.setattr(Figure, "savefig", partial_savefig)
    with pytest.raises((OSError, RuntimeError, ValueError)):
        export_chart_image(_dataset(), output, dpi=300)

    assert output.read_bytes() == b"previous-valid-output"
