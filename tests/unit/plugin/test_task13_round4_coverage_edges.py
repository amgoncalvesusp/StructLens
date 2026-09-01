"""Remaining GUI edge paths kept separate to keep test modules reviewable."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pytest
from PySide6.QtWidgets import QApplication

from structlens.application.chart_data import ChartDataset, ChartSeries, MatrixCell, MatrixDataset
from structlens.plugin.gui.main_panel import build_qt_panel
from structlens.plugin.gui.qt_reports import legacy_analysis

from ._task13_round4_fixtures import rich_report


@pytest.fixture(scope="module")
def application() -> QApplication:
    app = QApplication.instance() or QApplication([])
    yield app
    app.quit()


def test_site_input_and_refresh_status_error_paths(application, monkeypatch) -> None:
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    errors: list[str] = []
    monkeypatch.setattr(controller, "_show_error", errors.append)
    unavailable: list[str] = []
    site_unavailable = controller._site_unavailable
    monkeypatch.setattr(controller, "_site_unavailable", unavailable.append)
    controller.site_mode_combo = SimpleNamespace(currentData=lambda: None)
    controller._define_site_from_controls()
    controller.site_mode_combo = SimpleNamespace(currentData=lambda: "unsupported")
    controller._define_site_from_controls()
    controller.site_mode_combo = SimpleNamespace(currentData=lambda: "key_residues")
    controller.site_residues_edit.clear()
    controller._define_site_from_controls()
    assert len(unavailable) == 3

    monkeypatch.setattr(controller, "_site_unavailable", site_unavailable)
    controller._site_unavailable("explicit unavailable")
    assert errors[-1] == "explicit unavailable"
    controller.site_mode_combo = SimpleNamespace(currentData=lambda: "key_residues")
    controller.site_residues_edit.setText("A:999")
    fixture = Path(__file__).parents[2] / "fixtures" / "parsing" / "numbering_altloc.pdb"
    assert controller._load_source("reference", fixture)
    controller._define_site_from_controls()
    assert "positions" in controller.site_status_label.text().lower()

    class Launcher:
        def __init__(self, *_args, **_kwargs):
            pass

        def locate(self):
            return Path("pymol.exe")

    monkeypatch.setattr("structlens.plugin.gui.qt_exports_mixin.PyMOLLauncher", Launcher)
    controller._refresh_pymol_status()
    assert "pymol.exe" in controller.pymol_status.text()
    controller._set_combo_value(controller.mode_combo, "does-not-exist")
    controller.close()
    panel.deleteLater()


def test_visualization_mixin_handles_renderers_legacy_and_all_legend_modes(
    application, tmp_path: Path, monkeypatch
) -> None:
    request, report = rich_report(tmp_path)
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    assert controller._load_source("reference", Path(request.reference_selection.path))
    assert controller._load_source("target", Path(request.target_selection.path))
    controller._report_request = request
    controller._analysis_finished(report)

    controller._apply_preset("")
    controller._update_legend()
    controller.color_combo.setCurrentText("Cα displacement")
    controller._update_legend()
    controller.color_combo.setCurrentText("Backbone RMSD")
    controller._update_legend()
    controller.color_combo.setCurrentText("Mutation status")
    controller._update_legend()

    controller.set_chart_datasets(
        {
            "round-four-matrix": MatrixDataset(
                "round-four-matrix",
                "Round-four matrix",
                "row",
                "column",
                (MatrixCell("r", "c", 1.0, "1", "value"),),
                "Stored matrix values.",
            ),
            "round-four-chart": ChartDataset(
                "round-four-chart",
                "Round-four chart",
                "x",
                "y",
                "Å",
                (
                    ChartSeries("empty", ((1.0, None),)),
                    ChartSeries("observed", ((1.0, 1.0),)),
                ),
                "Stored chart values.",
            ),
        }
    )
    controller.chart_combo.setCurrentText("round-four-matrix")
    controller._render_selected_chart()
    controller.chart_combo.setCurrentText("round-four-chart")
    controller._render_selected_chart()
    monkeypatch.setattr(
        "structlens.plugin.gui.qt_visualization_mixin.load_chart_classes",
        lambda: None,
    )
    controller._render_selected_chart()
    controller._render_sequence_result(legacy_analysis(report))
    controller.close()
    panel.deleteLater()


def test_visualization_mixin_chart_rendering_and_dataset_fallbacks(application, tmp_path: Path) -> None:
    _, report = rich_report(tmp_path)
    panel = build_qt_panel(command=None)
    controller = panel._structlens_controller
    legacy = legacy_analysis(report)
    assert legacy is not None
    correspondences = tuple(
        replace(item, ca_displacement_angstrom=1.25, backbone_rmsd_angstrom=2.5) for item in legacy.correspondences
    )
    legacy = replace(legacy, correspondences=correspondences)
    controller.model = controller.model.with_analysis(legacy)

    controller._set_combo_value(controller.color_combo, "ca_displacement")
    controller._update_legend()
    assert "Cα displacement" in controller.legend_label.text()
    controller._set_combo_value(controller.color_combo, "backbone_rmsd")
    controller._update_legend()
    assert "Backbone RMSD" in controller.legend_label.text()

    chart = ChartDataset(
        "round-four-chart",
        "Round-four chart",
        "x",
        "y",
        "Å",
        (
            ChartSeries("observed", ((1.0, 1.0), (2.0, None))),
            ChartSeries("second", ((1.0, 2.0),)),
        ),
        "Stored chart values.",
    )
    controller._render_dataset(
        chart,
        controller.chart_preview_layout,
        controller.chart_preview_status,
        canvas_attribute="_chart_canvas",
        unavailable_message="not rendered",
    )
    matrix = MatrixDataset(
        "round-four-matrix",
        "Round-four matrix",
        "row",
        "column",
        (
            MatrixCell("r1", "c1", 1.0, "1", "value"),
            MatrixCell("r1", "c2", None, "—", "missing"),
            MatrixCell("r2", "c1", -1.0, "-1", "value"),
        ),
        "Stored matrix values.",
    )
    controller._render_dataset(
        matrix,
        controller.chart_preview_layout,
        controller.chart_preview_status,
        canvas_attribute="_chart_canvas",
        unavailable_message="not rendered",
    )

    controller._update_chart_explanation("Structural deviation profile")
    controller.chart_combo.setCurrentText("Structural deviation profile")
    assert controller._selected_chart_dataset(legacy) is not None
    controller._chart_datasets = {}
    controller.chart_combo.setCurrentText("MSA conservation profile")
    fallback = ChartDataset(
        "msa",
        "MSA",
        "column",
        "score",
        "fraction",
        (ChartSeries("alignment", ((1.0, 1.0),)),),
        "Stored alignment values.",
    )
    controller._msa_chart_dataset = fallback
    assert controller._selected_chart_dataset(legacy) is fallback
    controller.set_chart_datasets({})
    controller.close()
    panel.deleteLater()
