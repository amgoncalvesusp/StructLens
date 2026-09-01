"""Operate-mode Qt panel for the optional StructLens PyMOL plugin.

Impeccable direction contract — Evidence Bench
THESIS: make the residue correspondence table the primary workspace, refusing a
decorative dashboard of unexplained metrics.
OWN-WORLD: host-compatible graphite surfaces, cool blue action states, compact
system sans controls, and dense tabular data with explicit Å/fraction units.
STORY: the user loads two sources, chooses an alignment policy, runs one
reproducible comparison, then moves from evidence to a reversible PyMOL view.
FIRST VIEWPORT: a narrow workflow rail anchors a two-column Project page; the
header keeps status and Compare visible while the right canvas holds source
controls, chain choices, and the next action.
FORM: Operate-mode split workspace, chosen to keep eight scientific stages visible
without burying the task in a tab strip or card mosaic.
"""

from __future__ import annotations

from .model import GUI_SECTIONS, SCIENTIFIC_SECTIONS, WORKFLOW_HELP, StructLensPanelModel


def build_qt_panel(
    parent: object | None = None,
    *,
    command: object | None = None,
) -> object:
    """Build the full panel when a supported Qt binding is present."""

    from .qt_compat import load_qt

    qt = load_qt()
    if qt is None:
        raise RuntimeError(
            "StructLens GUI requires PySide6 or PyQt5; install structlens[gui] "
            "or launch it inside a PyMOL environment that provides Qt"
        )
    from .qt_panel import build_panel

    return build_panel(qt, parent=parent, command=command)


__all__ = [
    "GUI_SECTIONS",
    "SCIENTIFIC_SECTIONS",
    "WORKFLOW_HELP",
    "AnalysisReportController",
    "StructLensPanelModel",
    "build_qt_panel",
]

# Keep the public import location stable while the report boundary lives in a
# focused module.  The import is intentionally last so that the model above is
# fully defined before the controller references it.
from .report_controller import AnalysisReportController  # noqa: E402
