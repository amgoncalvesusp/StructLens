from structlens.plugin.gui.main_panel import (
    GUI_SECTIONS,
    WORKFLOW_HELP,
    StructLensPanelModel,
)


def test_gui_defines_six_scientific_sections_and_context_help() -> None:
    assert GUI_SECTIONS == (
        "Project",
        "Alignment",
        "Mutations",
        "Residues",
        "Visualization",
        "Results",
    )
    for key in (
        "Auto",
        "Sequence",
        "Structure",
        "Manual",
        "strict RMSD",
        "Cα displacement",
    ):
        assert key in WORKFLOW_HELP
        assert "what" in WORKFLOW_HELP[key].lower()


def test_panel_model_keeps_ui_state_separate_from_analysis_state() -> None:
    model = StructLensPanelModel()
    updated = model.with_status("Ready")
    assert updated.status == "Ready"
    assert updated.analysis is None
    assert model.status == "Choose a reference structure and a target to begin."


def test_panel_status_transition_leaves_cancelled_model_idle() -> None:
    model = StructLensPanelModel().with_busy("Comparing structures…")
    cancelled = model.with_status("Comparison cancelled.")

    assert model.busy is True
    assert cancelled.busy is False
    assert cancelled.status == "Comparison cancelled."
