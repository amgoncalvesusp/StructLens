from structlens.core.models import (
    AnalysisResult,
    CorrespondenceStatus,
    ResidueCorrespondence,
    ResidueId,
)
from structlens.integrations.pymol.adapter import PyMOLAdapter
from structlens.plugin.visualization.renderer import (
    ColorMode,
    HighlightFilter,
    VisualizationState,
)


class FakeCommand:
    def __init__(self) -> None:
        self.calls: list[tuple[str, tuple[object, ...]]] = []

    def __getattr__(self, name: str):
        def record(*args: object) -> None:
            self.calls.append((name, args))

        return record


def _result() -> AnalysisResult:
    reference = ResidueId("reference", "1", "A", "42", None, "ALA")
    target = ResidueId("target", "1", "A", "57", None, "ALA")
    item = ResidueCorrespondence(
        0,
        reference,
        target,
        "A",
        "A",
        CorrespondenceStatus.CONSERVED,
        ca_displacement_angstrom=1.4,
    )
    return AnalysisResult(
        reference_id="reference",
        target_id="target",
        correspondences=(item,),
        mutations=(),
        sequence_identity=1.0,
        sequence_coverage=1.0,
        alignment_decision="sequence-guided",
        mapped_residue_count=1,
    )


def test_adapter_applies_namespaced_view_and_resets_only_owned_selections() -> None:
    command = FakeCommand()
    adapter = PyMOLAdapter(command, project_id="GUI")

    adapter.apply(
        _result(),
        state=VisualizationState(
            highlight_filter=HighlightFilter.ALL,
            color_mode=ColorMode.CA_DISPLACEMENT,
            show_labels=True,
        ),
    )

    selected_names = [
        args[0]
        for method, args in command.calls
        if method == "select"
    ]
    created_names = [
        args[0]
        for method, args in command.calls
        if method == "create"
    ]
    assert created_names
    assert selected_names
    assert all(str(name).startswith("structlens_gui_") for name in selected_names)
    assert any(method == "label" for method, _ in command.calls)

    adapter.reset()

    deleted_names = [args[0] for method, args in command.calls if method == "delete"]
    assert set(selected_names).issubset(set(deleted_names))
    assert set(created_names).issubset(set(deleted_names))
    assert deleted_names
