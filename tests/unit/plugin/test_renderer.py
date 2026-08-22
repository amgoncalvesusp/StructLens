from structlens.core.models import (
    CorrespondenceStatus,
    ResidueCorrespondence,
    ResidueId,
)
from structlens.plugin.visualization.renderer import (
    HighlightFilter,
    VisualizationRenderer,
    VisualizationState,
)


def _item(
    index: int,
    status: CorrespondenceStatus,
    *,
    key: bool = False,
    outlier: bool = False,
):
    return ResidueCorrespondence(
        index,
        ResidueId("r", "1", "A", str(index), None, "ALA"),
        ResidueId("t", "1", "A", str(index), None, "ALA"),
        "A",
        "A",
        status,
        is_key_residue=key,
        is_outlier=outlier,
    )


def test_renderer_filters_mutations_and_outliers_without_changing_analysis() -> None:
    rows = (
        _item(1, CorrespondenceStatus.CONSERVED),
        _item(2, CorrespondenceStatus.SUBSTITUTION, key=True, outlier=True),
    )
    renderer = VisualizationRenderer()
    state = VisualizationState(highlight_filter=HighlightFilter.MUTATIONS)

    selected = renderer.filtered_correspondences(rows, state)

    assert [item.alignment_index for item in selected] == [2]
    assert renderer.apply_preset("Publication").highlight_filter is HighlightFilter.ALL
