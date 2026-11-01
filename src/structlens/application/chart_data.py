"""Pure scientific datasets consumed by chart widgets and exporters.

No values are recalculated by a renderer.  These functions only project the
authoritative correspondence and matrix models into labelled, unit-bearing
datasets suitable for Qt, matplotlib, or XLSX export.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, overload

from structlens.core.models import (
    AnalysisResult,
    MultipleStructureAnalysis,
    PairwiseMatrix,
    ReferenceVsManyAnalysis,
    ResidueId,
    TargetAnalysis,
)
from structlens.core.models.mutation import MutationKind
from structlens.core.reports import DistanceMapSnapshot


@dataclass(frozen=True, slots=True)
class ChartSeries:
    name: str
    points: tuple[tuple[float, float | None], ...]
    labels: tuple[str, ...] = ()
    metadata: tuple[Mapping[str, Any], ...] = ()


@dataclass(frozen=True, slots=True)
class ChartDataset:
    chart_id: str
    title: str
    x_label: str
    y_label: str
    unit: str | None
    series: tuple[ChartSeries, ...]
    interpretation: str


@dataclass(frozen=True, slots=True)
class MatrixCell:
    row: str
    column: str
    value: float | None
    text: str
    status: str
    metadata: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class DistanceChartCells(Sequence[MatrixCell]):
    """Adapt the matrix lazily instead of allocating an object for every pair at startup."""

    matrix: DistanceMapSnapshot
    has_measurements: bool = field(init=False, repr=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "has_measurements", any(any(row) for row in self.matrix.valid_mask))

    def __len__(self) -> int:
        return len(self.matrix.reference_positions) ** 2

    @overload
    def __getitem__(self, index: int) -> MatrixCell: ...
    @overload
    def __getitem__(self, index: slice) -> tuple[MatrixCell, ...]: ...

    def __getitem__(self, index: int | slice) -> MatrixCell | tuple[MatrixCell, ...]:
        if isinstance(index, slice):
            return tuple(self[i] for i in range(*index.indices(len(self))))
        labels = self.matrix.reference_positions
        if not -len(self) <= index < len(self):
            raise IndexError(index)
        r, c = divmod(index % len(self), len(labels))
        value = self.matrix.delta_angstrom[r][c] if self.matrix.valid_mask[r][c] else None
        return MatrixCell(labels[r], labels[c], value, "—" if value is None else f"{value:+.3f} Å",
                          "unavailable" if value is None else "measured")


@dataclass(frozen=True, slots=True)
class MatrixDataset:
    chart_id: str
    title: str
    row_label: str
    column_label: str
    cells: tuple[MatrixCell, ...] | DistanceChartCells
    interpretation: str


def structural_deviation_profile(
    analysis: AnalysisResult | ReferenceVsManyAnalysis,
    *,
    metric: str = "ca_displacement_angstrom",
    target_ids: Sequence[str] | None = None,
) -> ChartDataset:
    """Return one residue-position series per selected target."""

    results: dict[str, AnalysisResult | TargetAnalysis]
    if isinstance(analysis, AnalysisResult):
        results = {analysis.target_id: analysis}
    else:
        results = {
            target_id: target
            for target_id, target in analysis.targets.items()
            if target_ids is None or target_id in target_ids
        }
    labels = {
        "ca_displacement_angstrom": "Cα displacement (Å)",
        "backbone_rmsd_angstrom": "Backbone RMSD (Å)",
        "sidechain_rmsd_angstrom": "Side-chain RMSD (Å)",
        "all_heavy_atom_rmsd_angstrom": "Local RMSD (Å)",
    }
    if metric not in labels:
        raise ValueError(f"Unsupported structural deviation metric: {metric}")
    series: list[ChartSeries] = []
    for target_id, result in results.items():
        correspondences = result.correspondences if isinstance(result, AnalysisResult) else result.correspondence
        points: list[tuple[float, float | None]] = []
        point_labels: list[str] = []
        metadata: list[Mapping[str, Any]] = []
        for item in correspondences:
            position = _residue_position(item.reference)
            if position is None:
                continue
            points.append((float(position), getattr(item, metric)))
            point_labels.append(_residue_label(item.reference))
            metadata.append({"alignment_index": item.alignment_index, "target_id": target_id})
        series.append(ChartSeries(target_id, tuple(points), tuple(point_labels), tuple(metadata)))
    return ChartDataset(
        "structural_deviation_profile",
        "Structural deviation profile",
        "Reference residue position",
        labels[metric],
        "Å",
        tuple(series),
        "Each point is an explicit mapped residue metric; missing values remain empty.",
    )


def mutation_conservation_matrix(
    analysis: ReferenceVsManyAnalysis | AnalysisResult,
) -> MatrixDataset:
    """Project mutation status into a reference-position matrix."""

    if isinstance(analysis, AnalysisResult):
        rows = {analysis.target_id: (analysis.correspondences, analysis.mutations)}
    else:
        rows = {name: (target.correspondence, target.mutations) for name, target in analysis.targets.items()}
    scores = {MutationKind.CONSERVED: 1.0, MutationKind.SUBSTITUTION: 0.0}
    cells: list[MatrixCell] = []
    for row, (correspondences, mutations) in rows.items():
        by_index = {event.alignment_index: event for event in mutations}
        for item in correspondences:
            if item.reference is None:
                continue
            event = by_index.get(item.alignment_index)
            cells.append(
                MatrixCell(
                    row,
                    _residue_label(item.reference),
                    scores.get(event.kind) if event else None,
                    event.canonical_notation if event else "—",
                    event.kind.value if event else "unavailable",
                    {"alignment_index": item.alignment_index, "target_residue": _residue_label(item.target)},
                )
            )
    return MatrixDataset(
        "mutation_conservation_matrix",
        "Mutation / conservation matrix",
        "Structure",
        "Reference-aligned position",
        tuple(cells),
        "Blue marks conserved amino acids (1), orange substitutions (0), and grey unscored positions. "
        "Columns use reference positions; insertions without a reference position appear in the residue table.",
    )


def pairwise_similarity_heatmap(matrix: PairwiseMatrix) -> MatrixDataset:
    """Project a symmetric PairwiseMatrix without duplicating calculations."""

    cells: list[MatrixCell] = []
    for row in matrix.structure_ids:
        for column in matrix.structure_ids:
            value = matrix.value(row, column)
            cells.append(
                MatrixCell(
                    row,
                    column,
                    value,
                    "—" if value is None else f"{value:.3f}",
                    "diagonal" if row == column else "value",
                )
            )
    return MatrixDataset(
        "pairwise_similarity_heatmap",
        f"Pairwise {matrix.metric_name} heatmap",
        "Structure",
        "Structure",
        tuple(cells),
        "Off-diagonal values are calculated once and mirrored visually.",
    )


def sequence_structure_relationship(
    analyses: Iterable[AnalysisResult | ReferenceVsManyAnalysis],
    *,
    y_metric: str = "tm_score",
) -> ChartDataset:
    """Return sequence identity versus a structural metric scatter dataset."""

    expanded: list[AnalysisResult] = []
    for analysis in analyses:
        if isinstance(analysis, AnalysisResult):
            expanded.append(analysis)
        else:
            expanded.extend(_result_like_target(target) for target in analysis.targets.values())
    allowed = {"tm_score", "strict_rmsd_angstrom", "refined_rmsd_angstrom"}
    if y_metric not in allowed:
        raise ValueError(f"Unsupported sequence-structure metric: {y_metric}")
    points = tuple(
        (result.sequence_identity * 100.0, getattr(result, y_metric))
        for result in expanded
        if getattr(result, y_metric) is not None
    )
    labels = tuple(result.target_id for result in expanded if getattr(result, y_metric) is not None)
    return ChartDataset(
        "sequence_structure_relationship",
        "Sequence–structure relationship",
        "Sequence identity (%)",
        y_metric.replace("_", " ").title() + (" (Å)" if "rmsd" in y_metric else ""),
        "Å" if "rmsd" in y_metric else "score",
        (ChartSeries("pairs", points, labels),),
        "Each point is an analysis result; selecting it identifies the stored target pair.",
    )


def structural_conservation_profile(analysis: MultipleStructureAnalysis) -> ChartDataset:
    """Return Cα positional variability and coverage for each aligned position."""

    variability = tuple(
        (_residue_position(position.reference_residue) or position.alignment_index + 1,
         position.ca_positional_variability_angstrom)
        for position in analysis.aligned_positions
    )
    coverage = tuple(
        (_residue_position(position.reference_residue) or position.alignment_index + 1,
         position.coverage * 100.0)
        for position in analysis.aligned_positions
    )
    return ChartDataset(
        "structural_conservation_profile",
        "Structural conservation profile",
        "Reference-aligned residue position",
        "Cα positional variability (Å)",
        "Å",
        (ChartSeries("variability", variability), ChartSeries("coverage (%)", coverage)),
        "Coverage is shown separately so sparsely mapped positions are not treated as equally reliable.",
    )


def key_residue_comparison(
    analysis: ReferenceVsManyAnalysis | AnalysisResult,
    key_residues: Sequence[ResidueId],
    *,
    metric: str = "ca_displacement_angstrom",
) -> ChartDataset:
    """Return one series per target for explicitly selected reference residues."""

    results = (
        {analysis.target_id: analysis.correspondences}
        if isinstance(analysis, AnalysisResult)
        else {target_id: target.correspondence for target_id, target in analysis.targets.items()}
    )
    positions = tuple(float(_residue_position(residue) or index + 1) for index, residue in enumerate(key_residues))
    series = []
    for target_id, correspondences in results.items():
        by_residue = {item.reference: item for item in correspondences}
        values = tuple(getattr(by_residue[residue], metric) if residue in by_residue else None for residue in key_residues)
        series.append(ChartSeries(target_id, tuple(zip(positions, values, strict=True)), tuple(_residue_label(item) for item in key_residues)))
    return ChartDataset(
        "key_residue_comparison",
        "Key-residue comparison",
        "Key reference residue",
        metric.replace("_", " ").title() + " (Å)",
        "Å",
        tuple(series),
        "Only explicitly selected key residues are included; missing mappings remain empty.",
    )


def _residue_position(residue: ResidueId | None) -> int | None:
    if residue is None:
        return None
    try:
        return int(residue.auth_seq_id)
    except ValueError:
        return None


def _residue_label(residue: ResidueId | None) -> str:
    if residue is None:
        return "—"
    return f"{residue.chain_id}:{residue.auth_seq_id}{residue.insertion_code or ''}"


def _result_like_target(target: Any) -> AnalysisResult:
    """Build a lightweight result for scatter labels without recomputation."""

    return AnalysisResult(
        reference_id="reference",
        target_id=target.target_id,
        correspondences=target.correspondence,
        mutations=target.mutations,
        sequence_identity=target.sequence_metrics.identity,
        sequence_similarity=target.sequence_metrics.similarity,
        sequence_coverage=target.sequence_metrics.coverage,
        alignment_decision="imported multi-analysis",
        strict_rmsd_angstrom=(
            target.structural_metrics.strict_ca_rmsd_angstrom
            if target.structural_metrics is not None
            else None
        ),
    )


__all__ = [
    "ChartDataset",
    "ChartSeries",
    "MatrixCell",
    "MatrixDataset",
    "key_residue_comparison",
    "mutation_conservation_matrix",
    "pairwise_similarity_heatmap",
    "sequence_structure_relationship",
    "structural_conservation_profile",
    "structural_deviation_profile",
]
