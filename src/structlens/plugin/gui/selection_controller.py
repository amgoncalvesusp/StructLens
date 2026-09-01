"""Immutable correspondence and evidence focus for report-driven views."""

from __future__ import annotations

from dataclasses import dataclass, field, replace

from structlens.core.evidence import EvidenceCard
from structlens.core.models import AnalysisSelection, ResidueId
from structlens.core.reports import AnalysisReport, CorrespondenceSnapshot


@dataclass(frozen=True, slots=True)
class AnalysisSelectionController:
    """Bind UI focus to one report and return replacements for every change."""

    report: AnalysisReport | None = None
    selection: AnalysisSelection = field(default_factory=AnalysisSelection)

    @property
    def report_id(self) -> str | None:
        return self.report.report_id if self.report is not None else None

    def bind_report(self, report: AnalysisReport) -> AnalysisSelectionController:
        if not isinstance(report, AnalysisReport):
            raise TypeError("report must be an AnalysisReport")
        return replace(self, report=report, selection=AnalysisSelection())

    def clear_selection(self) -> AnalysisSelectionController:
        return replace(self, selection=AnalysisSelection())

    def select_alignment_index(self, alignment_index: int) -> AnalysisSelectionController:
        _validate_alignment_index(alignment_index)
        correspondence = self._correspondence(alignment_index)
        selection = AnalysisSelection(
            reference_residue=correspondence.reference,
            target_id=self.report.analysis.target_id if self.report and self.report.analysis else None,
            target_residue=correspondence.target,
            pair=(self.report.analysis.reference_id, self.report.analysis.target_id)
            if self.report is not None and self.report.analysis is not None
            else None,
        )
        return replace(self, selection=selection)

    def select_evidence_card(self, reference_residue: ResidueId) -> AnalysisSelectionController:
        if not isinstance(reference_residue, ResidueId):
            raise TypeError("reference_residue must be a ResidueId")
        card = next(
            (item for item in self._report_or_error().evidence_cards if item.residue_id == reference_residue), None
        )
        if card is None:
            raise KeyError(f"evidence card not found for {reference_residue}")
        return replace(self, selection=replace(self.selection, reference_residue=reference_residue))

    def lookup_correspondence(self, alignment_index: int | None = None) -> CorrespondenceSnapshot:
        if alignment_index is None:
            alignment_index = self._selected_alignment_index()
        _validate_alignment_index(alignment_index)
        return self._correspondence(alignment_index)

    def lookup_evidence_card(self, alignment_index: int | None = None) -> EvidenceCard:
        report = self._report_or_error()
        if alignment_index is None:
            alignment_index = self._selected_alignment_index()
        _validate_alignment_index(alignment_index)
        correspondence = self._correspondence(alignment_index)
        for card in report.evidence_cards:
            if card.residue_id == correspondence.reference:
                return card
        raise KeyError(f"evidence card not found for alignment index {alignment_index}")

    def _report_or_error(self) -> AnalysisReport:
        if self.report is None:
            raise ValueError("a report must be bound before selecting evidence")
        return self.report

    def _correspondence(self, alignment_index: int) -> CorrespondenceSnapshot:
        report = self._report_or_error()
        if report.analysis is None:
            raise ValueError("the bound report has no analysis correspondences")
        for item in report.analysis.correspondences:
            if item.alignment_index == alignment_index:
                return item
        raise IndexError(f"alignment index {alignment_index} is not present in the bound report")

    def _selected_alignment_index(self) -> int:
        report = self._report_or_error()
        reference = self.selection.reference_residue
        target = self.selection.target_residue
        if reference is None:
            raise ValueError("no correspondence is selected")
        if report.analysis is None:
            raise ValueError("the bound report has no analysis correspondences")
        for item in report.analysis.correspondences:
            if item.reference == reference and (target is None or item.target == target):
                return item.alignment_index
        raise ValueError("selected correspondence does not belong to the bound report")


__all__ = ["AnalysisSelectionController"]


def _validate_alignment_index(value: object) -> None:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError("alignment index must be an integer")
