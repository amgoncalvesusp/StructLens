"""Application boundary for the canonical report workflow."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Protocol

from structlens.application.dto import AnalysisReportRequest
from structlens.application.report_service import ReportService
from structlens.core.models import AnalysisSettings, ResidueId
from structlens.core.msa import MSASettings
from structlens.core.parsing import InputSelection, SourceSnapshot, StructureFormat, capture_snapshot
from structlens.core.reports import AnalysisReport
from structlens.core.sites import SiteDefinition

from .model import CanonicalReportBinding, StructLensPanelModel
from .presentation import ReportPresentation, present_report


class ReportView(Protocol):
    def show_report(self, report: AnalysisReport) -> None: ...

    def show_presentation(self, presentation: ReportPresentation) -> None: ...

    def show_status(self, status: str) -> None: ...


class ReportServiceProtocol(Protocol):
    """Minimal application boundary used by the GUI report controller."""

    def analyze(self, request: AnalysisReportRequest) -> AnalysisReport: ...


class AnalysisReportController:
    """Capture immutable inputs, execute one service, and deliver one report."""

    def __init__(self, service: ReportServiceProtocol | None = None, view: ReportView | None = None) -> None:
        self._service = service or ReportService()
        self._view = view
        self.model = StructLensPanelModel()
        self._last_request: AnalysisReportRequest | None = None
        self._pending_request: AnalysisReportRequest | None = None
        self._last_report: AnalysisReport | None = None
        self._last_presentation: ReportPresentation | None = None

    @property
    def last_request(self) -> AnalysisReportRequest | None:
        return self._last_request

    @property
    def last_report(self) -> AnalysisReport | None:
        return self._last_report

    @property
    def pending_request(self) -> AnalysisReportRequest | None:
        return self._pending_request

    @property
    def last_presentation(self) -> ReportPresentation | None:
        return self._last_presentation

    def build_request(
        self,
        reference_path: str | Path,
        target_path: str | Path,
        *,
        reference_chain_id: str | None = None,
        target_chain_id: str | None = None,
        reference_model_id: str | int | None = None,
        target_model_id: str | int | None = None,
        analysis_settings: AnalysisSettings | None = None,
        msa_settings: MSASettings | None = None,
        site_definitions: Sequence[SiteDefinition] = (),
        manual_pairs: Sequence[tuple[ResidueId, ResidueId]] = (),
    ) -> AnalysisReportRequest:
        reference_snapshot = capture_snapshot(reference_path)
        target_snapshot = capture_snapshot(target_path)
        reference_selection = self._selection(
            reference_snapshot, reference_chain_id, reference_path, reference_model_id
        )
        target_selection = self._selection(target_snapshot, target_chain_id, target_path, target_model_id)
        return self.build_request_from_snapshots(
            reference_snapshot,
            target_snapshot,
            reference_selection=reference_selection,
            target_selection=target_selection,
            analysis_settings=analysis_settings or AnalysisSettings(),
            msa_settings=msa_settings or MSASettings(),
            site_definitions=tuple(site_definitions),
            manual_pairs=tuple(manual_pairs),
        )

    def build_request_from_snapshots(
        self,
        reference_snapshot: SourceSnapshot,
        target_snapshot: SourceSnapshot,
        *,
        reference_selection: InputSelection,
        target_selection: InputSelection,
        analysis_settings: AnalysisSettings | None = None,
        msa_settings: MSASettings | None = None,
        site_definitions: Sequence[SiteDefinition] = (),
        manual_pairs: Sequence[tuple[ResidueId, ResidueId]] = (),
    ) -> AnalysisReportRequest:
        """Build a request without recapturing or reparsing immutable evidence."""

        request = AnalysisReportRequest(
            reference_snapshot,
            target_snapshot,
            reference_selection,
            target_selection,
            analysis_settings=analysis_settings or AnalysisSettings(),
            msa_settings=msa_settings or MSASettings(),
            site_definitions=tuple(site_definitions),
            manual_pairs=tuple(manual_pairs),
        )
        return request

    @staticmethod
    def verify_current_source(path: str | Path, snapshot: SourceSnapshot) -> None:
        """Fail if a mutable source no longer matches its loaded snapshot."""

        current = capture_snapshot(path)
        if current.content_id != snapshot.content_id or current.raw_sha256 != snapshot.raw_sha256:
            raise ValueError(f"source changed after it was loaded: {Path(path).name}")

    def analyze(self, request: AnalysisReportRequest) -> AnalysisReport:
        """Execute and deliver one request; retained as the headless adapter API."""

        self.model = self.model.with_busy("Building canonical analysis report…")
        report = self.execute(request)
        self.publish(report)
        return report

    def execute(self, request: AnalysisReportRequest) -> AnalysisReport:
        self._pending_request = request if isinstance(request, AnalysisReportRequest) else None
        report = self._service.analyze(request)
        if not isinstance(report, AnalysisReport):
            raise TypeError("report service must return an AnalysisReport")
        return report

    def discard_pending(self) -> None:
        """Discard an unaccepted request without touching canonical state."""

        self._pending_request = None

    def accept_report(
        self,
        report: AnalysisReport,
        request: AnalysisReportRequest | None = None,
    ) -> ReportPresentation:
        """Record a report after validating any supplied canonical request."""

        if not isinstance(report, AnalysisReport):
            raise TypeError("report must be an AnalysisReport")
        presentation = present_report(report)
        if request is not None:
            binding = CanonicalReportBinding.create(report, request, presentation)
            return self.accept_binding(binding)
        self._last_report = report
        self._last_presentation = presentation
        self._last_request = None
        self._pending_request = None
        self.model = self.model.with_report(report)
        return presentation

    def accept_binding(self, binding: CanonicalReportBinding) -> ReportPresentation:
        """Atomically accept one already verified canonical authority."""

        if not isinstance(binding, CanonicalReportBinding):
            raise TypeError("binding must be a CanonicalReportBinding")
        # Re-run validation at the acceptance boundary.  The values are
        # immutable, but this also guards hand-built lookalikes and stale
        # presentation/request combinations.
        CanonicalReportBinding.create(binding.report, binding.request, binding.presentation)
        self._last_request = binding.request
        self._pending_request = None
        self._last_report = binding.report
        self._last_presentation = binding.presentation
        self.model = self.model.with_binding(binding)
        return binding.presentation

    def publish(self, report: AnalysisReport) -> ReportPresentation:
        if not isinstance(report, AnalysisReport):
            raise TypeError("report must be an AnalysisReport")
        request = self._pending_request if isinstance(self._pending_request, AnalysisReportRequest) else None
        presentation = self.accept_report(report, request=request)
        if self._view is not None:
            self._view.show_report(report)
            show_presentation = getattr(self._view, "show_presentation", None)
            if show_presentation is not None:
                show_presentation(presentation)
            self._view.show_status(f"Analysis report complete · {report.report_id}")
        return presentation

    @staticmethod
    def _selection(
        snapshot: SourceSnapshot,
        chain_id: str | None,
        path: str | Path,
        model_id: str | int | None,
    ) -> InputSelection:
        fmt = StructureFormat(snapshot.logical_format)
        selected_model = model_id if model_id is not None else _default_model_id(snapshot, fmt)
        return InputSelection(
            snapshot.content_id,
            snapshot.display_name,
            fmt,
            selected_model,
            author_chain_ids=(chain_id,) if chain_id else (),
            path=str(path),
        )


def _default_model_id(snapshot: SourceSnapshot, fmt: StructureFormat) -> str:
    if fmt is not StructureFormat.PDB:
        return "1"
    for raw_line in snapshot.decompressed_bytes.splitlines():
        if raw_line[:6].strip().upper() == b"MODEL":
            value = raw_line[10:14].strip()
            if value:
                return value.decode("ascii", errors="strict")
    return "0"


__all__ = ["AnalysisReportController", "ReportServiceProtocol", "ReportView"]
