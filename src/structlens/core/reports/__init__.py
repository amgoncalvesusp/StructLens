"""Canonical report contracts for StructLens scientific outputs."""

from .models import (
    AnalysisReport,
    AnalysisSnapshot,
    CorrespondenceSnapshot,
    DistanceMapSnapshot,
    InputQualityBundle,
    SectionAvailability,
)
from .pockets import (
    PocketDetectionSnapshot,
    PocketLiningSnapshot,
    PocketReportSnapshot,
    PocketVolumeSnapshot,
)
from .schema import (
    SUPPORTED_REPORT_SCHEMA_VERSION,
    AnalysisReportSchemaError,
    load_analysis_report_schema,
    validate_analysis_report_payload,
    validate_report_schema,
)

__all__ = [
    "AnalysisReport",
    "AnalysisSnapshot",
    "CorrespondenceSnapshot",
    "DistanceMapSnapshot",
    "InputQualityBundle",
    "SectionAvailability",
    "PocketDetectionSnapshot",
    "PocketLiningSnapshot",
    "PocketReportSnapshot",
    "PocketVolumeSnapshot",
    "AnalysisReportSchemaError",
    "SUPPORTED_REPORT_SCHEMA_VERSION",
    "load_analysis_report_schema",
    "validate_analysis_report_payload",
    "validate_report_schema",
]
