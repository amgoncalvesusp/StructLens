"""Compatibility namespace for integrations that probe the former Qt adapter."""

from __future__ import annotations


def legacy_site_backend() -> object:
    """Fail closed; scientific site calculations are report-service owned."""

    raise RuntimeError("Legacy site calculations are unavailable in the report-first GUI")


__all__ = ["legacy_site_backend"]
