"""Version metadata shown by Help → About → Scientific Backends."""

from __future__ import annotations

from importlib import resources

from structlens import __version__


def backend_versions() -> dict[str, str]:
    values = {"StructLens": __version__, "bundle_schema": "3", "Python": "runtime"}
    for name, package in (("MUSCLE", "muscle"), ("US-align", "usalign"), ("FreeSASA", "freesasa")):
        try:
            values[name] = resources.files("structlens.resources").joinpath(package, "VERSION").read_text(encoding="utf-8").strip()
        except (FileNotFoundError, ModuleNotFoundError):
            values[name] = "unavailable"
    return values


__all__ = ["backend_versions"]
