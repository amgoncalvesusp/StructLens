"""Regression tests for the dependencies required by structural analysis."""

from __future__ import annotations

import re
import tomllib
from pathlib import Path

PROJECT_FILE = Path(__file__).parents[2] / "pyproject.toml"


def _dependency_spec(name: str, group: list[str]) -> str:
    prefix = f"{name}["
    for dependency in group:
        if dependency == name or dependency.startswith(f"{name}>") or dependency.startswith(f"{name}="):
            return dependency
        if dependency.startswith(prefix):
            return dependency
    raise AssertionError(f"{name!r} is not declared")


def _minimum_version(spec: str, name: str) -> tuple[int, ...]:
    match = re.search(rf"^{re.escape(name)}>=([0-9]+(?:\.[0-9]+)*)", spec)
    assert match, f"{name} must declare a lower bound: {spec!r}"
    return tuple(int(part) for part in match.group(1).split("."))


def test_structural_runtime_and_test_dependencies_are_declared() -> None:
    metadata = tomllib.loads(PROJECT_FILE.read_text(encoding="utf-8"))
    project = metadata["project"]
    runtime = project["dependencies"]
    development = project["optional-dependencies"]["dev"]

    scipy_spec = _dependency_spec("scipy", runtime)
    hypothesis_spec = _dependency_spec("hypothesis", development)

    assert _minimum_version(scipy_spec, "scipy") >= (1, 10)
    assert _minimum_version(hypothesis_spec, "hypothesis") >= (6, 90)


REPOSITORY = PROJECT_FILE.parent
RELEASE_WORKFLOW = REPOSITORY / ".github" / "workflows" / "release.yml"
LICENSING_POLICY = REPOSITORY / "docs" / "distribution-licensing.md"
MUSCLE_LICENSE = REPOSITORY / "licenses" / "MUSCLE-LICENSE.txt"


def _optional_dependencies() -> dict[str, list[str]]:
    metadata = tomllib.loads(PROJECT_FILE.read_text(encoding="utf-8"))
    groups: dict[str, list[str]] = metadata["project"]["optional-dependencies"]
    return groups


def test_default_gui_extra_is_the_lgpl_binding() -> None:
    """Official artifacts freeze PySide6; PyQt5 must never be the default path."""

    groups = _optional_dependencies()
    default = " ".join(groups["gui"])
    assert "PySide6" in default
    assert "PyQt5" not in default


def test_pyqt5_remains_an_explicitly_named_gpl_extra() -> None:
    groups = _optional_dependencies()
    assert "gui-pyqt5" in groups, "PyQt5 must stay available as an opt-in extra"
    assert "PyQt5" in " ".join(groups["gui-pyqt5"])


def test_release_workflow_does_not_freeze_pyqt5() -> None:
    """Freezing PyQt5 into a distributed binary triggers its GPL obligations."""

    workflow = RELEASE_WORKFLOW.read_text(encoding="utf-8")
    assert "--hidden-import PyQt5" not in workflow
    assert "--exclude-module PySide6" not in workflow


def test_licensing_policy_enumerates_every_redistributed_component() -> None:
    assert LICENSING_POLICY.is_file(), "docs/distribution-licensing.md must exist"
    policy = LICENSING_POLICY.read_text(encoding="utf-8")
    for component in ("StructLens", "PySide6", "PyQt5", "MUSCLE", "US-align", "FreeSASA"):
        assert component in policy, f"{component} is not covered by the licensing policy"
    for obligation in ("LGPL", "GPL-3", "corresponding source"):
        assert obligation in policy, f"the policy omits the {obligation} obligation"


def test_muscle_notice_carries_the_full_gpl_text_and_source_route() -> None:
    """A redistributed GPL-3 binary needs its licence text, not a promise of one."""

    notice = MUSCLE_LICENSE.read_text(encoding="utf-8")
    assert "GNU GENERAL PUBLIC LICENSE" in notice
    assert "Version 3" in notice
    assert len(notice.splitlines()) > 100, "the notice must be the complete GPL-3 text"


def test_backend_fetcher_retrieves_the_muscle_licence_and_source() -> None:
    """Redistributing the GPL-3 binary means fetching its licence and source too."""

    script = (REPOSITORY / "scripts" / "fetch_backends.py").read_text(encoding="utf-8")
    assert "MUSCLE_LICENSE_URL" in script, "the GPL-3 licence text must be fetched, not assumed"
    assert "MUSCLE_SOURCE_URL" in script, "the corresponding source archive must be fetched"
    assert "muscle-5.3-source.tar.gz" in script


def test_release_verifier_requires_the_muscle_compliance_files() -> None:
    verifier = (REPOSITORY / "scripts" / "verify_release_backends.py").read_text(encoding="utf-8")
    assert "LICENSE" in verifier
    assert "source" in verifier.lower()
