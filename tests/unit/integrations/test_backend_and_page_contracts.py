"""Cover the small backend, transform, and page modules the gates left at zero.

These modules are short but each carries a real branch: an unavailable backend,
a rejected transform, a missing PyMOL proxy method. None of them had test cover,
so a regression in any would have shipped silently.
"""

from __future__ import annotations

import importlib
from pathlib import Path

import numpy as np
import pytest

from structlens.integrations.freesasa.adapter import calculate_sasa
from structlens.integrations.muscle.executable import bundled_executable
from structlens.integrations.pymol.state_snapshot import PyMOLStateSnapshot
from structlens.integrations.pymol.transforms import apply_transform
from structlens.plugin.gui.widgets import PageDescriptor
from structlens.resources.backends import backend_versions

PAGE_MODULES = (
    "alignment_page",
    "mutations_page",
    "project_page",
    "residues_page",
    "results_page",
    "visualization_page",
)


class _RecordingCommand:
    """Stands in for the PyMOL command proxy."""

    def __init__(self) -> None:
        self.transforms: list[tuple[str, list[float]]] = []
        self.deleted: list[str] = []

    def transform_object(self, name: str, matrix: list[float]) -> None:
        self.transforms.append((name, matrix))

    def delete(self, name: str) -> None:
        self.deleted.append(name)


def test_apply_transform_sends_a_flat_twelve_value_matrix() -> None:
    command = _RecordingCommand()
    apply_transform(command, "target", np.eye(3), (1.0, 2.0, 3.0))
    name, matrix = command.transforms[0]
    assert name == "target"
    assert matrix == [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 3.0]


@pytest.mark.parametrize(
    ("rotation", "translation"),
    [(np.eye(2), (0.0, 0.0, 0.0)), (np.eye(3), (0.0, 0.0))],
)
def test_apply_transform_rejects_malformed_input(rotation: np.ndarray, translation: tuple[float, ...]) -> None:
    with pytest.raises(ValueError, match="rotation must be 3x3"):
        apply_transform(_RecordingCommand(), "target", rotation, translation)


def test_apply_transform_reports_a_proxy_without_transform_object() -> None:
    with pytest.raises(RuntimeError, match="does not expose transform_object"):
        apply_transform(object(), "target", np.eye(3), (0.0, 0.0, 0.0))


def test_state_snapshot_deletes_every_owned_name() -> None:
    command = _RecordingCommand()
    PyMOLStateSnapshot(("sel_a",), ("obj_a", "obj_b")).restore(command)
    assert command.deleted == ["sel_a", "obj_a", "obj_b"]


def test_state_snapshot_tolerates_a_proxy_without_delete() -> None:
    PyMOLStateSnapshot(("sel_a",), ()).restore(object())


def test_backend_versions_always_reports_every_backend() -> None:
    versions = backend_versions()
    assert versions["bundle_schema"] == "3"
    assert set(versions) >= {"StructLens", "MUSCLE", "US-align", "FreeSASA"}
    for name in ("MUSCLE", "US-align", "FreeSASA"):
        # A backend that is absent must say so rather than report a blank.
        assert versions[name]


def test_bundled_muscle_is_unavailable_on_unsupported_platforms() -> None:
    assert bundled_executable(system="Darwin", machine="arm64") is None
    assert bundled_executable(system="Linux", machine="riscv64") is None


@pytest.mark.parametrize(("system", "machine"), [("Windows", "AMD64"), ("Linux", "x86_64")])
def test_bundled_muscle_resolves_to_a_real_file_or_none(system: str, machine: str) -> None:
    """Source checkouts carry no binaries, so None is the honest answer here."""

    resolved = bundled_executable(system=system, machine=machine)
    assert resolved is None or resolved.is_file()


def test_calculate_sasa_returns_none_instead_of_zero_for_a_missing_file(tmp_path: Path) -> None:
    assert calculate_sasa(tmp_path / "absent.pdb") is None


def test_page_descriptors_are_named_and_documented() -> None:
    for module_name in PAGE_MODULES:
        module = importlib.import_module(f"structlens.plugin.gui.{module_name}")
        page = module.PAGE
        assert isinstance(page, PageDescriptor)
        assert page.name and page.purpose


def test_plugin_entrypoint_requires_a_pymol_host() -> None:
    from structlens.plugin.entrypoint import __init_plugin__

    with pytest.raises(RuntimeError, match="must be loaded from a PyMOL environment"):
        __init_plugin__()
