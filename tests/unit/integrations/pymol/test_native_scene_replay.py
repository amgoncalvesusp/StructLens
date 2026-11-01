"""Exercise scene replacement through the installed native PyMOL reader."""

import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from tests.unit.plugin._task13_round5_fixtures import multi_selected_report

from structlens.integrations.pymol.scene import report_scene_data, write_scene_data


def test_two_native_focuses_variants_and_failure_preserve_user_objects(tmp_path):
    executable = shutil.which("pymol")
    python = Path(sys.executable) if importlib.util.find_spec("pymol2") else None
    if python is None and executable and os.name == "nt":
        python = Path(executable).resolve().parents[1] / "python.exe"
    if python is None or not python.is_file():
        pytest.skip("native PyMOL Python is not installed")
    request, report = multi_selected_report(tmp_path)
    snapshots = (request.reference_snapshot, request.target_snapshot)
    first, second = report.analysis.correspondences[:2]
    initial = report_scene_data(report, snapshots, first.alignment_index)
    updated = report_scene_data(report, snapshots, second.alignment_index, show_target=False)
    write_scene_data(tmp_path / "first.pml", initial)
    write_scene_data(tmp_path / "second.pml", updated)
    (tmp_path / "data.json").write_text(json.dumps(updated))
    runtime = Path(__file__).parents[4] / "src/structlens/integrations/pymol/scene_runtime.py"
    script = r"""
import json, sys, runpy
from pathlib import Path
import pymol2
from pymol import cmd
root = Path(sys.argv[1])
pm = pymol2.SingletonPyMOL()
pm.start()
try:
    cmd.pseudoatom("User_object")
    cmd.load(str(root / "first.pml"))
    assert len(cmd.get_names("objects")) == 3
    cmd.load(str(root / "second.pml"))
    assert len(cmd.get_names("objects")) == 3, cmd.get_names("objects")
    assert "User_object" in cmd.get_names("objects", enabled_only=1)
    assert not any("_target" in n for n in cmd.get_names("objects", enabled_only=1))
    assert len(cmd.get_scene_list()) == 2
    data = json.loads((root / "data.json").read_text())
    data["report"] = "c" * 64
    render = runpy.run_path(sys.argv[2])["render_scene"]
    render(cmd, data)
    assert len(cmd.get_names("objects")) == 3
    assert all("cccccccc" in n for n in cmd.get_names("objects") if n != "User_object")
    before = (sorted(cmd.get_names("objects", enabled_only=1)), list(cmd.get_scene_list()), tuple(cmd.get_view()))
    data["sources"][1]["chain"] = "does not exist; delete all"
    try:
        render(cmd, data)
    except ValueError:
        pass
    else:
        raise AssertionError("malformed selection succeeded")
    assert before[0] == sorted(cmd.get_names("objects", enabled_only=1))
    assert before[1] == list(cmd.get_scene_list())
    assert all(abs(a-b) < 0.00001 for a,b in zip(before[2],cmd.get_view()))
finally:
    pm.stop()
"""
    env = dict(os.environ)
    env.pop("PYTHONPATH", None)
    if os.name == "nt":
        env["PATH"] = os.pathsep.join((str(python.parent / "Library/bin"), str(python.parent), env["PATH"]))
    result = subprocess.run([str(python), "-c", script, str(tmp_path), str(runtime)], capture_output=True, text=True, env=env, timeout=90)
    assert result.returncode == 0, (result.stdout + result.stderr)[-4000:]
