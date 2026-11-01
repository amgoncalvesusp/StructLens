"""Exercise the CLI, science report and Qt chart installed by one native package."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("cli", type=Path)
    parser.add_argument("gui", type=Path)
    args = parser.parse_args()
    structure = Path(__file__).resolve().parents[1] / "tests" / "fixtures" / "parsing" / "enriched.pdb"
    env = {**os.environ, "QT_QPA_PLATFORM": "offscreen", "QT_API": "pyside6", "PYTHONIOENCODING": "utf-8"}
    env.pop("PYTHONPATH", None)
    with TemporaryDirectory() as scratch:
        report = Path(scratch) / "canonical.json"
        subprocess.run(
            [str(args.cli.resolve()), "compare", str(structure), str(structure), "--mode", "sequence",
             "--report", str(report)], cwd=scratch, env=env, timeout=180, check=True,
        )
        analysis = json.loads(report.read_text(encoding="utf-8"))["analysis"]
        if not analysis or analysis["sequence_identity"] != 1 or analysis["mapped_residue_count"] < 1:
            raise RuntimeError("frozen CLI failed the known self-comparison")
        subprocess.run(
            [str(args.gui.resolve()), "--smoke-report", str(report), str(structure)],
            cwd=scratch, env=env, timeout=180, check=True,
        )


if __name__ == "__main__":
    main()
