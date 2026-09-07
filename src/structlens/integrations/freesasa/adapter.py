"""FreeSASA adapter with no silent zero substitution."""

from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory


class FreeSASAAdapter:
    def __init__(self) -> None:
        try:
            import freesasa  # type: ignore[import-not-found]
        except ImportError as error:
            raise RuntimeError("FreeSASA is not installed in this runtime") from error
        self._freesasa = freesasa

    def calculate_file(self, path: str | Path) -> float:
        structure = self._freesasa.Structure(str(path))
        result = self._freesasa.calc(structure)
        return float(result.totalArea())

    def calculate_pdb(self, pdb_text: str) -> float:
        # FreeSASA accepts a filename, not PDB text. Close the file before the
        # native library opens it, including on Windows.
        with TemporaryDirectory(prefix="structlens-sasa-") as directory:
            path = Path(directory) / "structure.pdb"
            path.write_text(pdb_text, encoding="utf-8")
            return self.calculate_file(path)


def calculate_sasa(path: str | Path) -> float | None:
    try:
        return FreeSASAAdapter().calculate_file(path)
    except (ImportError, OSError, RuntimeError, ValueError):
        return None


__all__ = ["FreeSASAAdapter", "calculate_sasa"]
