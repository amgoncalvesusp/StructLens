"""FreeSASA adapter with no silent zero substitution."""

from __future__ import annotations

from pathlib import Path


class FreeSASAAdapter:
    def __init__(self) -> None:
        try:
            import freesasa  # type: ignore[import-not-found]
        except ImportError as error:
            raise RuntimeError("FreeSASA is not installed in this runtime") from error
        self._freesasa = freesasa

    def calculate_file(self, path: str | Path) -> float:
        result = self._freesasa.calc(str(path))
        return float(result.totalArea())

    def calculate_pdb(self, pdb_text: str) -> float:
        structure = self._freesasa.Structure(pdb_text)
        return float(self._freesasa.calc(structure).totalArea())


def calculate_sasa(path: str | Path) -> float | None:
    try:
        return FreeSASAAdapter().calculate_file(path)
    except (ImportError, OSError, RuntimeError, ValueError):
        return None


__all__ = ["FreeSASAAdapter", "calculate_sasa"]
