# StructLens

StructLens is an English-only desktop analysis application for integrated protein sequence and structure comparison. v0.3 adds bundled MSA/conservation, structural interactions, site geometry, distance-difference maps, displacement vectors, and Residue Evidence Cards while preserving an explicit residue correspondence table.

Author: Adriano Marques Gonçalves (UNIARA)

See [docs/v0.3.md](docs/v0.3.md) for the scientific methods, limitations,
offline installation artifacts, reproducibility notes, and PyMOL plugin workflow.

## Installers and releases

Tagged releases publish a Python wheel plus native setup installers for Windows x64 and Linux x86_64. On Windows, double-click `StructLens-Setup.exe`; on Linux, open `StructLens-Setup.deb` with the system software/package installer. Both setups install the standalone CLI and desktop GUI for the main StructLens application.

## Install

```bash
python -m pip install -e .
```

PDB, mmCIF/CIF (including gzip-compressed PDB), and FASTA inputs are supported. NumPy, Biopython, openpyxl, and Pillow are installed as runtime dependencies for analysis and export. Install `structlens[charts]` for publication chart images (matplotlib). Structural mapping resolves a custom US-align path first, then a tested bundled binary when present, and finally PATH only as a developer/diagnostic fallback; StructLens never downloads binaries automatically.

## First comparison

```bash
structlens compare reference.pdb target.pdb --mode auto --output comparison.xlsx
```

The command reports sequence identity, coverage, strict Cα RMSD, optional refined RMSD, mapped residues, mutation count, and the selected mapping branch. Use `--csv` or `--json` for additional exports.

## Desktop GUI

Run `structlens-gui` after installing the `gui` extra (PyQt5), or launch StructLens from the Windows/Linux setup shortcut. The standalone Evidence Bench loads PDB/mmCIF files, runs comparisons without freezing the window, exposes Project, Sequences, Structures, Residues, Charts, PyMOL, Results, and Export areas, and writes XLSX/CSV/JSON plus validated `.structlens-pymol` bundles. PySide6 users can opt in with `pip install "structlens[gui-pyside6]"`.

## Scientific terminology

`Cα displacement` is the distance between one aligned Cα pair. It is not a residue RMSD. Backbone RMSD uses matched N/Cα/C/O atoms; side-chain RMSD uses matched side-chain heavy atoms and valid symmetry permutations. Strict RMSD preserves every eligible pair. Refined RMSD is a separate result with an explicit cutoff and excluded alignment indices; no outlier is silently deleted.

BLOSUM62, Grantham distance, structural displacement, and mutation classes are descriptors only. StructLens does not infer pathogenicity, stability, catalysis, or function.

## PyMOL integration

PyMOL is a separate visualization product in `https://github.com/amgoncalvesusp/pymol-plugins/` under `structlens-pymol/`. This repository owns the scientific analysis, bundle schema, bundle writer, and optional external launch handoff. The desktop app can always export a validated `.structlens-pymol` file; when a PyMOL executable is configured, `Open in PyMOL` launches the external application with that bundle.

## Exports and limitations

Every scientific table exposed by the application has XLSX/CSV/JSON export paths. Chart datasets can be exported to XLSX or publication JPEG/TIFF at real 300/600 dpi when the `charts` extra is installed. Online structure downloads, energetic/function prediction, trajectories, cloud services, and live IPC synchronization remain outside this revision.
