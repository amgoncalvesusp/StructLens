# StructLens

StructLens is an English-only Python engine and optional PyMOL plugin for reproducible protein sequence/structure comparison. It keeps an explicit residue correspondence table, detects descriptive mutation classes, calculates strict and optional refined structural metrics, and exports scientific tables.

Author: Adriano Marques Gonçalves (UNIARA)

## Installers and releases

Tagged releases publish a Python wheel plus native setup installers for Windows x64 and Linux x86_64. On Windows, double-click `StructLens-Setup.exe`; on Linux, open `StructLens-Setup.deb` with the system software/package installer. Both setups install the standalone CLI and include the wheel and PyMOL plugin installer files.

For PyMOL, run the included plugin installer script against the Python executable used by your PyMOL installation, then load `structlens.plugin.entrypoint` through Plugin > Install Plugin.

## Install

```bash
python -m pip install -e .
```

PDB, mmCIF/CIF (including gzip-compressed PDB), and FASTA inputs are supported. NumPy, Biopython, openpyxl, and Pillow are installed as runtime dependencies for analysis and export. Structure-guided mapping uses an explicitly configured US-align executable; StructLens never downloads binaries automatically.

## First comparison

```bash
structlens compare reference.pdb target.pdb --mode auto --output comparison.xlsx
```

The command reports sequence identity, coverage, strict Cα RMSD, optional refined RMSD, mapped residues, mutation count, and the selected mapping branch. Use `--csv` or `--json` for additional exports.

## Scientific terminology

`Cα displacement` is the distance between one aligned Cα pair. It is not a residue RMSD. Backbone RMSD uses matched N/Cα/C/O atoms; side-chain RMSD uses matched side-chain heavy atoms and valid symmetry permutations. Strict RMSD preserves every eligible pair. Refined RMSD is a separate result with an explicit cutoff and excluded alignment indices; no outlier is silently deleted.

BLOSUM62, Grantham distance, structural displacement, and mutation classes are descriptors only. StructLens does not infer pathogenicity, stability, catalysis, or function.

## PyMOL plugin

Load `structlens.plugin.entrypoint` from the PyMOL plugin manager. The Evidence Bench GUI provides Project, Alignment, Mutations, Residues, Visualization, and Results stages: choose files or existing PyMOL objects, select chains, run a background comparison, inspect double-clickable tables, apply a preset, and export the result. The plugin is namespaced and reversible: it only deletes selections it created. PyMOL is a visualization backend, never the scientific state. If Qt/PyMOL is unavailable, the core and CLI remain usable.

The panel sections are Project, Alignment, Mutations, Residues, Visualization, and Results. Each workflow option has contextual English help. Visualization presets are Minimal, Publication, Mutation focus, Structural deviation, Active site, and Presentation.

## Exports and limitations

Every scientific table exposed by the application has XLSX/CSV/JSON export paths. The core export service currently provides publication-ready tabular output; raster JPEG/TIFF rendering requires a host molecular renderer and is intentionally isolated from the scientific engine. Online structure downloads, energetic/function prediction, trajectories, cloud services, and PDF/HTML reports are deferred v0.1 features.
