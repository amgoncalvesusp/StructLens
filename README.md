# StructLens

StructLens compares protein sequences and structures by explicit residue correspondence. The desktop application links each mutation to its evidence, lets you inspect ligand contacts and distance changes, and compares multiple variants against a fixed reference.

Author: Adriano Marques Gonçalves (UNIARA)

See the [v0.5.0 release notes](RELEASE_NOTES_v0.5.0.md) for the current workflows and the [scientific methods](docs/scientific-methods.md) for measurement definitions.

## Installers and releases

Download the x86-64 Windows Setup, Linux Debian package or Linux AppImage from [GitHub Releases](https://github.com/amgoncalvesusp/StructLens/releases). Each asset includes a separate SHA-256 checksum.

## Install

```bash
python -m pip install -e ".[gui,charts]"
```

PDB, mmCIF/CIF (including gzip-compressed PDB), and FASTA inputs are supported. NumPy, Biopython, openpyxl, and Pillow are installed as runtime dependencies for analysis and export. Install `structlens[charts]` for publication chart images (matplotlib). Structural mapping resolves a custom US-align path first, then a tested bundled binary when present, and finally PATH only as a developer/diagnostic fallback; StructLens never downloads binaries automatically.

## First comparison

```bash
structlens compare reference.pdb target.pdb --mode auto --output comparison.xlsx
```

The command reports sequence identity, coverage, strict Cα RMSD, optional refined RMSD, mapped residues, mutation count, and the selected mapping branch. Use `--csv` or `--json` for additional exports.

## Desktop GUI

Start `structlens-gui` from a source installation, or open StructLens from the installed shortcut. Official desktop builds use PySide6.

1. Load a reference and target in **Project**, choose the model and chain, then click **Compare structures**. For multiple variants, use the batch dialog in **Project**.
2. Click a mutation count in **Results** to filter **Residues** and open the evidence card. Search for an author position such as `A:103B`, or select the corresponding column in **Sequences**.
3. In **Sites & pockets**, pick the ligand instance, set a radius and preview the residue list. Use the contact map or table to inspect a residue. DockLens results require the matching receptor and original pose file.
4. In **Charts**, select the distance-difference heatmap to inspect a measured pair or export its matrix. In **PyMOL**, save or open a scene for the last verified comparison.

For an interface within one complex, load its coordinates as the reference and open **Project → Subunit interfaces…**. Select two protein or peptide chains in the same model. Inspect the closest atom pair for each contacting residue pair, rank residues by the number of distinct partners on the opposite chain, and export the evidence as CSV. The default contact cutoff is 5.0 Å.

## Scientific terminology

`Cα displacement` is the distance between one aligned Cα pair. It is not a residue RMSD. Backbone RMSD uses matched N/Cα/C/O atoms; side-chain RMSD uses matched side-chain heavy atoms and valid symmetry permutations. Strict RMSD preserves every eligible pair. Refined RMSD is a separate result with an explicit cutoff and excluded alignment indices; no outlier is silently deleted.

BLOSUM62, Grantham distance, structural displacement, and mutation classes are descriptors only. StructLens does not infer pathogenicity, stability, catalysis, or function.

## PyMOL integration

For verified canonical comparisons, StructLens writes a self-contained `.pml` file with captured coordinates, the computed fit and named views. Open it in PyMOL or configure a PyMOL executable in the desktop application. The legacy `.structlens-pymol` workflow uses the [companion plugin](https://github.com/amgoncalvesusp/pymol-plugins/).

## Exports and limitations

Export the completed analysis from **Export**. Chart datasets support XLSX and, with the chart renderer installed, publication JPEG/TIFF. Displayed residue distances are in Å. Missing measurements stay unavailable. A DockLens docking pose and a deposited ligand have different coordinates and retain separate provenance.
