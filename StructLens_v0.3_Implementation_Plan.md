# STRUCTLENS v0.3.0 — Implementation Plan

This file records the v0.3.0 plan supplied in the implementation request. The
complete normative text is the user message in the task thread; this local
ledger copy keeps the execution order and release gates durable across context
compaction.

## Binding scientific capabilities

1. Multiple sequence alignment, alignment conservation, insertion columns,
   sequence logo, and provenance.
2. Structural neighborhoods and reference-normalized interaction differences.
3. Key-residue, ligand-radius, and residue-radius site definitions with
   global-frame and site-fitted RMSD, SASA, and atomic-envelope volume.
4. Distance-difference matrices and transformed structural displacement vectors.
5. Residue Evidence Card with sequence, structure, interactions, site, and
   evidence-quality sections; never an impact/damage/function score.

## Required contracts

- Preserve `ResidueId` source numbering and authoritative correspondence.
- Represent missing/unavailable values as `None`, never as zero.
- MSA column → sequence index → source `ResidueId`; preserve reference-gap and
  reference-relative insertion columns.
- Entropy excludes gaps and ambiguous residues from canonical amino-acid
  counts; gap and ambiguous fractions remain separate.
- Interaction comparison maps target residues to reference positions first.
- Heavy-atom hydrogen bonds are putative (`heavy_atom_geometry`) without
  explicit hydrogens.
- Site atomic-envelope volume is a convex-hull envelope, not cavity volume.
- Distance differences are internal and require no superposition.
- Vectors apply exactly one target→reference transform and match stored Cα
  displacement magnitudes.
- PyMOL consumes bundle results and never recomputes science.

## Required work order

1. Freeze architecture, v0.3 version, shared domain types, and service
   protocols.
2. Bundle MUSCLE and verify redistribution metadata.
3. Parallel-safe wave: MSA/conservation; interactions; sites; difference maps.
4. Review and integrate the four scientific modules.
5. Implement Evidence Card, schema v3, XLSX/charts, GUI, bundle, and plugin.
6. Build Windows Setup, Linux Debian, Linux AppImage, and plugin ZIP.
7. Run cross-platform parity, end-to-end fixture, documentation, and release
   automation checks.

## GUI navigation

`Project`, `Sequences`, `Structures`, `Residues`, `Sites`, `Charts`, `PyMOL`,
`Results`, `Export`.

Before GUI changes, use the available `impeccable` skill and repeat its final
critique after implementation.

## Release gates

- pytest, Ruff, and mypy pass.
- Golden tests cover MSA, conservation, interactions, sites, distance maps,
  vectors, and Evidence Cards.
- Windows/Linux parity and clean installs pass.
- US-align, MUSCLE, and FreeSASA work offline from installers.
- XLSX, JPEG, and genuine 300/600 dpi TIFF exports pass.
- `.structlens-pymol` v0.3 exports and imports in the plugin.
- Plugin interactions, sites, vectors, Evidence Inspector, reset, and image
  exports pass.
- Documentation, licenses, citations, checksums, and release notes are
  complete.

## Branches and artifacts

- Main branch: `feat/structlens-v0.3.0`, release tag `v0.3.0`.
- Plugin branch: `feat/structlens-pymol-v0.3.0`, tag
  `structlens-pymol-v0.3.0`.
- Main artifacts: `StructLens-v0.3.0-Windows-x86_64-Setup.exe`,
  `structlens_0.3.0_amd64.deb`,
  `StructLens-v0.3.0-Linux-x86_64.AppImage`, with SHA-256 files.
- Plugin artifacts: `StructLens-PyMOL-v0.3.0.zip` and
  `StructLens-PyMOL-v0.3.0.sha256`.
