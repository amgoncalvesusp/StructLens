# StructLens v0.3.1

## Fixes

- Comparison preserves the selected reference and target chains.
- Results and Export navigation now opens the matching page.
- Reloading same-named PyMOL files replaces only StructLens-owned objects.
- Multiple-structure alignment keeps reference residue identity across target insertions and deletions.
- Refined RMSD, inlier counts, and excluded indices now use one final fit.
- Sequence and structure analyses export the transform used by their geometry metrics.

## Validation

The release passes the full test suite, Ruff, mypy, and offscreen Qt checks.
