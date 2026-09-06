# StructLens v0.4.0

## Highlights

- Canonical pairwise reports now retain typed pocket detections, volume measurements, matching states, diagnostics, and provenance.
- Quality evidence is visible in the GUI with native error, warning, and informational diagnostics.
- Sites & Pockets evidence shows candidate identity, alpha-sphere counts, lining counts, matched/ambiguous/unmatched states, coarse/fine free volume, and sensitivity in Å³.
- Canonical request provenance is checked before report binding or export.
- PyMOL object and selection ownership is collision-safe; temporary coordinate files are removed immediately after use.
- Raw reference and target interaction records remain visible even when no difference is detected.

## Scientific limits

StructLens reports descriptive structural evidence. It does not predict function,
pathogenicity, stability, binding affinity, free energy, or druggability. Atomic
envelope volume and pocket free volume are distinct measurements. Pocket results
depend on the selected structure, coordinate quality, atom scope, radii table,
grid settings, and explicit provenance retained in the report.

## Packaging

The GitHub tag workflow builds the Windows setup, Linux Debian package, and Linux
AppImage. MUSCLE and other external backends remain separately licensed components;
see `THIRD_PARTY_NOTICES.txt` and `licenses/`.
