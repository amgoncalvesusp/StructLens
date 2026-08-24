# Multi-structure analysis

StructLens v0.2 keeps pairwise analysis and adds three explicit topologies:

- **Reference vs many** runs an independent alignment for every target. Each target keeps its own correspondence rows, mutations, structural metrics, and transform while all display positions remain anchored to the reference.
- **All vs all** evaluates each unordered pair once (`N × (N - 1) / 2`). `PairwiseMatrix` stores one canonical value per pair and mirrors it only when read or plotted. Missing/failed pairs remain `None`.
- **Multiple structure alignment** exposes common aligned positions, coverage, and Cα positional variability. The public terminology is positional variability/structural variability, never trajectory RMSF.

The pairwise result remains the authoritative scientific record. Multi-analysis wrappers do not recalculate values and preserve target residue locators, including chains and insertion codes. `AnalysisSelection` is the shared future-facing selection contract for residues, charts, and PyMOL export.

Sequence identity, similarity, coverage, BLOSUM62, Grantham distance, and mutation classes are descriptors. They do not imply function, stability, pathogenicity, or energetic effects.
