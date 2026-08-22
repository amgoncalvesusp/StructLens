# Scientific methods

## Mapping

Sequence mode uses Biopython's global PairwiseAligner with an explicit BLOSUM62 matrix and serialized gap-open/gap-extend settings. AUTO reports the identity and reference coverage used to choose sequence-guided mapping; below the configured thresholds it requires the US-align structural adapter. Manual mappings are explicit and may be marked locked.

## Geometry

Kabsch fits target row-vector coordinates onto reference coordinates and explicitly corrects reflections. Strict Cα RMSD uses all eligible mapped Cα pairs. Refinement iterates a configured Å cutoff and retains excluded correspondences with `is_outlier=True`. A one-pair Cα value is displayed as Cα displacement, never RMSD.

Per-residue backbone RMSD requires N, Cα, C, and O in both residues. Side-chain and all-heavy-atom RMSDs require complete matching heavy-atom name sets; valid ASP/GLU/ARG/PHE/TYR symmetry permutations are evaluated and the minimum chemically valid RMSD is used. Missing atoms produce an unavailable metric.

## Mutation descriptors

Mutation classification happens after mapping and distinguishes conserved, substitution, insertion, deletion, and non-standard residues. BLOSUM62 and Grantham are embedded, test-verified descriptors. Physicochemical categories are descriptive and carry no functional claim.

## Reproducibility

Project JSON stores schema version, source paths, mapping settings, key residues, visualization state, correspondence tables, mutation descriptors, exclusions, and provenance. Source numbering preserves author sequence IDs, insertion codes, and mmCIF label sequence IDs.
