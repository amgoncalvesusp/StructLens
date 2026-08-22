# Validation

The scientific unit suite covers identity, substitutions, insertions/deletions, numbering offsets, PDB insertion codes, mmCIF author/label numbering, alternate locations, gzip PDB input, rigid-body invariance, reflection-safe Kabsch fitting, missing atoms, symmetric side-chain names, refinement outliers, mutation descriptors, and US-align parser fixtures.

Numerical comparisons use NumPy floating-point tolerances appropriate for Å coordinates. A self-comparison is expected to have zero mutations and RMSD at floating-point zero. Raw US-align text is parsed into explicit pairs and provenance; stdout is never persisted as scientific state.
