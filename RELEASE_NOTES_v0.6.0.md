# StructLens v0.6.0

Open **Project → Subunit interfaces…**, select two chains from one model of the reference structure, and calculate their contacts. Select a residue in the ranking to inspect its partners, nearest atoms and distances. Export the evidence to CSV.

The ranking counts distinct contacting residues. Salt bridges, hydrogen bonds and hydrophobic contacts are geometric candidates. These measurements do not establish binding energy or the energetic contribution of a residue. To study a biological assembly, load coordinates containing the relevant subunits.

Verified with HIV-1 protease and its tripeptide in PDB entry 1A30. At 5 Å, chains A and B have 138 contacting residue pairs. Chains A and C have 22, including 12 chain A residues contacting peptide Asp C:507. PDB and mmCIF inputs gave identical counts.
