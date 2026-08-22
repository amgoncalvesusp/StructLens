# StructLens architecture

StructLens is divided into a PyMOL-independent scientific core, application services, external adapters, a CLI, and a thin plugin UI. `ResidueCorrespondence` is the authoritative map. Parsing produces normalized `ProteinStructure`/`ProteinChain` records with source numbering and atom coordinates. Sequence and structural engines propose mappings; analysis services validate them, calculate geometry, detect mutation descriptors, and serialize `AnalysisResult`. PyMOL receives already-computed results and only applies namespaced selections, transforms, and representations.

The plugin is optional: importing `structlens.core` and `structlens.application` never imports PyMOL or Qt. GUI pages are organized as Project, Alignment, Mutations, Residues, Visualization, and Results. UI state is separate from scientific state, and visualization reset affects only StructLens-owned names.

## Impeccable gate

`impeccable` was discovered at `C:\Users\adria\.codex\skills\impeccable\SKILL.md` and invoked for the GUI architecture pass. The visual direction is recorded in `DESIGN.md`. A second critique/polish pass is recorded after the optional GUI surface is complete.
