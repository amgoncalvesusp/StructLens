# StructLens architecture

StructLens is divided into a PyMOL-independent scientific core, application services, external adapters, a CLI, and a thin desktop UI. `ResidueCorrespondence` is the authoritative map. Parsing produces normalized `ProteinStructure`/`ProteinChain` records with source numbering and atom coordinates. Sequence and structural engines propose mappings; analysis services validate them, calculate geometry, detect mutation descriptors, and serialize `AnalysisResult`. Multi-structure services build `ReferenceVsManyAnalysis`, `AllVsAllAnalysis`, and `MultipleStructureAnalysis` without renumbering source residues. PyMOL receives already-computed results through a validated data-only `.structlens-pymol` bundle and only applies namespaced selections, transforms, and representations.

The plugin is optional: importing `structlens.core` and `structlens.application` never imports PyMOL or Qt. The Qt surface is an operate-mode Evidence Bench with the v0.2 scientific labels Project, Sequences, Structures, Residues, Charts, and Results; legacy host identifiers remain source-compatible. It loads files, chains, or existing PyMOL objects, runs comparisons in a background executor with progress/cancellation, exposes non-editable evidence tables, persists projects, and exports the current result or a `.structlens-pymol` bundle. UI state is separate from scientific state. The injected PyMOL adapter copies selected atoms into StructLens-owned temporary objects before applying representations/colors/labels, so reset or panel close can remove only StructLens-owned state without mutating the user's original objects.

## v0.2 boundaries

- The bundled-backend resolver selects custom executable, packaged platform binary, then PATH only as a diagnostic fallback. Source distributions contain the upstream citation/license metadata; release maintainers must add and checksum a redistributable binary after reviewing upstream terms.
- `structlens.application.chart_data` is renderer-independent. The six chart datasets are the only source for Qt/matplotlib/XLSX values, preserving units and missing data.
- `.structlens-pymol` is a ZIP data contract, never an executable plugin payload. `validate_pymol_bundle` rejects traversal, duplicate entries, executable suffixes, malformed manifests, unknown structure IDs, and future schema majors.
- StructLens-PyMOL is a separate product in `amgoncalvesusp/pymol-plugins`; this repository contains only the schema, writer, docs, and optional launch handoff.

## Impeccable gate

`impeccable` was discovered at `C:\Users\adria\.codex\skills\impeccable\SKILL.md` and invoked for the GUI architecture pass. The visual direction is recorded in `DESIGN.md`. A second critique/polish pass is recorded after the optional GUI surface is complete.
