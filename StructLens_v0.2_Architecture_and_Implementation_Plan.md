# StructLens v0.2 Architecture and Implementation Plan

> **For Luna Extra High / agentic implementation**
>
> **REQUIRED SUPERPOWERS WORKFLOW:** use `superpowers:using-git-worktrees`, `superpowers:test-driven-development`, `superpowers:subagent-driven-development`, `superpowers:requesting-code-review`, and `superpowers:verification-before-completion`.
>
> This plan is an approved architectural revision of the existing StructLens project. It supersedes conflicting decisions in earlier plans where PyMOL was treated as an internal/embedded visualization backend.
>
> **GUI requirement:** before changing the StructLens desktop GUI or creating the StructLens-PyMOL GUI, discover whether the `impeccable` skill is available. If available, use it before GUI implementation and again for the final GUI critique/polish pass. If unavailable, record that fact and follow the GUI requirements in this document.

**Goal:** evolve StructLens from a two-structure comparison application into an integrated sequence-and-structure analysis platform with bundled US-align, multi-structure analyses, scientific charts, a reproducible PyMOL interchange format, and a separately released StructLens-PyMOL visualization plugin.

**Main implementation model:** Luna Extra High.

**User-facing language:** English only.

---

# 1. Product Identity

The product must no longer present itself primarily as:

```text
Protein structure comparison
```

Use a broader scientific identity such as:

```text
Integrated protein sequence and structure comparison
```

or, where space is limited:

```text
Residue-level sequence and structural analysis
```

The product differentiator is the integration of:

```text
sequence comparison
        +
residue correspondence
        +
mutation analysis
        +
multiple structural comparison
        +
residue-level structural metrics
        +
comparative charts
        +
PyMOL exploration
```

StructLens must not become merely a graphical wrapper around US-align.

US-align is an important structural-alignment backend. StructLens owns:

- project state;
- input normalization;
- sequence analysis;
- residue correspondence;
- mutation analysis;
- structural metrics;
- multi-structure summary metrics;
- charts;
- exports;
- provenance;
- PyMOL interchange data.

---

# 2. Approved High-Level Architecture

Implement two separate products with a versioned interchange contract.

```text
                         STRUCTLENS
                  Desktop analysis software
                            │
        ┌───────────────────┼───────────────────┐
        │                   │                   │
     Sequence            Structure           Residue
     Analysis            Analysis            Analysis
        │                   │                   │
        ├──────────────┬────┴────┬──────────────┤
                       │         │
                     Charts   Results
                       │         │
        ┌──────────────┴────┬────┴──────────────┐
        │                   │                   │
      XLSX              JPEG / TIFF       .structlens-pymol
                                                │
                                                ▼
                                      STRUCTLENS-PYMOL
                                        separate plugin
                                                │
                                                ▼
                                              PyMOL
```

The PyMOL plugin is **not** part of the StructLens Python environment.

Do not require the PyMOL Python interpreter to import the complete StructLens application.

Do not require StructLens to import PyMOL.

The interchange bundle is the stable contract between them.

---

# 3. Repository Boundaries

## 3.1 Main StructLens repository

The existing StructLens repository remains responsible for:

- desktop GUI;
- sequence analysis;
- structure analysis;
- US-align integration;
- multi-structure analysis;
- mutations;
- residue metrics;
- charts;
- XLSX/CSV/TSV/JSON export;
- JPEG/TIFF publication export;
- creation of `.structlens-pymol` bundles;
- optional launching of the external PyMOL application.

The main repository must **not** contain the distributable PyMOL plugin as its canonical source.

It may contain:

```text
docs/pymol-integration.md
schemas/structlens-pymol/
tests/interchange/
```

but the plugin source belongs in the plugin repository below.

## 3.2 PyMOL plugin repository

Canonical repository:

```text
https://github.com/amgoncalvesusp/pymol-plugins
```

Create a new top-level directory:

```text
structlens-pymol/
```

Do not place the new plugin inside the existing `noncovalent-interactions` plugin.

The repository becomes conceptually:

```text
pymol-plugins/
├── README.md
├── LICENSE
├── noncovalent-interactions/
│   └── ...
└── structlens-pymol/
    ├── README.md
    ├── CHANGELOG.md
    ├── THIRD_PARTY_NOTICES.md
    ├── structlens_pymol_plugin/
    │   ├── __init__.py
    │   ├── plugin.py
    │   ├── bundle_reader.py
    │   ├── schema.py
    │   ├── controller.py
    │   ├── selections.py
    │   ├── visualization.py
    │   ├── colors.py
    │   ├── presets.py
    │   ├── commands.py
    │   ├── errors.py
    │   ├── gui/
    │   │   ├── main_dialog.py
    │   │   ├── analysis_panel.py
    │   │   ├── highlight_panel.py
    │   │   └── help_panel.py
    │   └── assets/
    │       └── structlens_icon.png
    ├── tests/
    │   ├── unit/
    │   ├── integration/
    │   └── fixtures/
    ├── scripts/
    │   └── build_release_zip.py
    └── dist/
        └── .gitkeep
```

The root `README.md` of `pymol-plugins` must gain a new plugin row for StructLens-PyMOL.

---

# 4. Release Strategy for StructLens-PyMOL

GitHub Releases are repository-wide, so StructLens-PyMOL must use its own tag namespace.

Required semantic version tag format:

```text
structlens-pymol-vMAJOR.MINOR.PATCH
```

First public release example:

```text
structlens-pymol-v0.1.0
```

Release title:

```text
StructLens-PyMOL v0.1.0
```

Required release asset:

```text
StructLens-PyMOL-v0.1.0.zip
```

Recommended additional asset:

```text
StructLens-PyMOL-v0.1.0.sha256
```

Do not use a generic repository tag such as:

```text
v0.1.0
```

because the repository contains multiple independently versioned plugins.

Each StructLens-PyMOL release must include:

- installable ZIP;
- release notes;
- supported PyMOL versions tested;
- supported `.structlens-pymol` schema versions;
- installation instructions;
- upgrade instructions;
- known limitations.

The plugin `CHANGELOG.md` tracks plugin versions independently from the main StructLens application's version.

---

# 5. US-align Must Be Bundled With StructLens

The current behavior that requires users to put US-align on `PATH` is not acceptable as the default workflow.

## 5.1 User experience

Replace a visible primary-flow field such as:

```text
US-align executable: USalign on PATH
```

with a status indicator:

```text
US-align
Bundled backend · Ready
```

An advanced settings section may expose:

```text
Use custom US-align executable...
```

but ordinary users must not need this.

## 5.2 Resource layout

Main StructLens package/bundle:

```text
resources/
└── usalign/
    ├── VERSION
    ├── LICENSE
    ├── CITATION.md
    ├── windows-x64/
    │   └── USalign.exe
    ├── linux-x64/
    │   └── USalign
    ├── linux-arm64/
    │   └── USalign
    ├── macos-x64/
    │   └── USalign
    └── macos-arm64/
        └── USalign
```

Only package platforms that are actually supported and tested by the release.

Do not ship an untested binary merely because an architecture directory exists.

## 5.3 Backend selection policy

Resolution order:

1. explicit user-selected custom executable;
2. bundled platform executable;
3. PATH fallback only for developer/diagnostic builds.

The normal release must succeed without PATH configuration.

## 5.4 Provenance

Each structural analysis records:

```text
backend = US-align
US-align version
binary source = bundled | custom | PATH
platform
command options
```

## 5.5 Compressed input

StructLens owns temporary input normalization.

If the structural backend cannot consume a compressed input reliably on the current platform:

```text
*.pdb.gz
*.cif.gz
*.mmcif.gz
```

StructLens extracts it into an application-owned temporary directory before invoking US-align.

The temporary directory must be cleaned up safely.

## 5.6 Distribution compliance

Before a release containing US-align binaries:

- verify the current upstream redistribution/license terms;
- preserve the required license/copyright notices;
- preserve required citation information;
- add the dependency to `THIRD_PARTY_NOTICES.md`;
- record the exact vendored version and checksum.

Never silently replace a bundled US-align version during runtime.

---

# 6. Comparison Modes: More Than Two Structures

Multi-structure comparison is now a first-class feature.

Create:

```python
class ComparisonMode(str, Enum):
    PAIRWISE = "pairwise"
    REFERENCE_VS_MANY = "reference_vs_many"
    ALL_VS_ALL = "all_vs_all"
    MULTIPLE_STRUCTURE_ALIGNMENT = "multiple_structure_alignment"
```

The GUI must explain each mode before the user selects it.

---

# 7. Pairwise Mode

Retain the existing 1 × 1 analysis.

Use it for:

- detailed residue inspection;
- WT vs mutant;
- deep local analysis;
- publication views of one pair.

All existing pairwise metrics remain available.

---

# 8. Reference vs Many Mode

This becomes the recommended multi-protein workflow.

Example:

```text
Reference
CTX-M-15

Targets
SHV-1
TEM-1
KPC-2
Mutant-E166A
Mutant-G238S
```

Each target has its own:

- sequence alignment to the reference;
- residue correspondence table;
- rigid-body transformation;
- strict RMSD;
- refined RMSD;
- TM-score where available;
- mutation events;
- key-residue mapping;
- per-residue metrics.

All results remain anchored to the reference numbering.

This is important for biologically meaningful statements such as:

```text
Reference E166 corresponds to:
SHV-1 E166
TEM-1 E166
Protein-X E164
Mutant E166A
```

Do not collapse target numbering into reference numbering internally.

---

# 9. All vs All Mode

For N selected structures calculate all unique pairs:

```text
N × (N - 1) / 2
```

Metrics:

- sequence identity;
- sequence similarity;
- alignment coverage;
- TM-score where available;
- strict Cα RMSD;
- refined Cα RMSD if enabled;
- key-site RMSD if a key-residue group is selected.

Store results in a symmetric matrix representation without duplicating calculations.

Example:

```text
         A      B      C      D
A        —
B       ...
C       ...    ...
D       ...    ...    ...
```

All-vs-all mode is optimized for:

- heatmaps;
- clustering;
- identifying structural groups;
- identifying unusual proteins.

Do not compute expensive side-chain metrics for every pair by default.

They may be requested for selected pairs later.

---

# 10. Multiple Structure Alignment Mode

Use US-align multiple-structure alignment capabilities as the structural backend.

This mode is conceptually distinct from independently fitting each target to the reference.

Output must include:

- common aligned structural positions;
- participating structures;
- per-position coverage;
- per-position positional variability;
- per-structure deviation from the consensus/reference frame;
- mapping back to original residue locators.

Do not renumber source residues.

## 10.1 Terminology

Do not call multi-structure positional variability `RMSF`.

Use:

```text
Cα positional variability (Å)
```

or:

```text
Structural variability (Å)
```

because RMSF is strongly associated with trajectory fluctuations in molecular dynamics.

---

# 11. Sequence Analysis Must Be a First-Class Area

The application must visibly distinguish:

```text
SEQUENCE ANALYSIS
```

from:

```text
STRUCTURE ANALYSIS
```

and provide an integrated layer:

```text
SEQUENCE + STRUCTURE
```

The initial multi-protein implementation may remain reference-centered rather than introducing an additional multiple-sequence-alignment executable.

Required sequence functions:

- reference vs target sequence alignment;
- reference vs many mutation table;
- sequence identity;
- sequence similarity;
- coverage;
- substitution/insertion/deletion detection;
- physicochemical mutation class;
- BLOSUM62;
- Grantham distance;
- reference-position conservation summary;
- all-vs-all sequence identity matrix.

A true general-purpose multiple sequence alignment engine is deferred unless separately approved.

---

# 12. Integrated Sequence–Structure Questions

StructLens should make the following analyses easy to answer:

```text
Is a residue conserved in sequence but structurally displaced?
```

```text
Does a conservative substitution coincide with a large local structural change?
```

```text
Which positions are mutation hotspots but structurally conserved?
```

```text
Which positions are sequence-conserved but structurally variable across homologs?
```

```text
Which proteins have low sequence identity but high structural similarity?
```

These analyses are core differentiators.

---

# 13. New Main Navigation

Replace the current workflow emphasis:

```text
Project
Alignment
Mutations
Residues
Visualization
Results
```

with:

```text
Project
Sequences
Structures
Residues
Charts
PyMOL
Results
Export
```

Pairwise/multi-structure settings belong inside `Structures`.

Alignment policy is not removed; it becomes part of the relevant sequence/structure workflow instead of the primary identity of the application.

---

# 14. Project Page

Purpose text:

```text
Import and organize the protein structures and sequences you want to compare.
```

Required controls:

- add structure files;
- add FASTA sequence files;
- identify chains;
- select reference;
- select one or many targets;
- name structures;
- remove structures;
- inspect missing-chain/missing-sequence warnings;
- choose comparison mode.

Comparison mode cards must explain:

- what the mode does;
- when to use it;
- computational implications;
- main outputs.

---

# 15. Sequences Page

Purpose text:

```text
Compare amino-acid sequences, map equivalent positions, and inspect mutations and conservation.
```

Required views:

- pairwise aligned sequence;
- reference-vs-many mutation/conservation table;
- mutation list;
- insertion/deletion list;
- reference-position conservation;
- filters for key residues;
- sequence identity/similarity/coverage summaries.

Selection of a residue in the sequence view must update the shared application selection state.

The same selection must be usable by:

- Residues;
- Charts;
- PyMOL export presets.

---

# 16. Structures Page

Purpose text:

```text
Compare protein folds, superpositions, structural similarity, and residue-level deviations.
```

Required sections:

- comparison mode;
- alignment strategy;
- US-align bundled-backend status;
- Auto/Sequence-guided/Structure-guided/Manual mapping;
- strict vs refined RMSD settings;
- target summary;
- multi-structure analysis controls;
- selected-pair detail.

Do not expose US-align executable configuration in the primary workflow.

---

# 17. Residues Page

Purpose text:

```text
Inspect explicit residue correspondence and connect sequence changes to structural differences.
```

Required table columns:

- reference residue;
- target residue;
- amino-acid status;
- mutation notation;
- sequence score;
- Cα displacement (Å);
- backbone RMSD (Å);
- side-chain RMSD (Å);
- local RMSD (Å);
- outlier;
- key-residue group.

For Reference vs Many, allow:

```text
Target selector
```

and optionally:

```text
Compare this residue across all targets
```

which opens a multi-target residue summary.

---

# 18. Charts Page

Charts are a dedicated scientific-analysis area, not decoration.

All charts must:

- use the authoritative analysis state;
- be exportable;
- expose exact values on selection/hover;
- allow table/data export to XLSX;
- support publication-quality JPEG/TIFF output;
- use units in axes;
- not rely on color alone;
- explain what the chart means;
- explain how to interpret it;
- link selected residues/pairs back to the application selection state.

Implement the charts below.

---

# 19. Chart 1 — Structural Deviation Profile

X-axis:

```text
Reference residue position
```

Y-axis selectable:

- Cα displacement (Å);
- backbone RMSD (Å);
- side-chain RMSD (Å);
- local RMSD (Å).

For pairwise analysis:

```text
one target = one profile
```

For Reference vs Many:

- display a manageable number of selected target lines;
- when targets exceed a readable threshold, default to aggregate statistics;
- allow the user to choose individual targets.

Aggregate mode may show:

- median;
- interquartile range;
- selected percentile envelope.

Do not hide how the aggregate was calculated.

Clicking a residue point selects that residue across the linked analysis views.

---

# 20. Chart 2 — Mutation / Conservation Matrix

Rows:

```text
structures
```

Columns:

```text
reference-aligned positions
```

Cell content/color can represent:

- conserved;
- substitution;
- insertion/deletion boundary;
- non-standard residue;
- physicochemical mutation class;
- Grantham category.

Always retain text/tooltip identity such as:

```text
E166A
Reference A:GLU166
Target B:ALA164
```

Clicking a cell selects the mapped residue pair.

---

# 21. Chart 3 — Pairwise Similarity Heatmap

Supported matrices:

- TM-score;
- strict Cα RMSD;
- refined Cα RMSD;
- sequence identity;
- sequence similarity;
- key-site RMSD.

For asymmetric metrics, do not force a symmetric display.

For the listed symmetric metrics, store one calculated value per pair and mirror it visually.

Allow ordering by:

- input order;
- reference similarity;
- hierarchical clustering when implemented.

If clustering is added, document the distance transformation and linkage method.

---

# 22. Chart 4 — Sequence–Structure Relationship

Scatter plot.

Default:

```text
X = sequence identity (%)
Y = TM-score
```

Alternative Y values:

- strict Cα RMSD;
- refined Cα RMSD;
- key-site RMSD.

Each point is a target or structure pair, depending on analysis mode.

Clicking a point opens the selected pair in:

- Results;
- Residues;
- PyMOL preparation.

This chart is specifically intended to reveal cases such as:

```text
low sequence identity + high structural similarity
```

---

# 23. Chart 5 — Structural Conservation Profile

For a multiple-structure or reference-vs-many comparison:

X-axis:

```text
Reference-aligned residue position
```

Y-axis:

```text
Cα positional variability (Å)
```

Also expose:

```text
position coverage (%)
```

A position with low coverage must not be presented as equally reliable to a position mapped across all proteins.

Allow an optional companion sequence-conservation track.

---

# 24. Chart 6 — Key-Residue Comparison

Rows or X-axis:

```text
key reference residues
```

Series:

```text
targets
```

Metric selectable:

- Cα displacement;
- backbone RMSD;
- side-chain RMSD;
- local RMSD.

This is designed for:

- catalytic residues;
- binding-site residues;
- mutation hotspots;
- experimentally relevant positions.

---

# 25. Chart Export

Every chart must export:

```text
XLSX
JPEG
TIFF
```

XLSX exports the underlying plotted data, not screenshots.

Image options:

- 300 dpi;
- 600 dpi;
- publication width/height;
- white background;
- current theme background where appropriate;
- legend on/off;
- title on/off;
- annotations on/off.

For heatmaps and matrix plots, ensure exported labels remain legible at requested publication size.

---

# 26. PyMOL Page in Main StructLens

Purpose text:

```text
Prepare the current StructLens analysis for interactive 3D exploration in PyMOL.
```

This page makes the integration explicit.

Required state:

```text
PyMOL integration status
```

Possible statuses:

```text
Plugin not detected / not configured
PyMOL executable not configured
Ready to export
Ready to open in PyMOL
Bundle created
```

Required actions:

```text
Open in PyMOL
Export for PyMOL
Plugin installation instructions
```

Never display vague copy such as:

```text
3D rendering is available in the PyMOL plugin
```

without an actionable path.

---

# 27. The `.structlens-pymol` Interchange Bundle

File extension:

```text
.structlens-pymol
```

Physical format:

```text
ZIP container
```

The user should not need to unzip it manually.

Required internal layout:

```text
analysis.structlens-pymol
├── manifest.json
├── provenance.json
├── structures/
│   ├── reference.<pdb|cif>
│   ├── target_0001.<pdb|cif>
│   ├── target_0002.<pdb|cif>
│   └── ...
├── analysis/
│   ├── summary.json
│   ├── correspondence.json
│   ├── mutations.json
│   ├── key_residues.json
│   ├── outliers.json
│   └── metrics.json
├── transforms/
│   └── transforms.json
└── visualization/
    ├── default_view.json
    └── presets.json
```

Do not put arbitrary executable Python code inside the interchange bundle.

The plugin treats it as data.

---

# 28. Bundle Manifest Contract

Minimum manifest fields:

```json
{
  "format": "structlens-pymol",
  "schema_version": "1.0",
  "created_by": {
    "application": "StructLens",
    "version": "..."
  },
  "plugin_compatibility": {
    "minimum_version": "0.1.0"
  },
  "analysis_id": "...",
  "comparison_mode": "reference_vs_many",
  "reference_id": "...",
  "target_ids": ["...", "..."]
}
```

Add deterministic structure identifiers.

Do not use filenames alone as scientific identifiers.

---

# 29. Bundle Validation

Both applications implement the same schema contract independently.

Main StructLens:

```text
BundleWriter
BundleValidator
```

StructLens-PyMOL:

```text
BundleReader
BundleValidator
```

The plugin must reject:

- unsupported future major schema version;
- missing manifest;
- malformed JSON;
- missing referenced structure;
- path traversal entries;
- duplicate structure IDs;
- correspondence referring to unknown structure/residue;
- corrupt ZIP.

Never extract bundle entries with unvalidated relative paths.

---

# 30. Structure Fidelity in the Bundle

Preserve:

- chain identifiers;
- author residue numbering when available;
- insertion codes;
- residue names;
- coordinate identity;
- model identity where supported.

Do not renumber structures merely for PyMOL convenience.

If the main application generated a fitted/superposed coordinate snapshot, record that fact separately from the original coordinates.

Prefer representing transformations explicitly so the plugin can reproduce the StructLens superposition.

---

# 31. `Export for PyMOL`

This action must always work even if PyMOL is not installed.

Flow:

```text
current analysis
      ↓
validate analysis completeness
      ↓
create .structlens-pymol bundle
      ↓
validate the generated bundle
      ↓
write final file atomically
```

If validation fails, do not leave a file that appears complete.

Suggested default filename:

```text
StructLens_<reference>_<analysis-mode>.structlens-pymol
```

---

# 32. `Open in PyMOL`

This action is optional but should provide a direct user experience when PyMOL is installed.

Flow:

```text
current analysis
      ↓
create validated temporary .structlens-pymol bundle
      ↓
locate configured PyMOL executable
      ↓
launch external PyMOL process
      ↓
invoke StructLens-PyMOL import through a supported PyMOL startup mechanism
      ↓
plugin opens bundle
      ↓
structures + selections + view appear
```

Implementation requirements:

- inspect the installed PyMOL command-line/startup capabilities on each supported platform;
- use a supported startup script/command mechanism;
- quote paths safely;
- do not use shell-string interpolation;
- retain `Export for PyMOL` as the fallback if automatic launch is unavailable;
- surface a clear error if the plugin is not installed.

Acceptance is behavioral:

```text
Click Open in PyMOL
→ PyMOL starts
→ StructLens-PyMOL receives the bundle
→ reference and targets appear
→ saved alignment transform is reproduced
→ StructLens selections are available
```

The exact low-level launch flag may be platform-specific and should be isolated behind:

```python
class PyMOLLauncher:
    def launch_bundle(self, bundle_path: Path) -> LaunchResult:
        ...
```

---

# 33. No Live IPC in This Milestone

Do **not** implement:

- localhost HTTP server;
- WebSocket link;
- socket daemon;
- continuous cross-process synchronization.

Version 0.2 communication is:

```text
versioned bundle + optional launch handoff
```

A future Live Link may be planned separately.

This keeps the first integration:

- reproducible;
- debuggable;
- portable;
- firewall-independent;
- independent of Python interpreter versions.

---

# 34. StructLens-PyMOL Responsibilities

The plugin is a visualization and exploration client for StructLens results.

It does **not** recompute the main scientific analysis.

It may calculate trivial display-only geometry, but must not silently replace:

- residue mapping;
- mutation calls;
- RMSD;
- TM-score;
- structural alignment;
- outlier classification.

Those values come from the bundle.

---

# 35. StructLens-PyMOL Import Flow

When opening a bundle:

1. validate ZIP safety;
2. validate schema;
3. read manifest;
4. check compatibility;
5. load reference;
6. load selected targets;
7. apply recorded transformations;
8. create StructLens namespaced selections;
9. apply default visualization preset;
10. populate plugin GUI;
11. report any unsupported optional feature without failing the entire analysis.

---

# 36. PyMOL Object Naming

All plugin-created objects/selections use a namespace.

Example:

```text
SL_<analysis>_REF
SL_<analysis>_T001
SL_<analysis>_T002
SL_<analysis>_mutations
SL_<analysis>_key_residues
SL_<analysis>_mutated_key_residues
SL_<analysis>_outliers
SL_<analysis>_insertions
SL_<analysis>_deletions
```

Per-target selections:

```text
SL_<analysis>_T001_mutations
SL_<analysis>_T001_key_residues
SL_<analysis>_T001_outliers
```

Sanitize identifiers before inserting them into PyMOL names.

Never overwrite unrelated user objects.

---

# 37. Plugin Visualization Modes

The plugin must allow explicit separation of:

```text
what to highlight
```

and:

```text
how to show it
```

Highlight filters:

- mutations;
- key residues;
- mutated key residues;
- outliers;
- insertions/deletions;
- high Cα displacement;
- high backbone RMSD;
- selected residue;
- selected local environment.

Representations:

- sticks;
- spheres;
- sticks + spheres;
- cartoon + sticks;
- surface patch;
- labels.

Color modes:

- reference vs target;
- mutation status;
- physicochemical class;
- structural displacement;
- backbone RMSD;
- side-chain RMSD;
- outlier;
- key-residue group.

---

# 38. Plugin GUI

Keep the plugin intentionally smaller than the main StructLens GUI.

Suggested layout:

```text
STRUCTLENS-PYMOL
────────────────────────────────

Analysis
[ analysis name ]

Target
[ SHV-1 ▼ ]

Display
☑ Reference
☑ Target

Highlight
[ Mutations ▼ ]

Representation
[ Sticks ▼ ]

Color by
[ Mutation status ▼ ]

☑ Labels
☐ Local environment
Radius [ 5.0 Å ]

[ Focus selection ]
[ Apply ]
[ Reset StructLens view ]
```

Secondary sections:

```text
Key residues
Mutations
Outliers
Help
```

The plugin must make clear:

```text
Scientific calculations were generated by StructLens.
This plugin visualizes the imported analysis.
```

---

# 39. Residue-of-Interest Differentiation in PyMOL

When the bundle includes key residues, the plugin must create distinct categories.

At minimum:

```text
key residue, conserved
key residue, mutated
mutation outside key site
outlier
insertion/deletion context
```

Do not rely only on color.

Combine:

- color;
- representation;
- selection name;
- optional label.

Example preset:

```text
Reference = cartoon
Target = cartoon
Conserved key residues = sticks
Mutated key residues = sticks + spheres + labels
Other mutations = sticks
Outliers = optional separate selection
```

---

# 40. Multiple Targets in PyMOL

The plugin must support bundles with more than one target.

Target control:

```text
Active target
[ SHV-1 ▼ ]
```

Display modes:

```text
Active target only
Selected targets
All targets
```

When all targets are shown, prevent visual chaos by default:

- reference remains visually dominant;
- inactive targets use restrained representations;
- key/mutation highlights apply to the active target unless user explicitly requests all.

---

# 41. PyMOL Presets

Ship plugin presets:

```text
Overview
Mutation Focus
Key Residues
Mutated Key Residues
Structural Deviation
Outliers
Publication
```

`Publication` should:

- hide unnecessary labels;
- use clean cartoon/stick representations;
- choose a publication-appropriate background;
- retain explicit user override;
- not modify scientific data.

---

# 42. Plugin Publication Image Export

Even though StructLens main already exports publication figures, the plugin must allow the user to export the **current PyMOL 3D scene**.

Required:

```text
JPEG
TIFF
```

Presets:

```text
300 dpi
600 dpi
```

Use actual PyMOL ray/render capabilities where appropriate.

Do not label a screenshot upscale as 600 dpi.

Expose:

- pixel dimensions;
- physical size;
- background;
- ray tracing on/off;
- labels on/off.

---

# 43. Plugin Installation Documentation

`structlens-pymol/README.md` must include:

```text
Installation
```

with the GitHub Release workflow.

Required user instructions:

```text
1. Open the Releases page of amgoncalvesusp/pymol-plugins.
2. Open the latest release named “StructLens-PyMOL”.
3. Download StructLens-PyMOL-vX.Y.Z.zip.
4. Open PyMOL.
5. Open Plugin → Plugin Manager.
6. Choose Install New Plugin.
7. Select the downloaded ZIP.
8. Restart PyMOL if requested.
9. Open Plugin → StructLens-PyMOL.
```

Also document:

```text
Open StructLens analysis
```

workflow:

```text
Plugin → StructLens-PyMOL → Open Analysis
```

and select:

```text
*.structlens-pymol
```

---

# 44. Root `pymol-plugins` README Update

Add StructLens-PyMOL to the existing plugin table.

Example description:

```text
StructLens-PyMOL | Visualizes StructLens sequence/structure comparison results in PyMOL, including mapped residues, mutations, key residues, outliers, and structural deviations.
```

Keep each plugin's detailed installation inside its own directory README.

---

# 45. Build the Plugin ZIP Reproducibly

Create:

```text
structlens-pymol/scripts/build_release_zip.py
```

Input:

```text
plugin source tree
version
```

Output:

```text
dist/StructLens-PyMOL-vX.Y.Z.zip
```

Build requirements:

- deterministic file ordering;
- no `__pycache__`;
- no `.pyc`;
- no tests;
- no local virtualenv;
- no Git metadata;
- include plugin package;
- include needed icon/assets;
- include minimal license/readme metadata required for installation.

Test the ZIP by extracting it to a temporary directory and validating expected entrypoints.

---

# 46. Plugin Version Metadata

Expose:

```python
__version__ = "0.1.0"
SUPPORTED_BUNDLE_SCHEMA_MAJOR = 1
```

The plugin About/Help section displays:

```text
StructLens-PyMOL 0.1.0
Bundle schema 1.x
```

---

# 47. Release Automation

Prefer a dedicated GitHub Actions workflow scoped to tags:

```text
structlens-pymol-v*
```

Suggested workflow:

```text
tag pushed
    ↓
run plugin unit tests
    ↓
build plugin ZIP
    ↓
validate ZIP
    ↓
calculate SHA-256
    ↓
create GitHub Release
    ↓
attach ZIP
    ↓
attach checksum
```

Do not trigger the StructLens-PyMOL release workflow for tags belonging to other plugins.

---

# 48. Main StructLens Export Page

Purpose text:

```text
Export scientific data, publication figures, project files, or PyMOL visualization bundles.
```

Required sections:

```text
Data
    XLSX
    CSV
    TSV
    JSON

Figures
    JPEG
    TIFF

Project
    .structlens.json

PyMOL
    .structlens-pymol
```

Each export option explains:

- what is included;
- appropriate use;
- units/resolution;
- whether the output is re-importable.

---

# 49. Main StructLens XLSX Expansion for Multi-Structure Analysis

For multi-target projects, workbook sheets include as applicable:

```text
Summary
Target Summary
Pairwise Matrix
Sequence Metrics
Structure Metrics
Residue Mapping
Mutation Matrix
Key Residues
Outliers
Per-Residue Metrics
Structural Variability
Settings
Provenance
```

Use real numeric cells.

Units stay in headers.

Do not store heatmap screenshots in place of scientific matrix values.

---

# 50. Data Model Additions

Create stable multi-analysis types.

```python
@dataclass(slots=True)
class TargetAnalysis:
    target_id: str
    correspondence: tuple[ResidueCorrespondence, ...]
    mutations: tuple[MutationEvent, ...]
    sequence_metrics: SequenceMetrics
    structural_metrics: StructuralMetrics
    transform: StructuralTransform
```

```python
@dataclass(slots=True)
class ReferenceVsManyAnalysis:
    reference_id: str
    targets: dict[str, TargetAnalysis]
```

```python
@dataclass(slots=True)
class PairwiseMatrix:
    metric_name: str
    structure_ids: tuple[str, ...]
    values: dict[tuple[str, str], float | None]
    unit: str | None
```

```python
@dataclass(slots=True)
class MultipleStructureAnalysis:
    structure_ids: tuple[str, ...]
    aligned_positions: tuple[MultiStructurePosition, ...]
    transforms: dict[str, StructuralTransform]
```

---

# 51. Shared Selection State in Main App

Create a single application-level selection concept.

```python
@dataclass(slots=True)
class AnalysisSelection:
    reference_residue: ResidueId | None
    target_id: str | None
    target_residue: ResidueId | None
    pair: tuple[str, str] | None
    key_group: str | None
```

Pages subscribe to selection changes.

This allows:

```text
click chart
→ residue selected
→ Residues page updates
→ PyMOL export preset knows selected residue
```

without direct page-to-page coupling.

---

# 52. Error Handling

New typed errors:

```python
class BundledBackendUnavailableError(StructLensError): ...
class UnsupportedPlatformError(StructLensError): ...
class PyMOLNotConfiguredError(StructLensError): ...
class PyMOLPluginUnavailableError(StructLensError): ...
class BundleValidationError(StructLensError): ...
class BundleCompatibilityError(StructLensError): ...
class UnsafeBundleError(StructLensError): ...
class MultiStructureAlignmentError(StructLensError): ...
```

Do not show users stack traces by default.

GUI messages must be actionable.

Bad:

```text
USalign executable was not found.
```

Good in a packaged release:

```text
The bundled structural-alignment backend is unavailable for this installation.
Open Diagnostics to verify the StructLens installation or choose a custom
US-align executable in Advanced Settings.
```

---

# 53. Main StructLens Implementation Tasks

## Task S1 — Freeze new architecture

Create/update:

```text
docs/architecture.md
docs/pymol-integration.md
docs/multi-structure-analysis.md
```

Document:

- desktop/plugin separation;
- bundle contract;
- comparison modes;
- chart architecture;
- US-align bundling;
- plugin release repository.

Commit:

```text
docs: define StructLens multi-analysis and PyMOL architecture
```

---

## Task S2 — Bundle US-align

TDD first.

Test:

- platform detection;
- bundled executable discovery;
- executable permission handling;
- custom override;
- provenance;
- failure message;
- gzip temporary normalization;
- cleanup.

Acceptance:

```text
fresh packaged StructLens install
→ Structure mode
→ no PATH configuration
→ US-align runs
```

Commit:

```text
feat: bundle US-align structural backend
```

---

## Task S3 — Multi-structure domain model

Implement:

- `ComparisonMode`;
- `TargetAnalysis`;
- `ReferenceVsManyAnalysis`;
- `PairwiseMatrix`;
- `MultipleStructureAnalysis`.

Tests:

- deterministic IDs;
- round-trip project serialization;
- no result overwriting between targets.

Commit:

```text
feat: add multi-structure analysis model
```

---

## Task S4 — Reference vs Many

Tests:

- one reference + 3 targets;
- independent transformations;
- independent mutation sets;
- reference numbering preserved;
- key-residue mapping.

Commit:

```text
feat: compare one reference against multiple proteins
```

---

## Task S5 — All vs All

Tests:

- correct pair count;
- no duplicate calculations;
- symmetric matrix output;
- missing/failed pair represented explicitly.

Commit:

```text
feat: add all-vs-all sequence and structure matrices
```

---

## Task S6 — Multiple Structure Alignment

Integrate US-align multiple-structure mode.

Tests:

- parse multi-structure output fixture;
- preserve structure IDs;
- map aligned columns back to original residues;
- compute coverage;
- compute Cα positional variability.

Commit:

```text
feat: add multiple structural alignment
```

---

## Task S7 — Sequence area

Refactor GUI/application services so sequence analysis is visible as its own scientific area.

Tests:

- reference-vs-many mutation matrix;
- all-vs-all sequence identity;
- residue selection propagation.

Commit:

```text
feat: promote sequence analysis to first-class workflow
```

---

## Task S8 — New navigation

Implement:

```text
Project
Sequences
Structures
Residues
Charts
PyMOL
Results
Export
```

Apply Impeccable first if available.

Commit:

```text
refactor: reorganize StructLens scientific workflow
```

---

## Task S9 — Charts foundation

Create chart data models independent from widgets.

Charts consume application data and return explicit plot datasets.

Do not let chart widgets recalculate scientific values.

Commit:

```text
feat: add reusable scientific chart data layer
```

---

## Task S10 — Six required charts

Implement and test:

1. structural deviation profile;
2. mutation/conservation matrix;
3. pairwise similarity heatmap;
4. sequence–structure scatter;
5. structural conservation profile;
6. key-residue comparison.

Acceptance:

- point/cell selection updates shared selection;
- units present;
- underlying values exportable;
- no invented values.

Commit:

```text
feat: add comparative sequence-structure charts
```

---

## Task S11 — Chart publication export

Implement:

```text
XLSX
JPEG
TIFF
300 dpi
600 dpi
```

Tests:

- exact data round-trip;
- dimensions;
- DPI metadata where supported;
- readable files;
- no screenshot upscaling.

Commit:

```text
feat: export publication-ready comparative charts
```

---

## Task S12 — Bundle schema

Create JSON schema/documented equivalent for:

```text
.structlens-pymol schema 1.0
```

Add:

- valid fixtures;
- invalid fixtures;
- security fixtures.

Commit:

```text
feat: define StructLens-PyMOL interchange schema
```

---

## Task S13 — Bundle writer

Implement:

```python
write_pymol_bundle(...)
validate_pymol_bundle(...)
```

Security:

- atomic output;
- deterministic manifest;
- no executable payload;
- no unsafe paths.

Commit:

```text
feat: export validated PyMOL analysis bundles
```

---

## Task S14 — PyMOL page

Implement:

```text
Open in PyMOL
Export for PyMOL
Plugin installation instructions
```

Use explicit integration status.

Commit:

```text
feat: add explicit PyMOL integration workflow
```

---

## Task S15 — External PyMOL launcher

Implement `PyMOLLauncher`.

Test on each supported OS in CI/manual matrix where feasible.

Acceptance is end-to-end launch behavior.

Commit:

```text
feat: open StructLens analyses in external PyMOL
```

---

# 54. StructLens-PyMOL Implementation Tasks

All plugin implementation happens in:

```text
amgoncalvesusp/pymol-plugins
```

under:

```text
structlens-pymol/
```

Use a separate worktree/branch.

Suggested branch:

```text
feat/structlens-pymol
```

---

## Task P1 — Plugin skeleton

Create:

```text
structlens-pymol/
```

with package, README, CHANGELOG, tests, assets, build script.

Add root README row.

Commit:

```text
feat(structlens-pymol): add plugin skeleton
```

---

## Task P2 — Bundle validator and reader

No PyMOL import in the core bundle parser.

Test ZIP traversal attacks such as:

```text
../../evil.py
```

Test malformed manifest.

Test unsupported schema major.

Commit:

```text
feat(structlens-pymol): read validated StructLens bundles
```

---

## Task P3 — PyMOL adapter

Implement object loading, transformations, safe namespaced selections.

Tests where possible with mocked `cmd`.

Commit:

```text
feat(structlens-pymol): load analysis structures into PyMOL
```

---

## Task P4 — Semantic selections

Create:

- mutations;
- key residues;
- mutated key residues;
- outliers;
- insertions/deletions;
- selected residue;
- local environment.

Commit:

```text
feat(structlens-pymol): create residue-analysis selections
```

---

## Task P5 — Plugin GUI

Apply Impeccable if available.

Implement the compact visualization client described above.

Commit:

```text
feat(structlens-pymol): add interactive visualization panel
```

---

## Task P6 — Visualization renderer

Implement highlight/color/representation independence.

Presets:

- Overview;
- Mutation Focus;
- Key Residues;
- Mutated Key Residues;
- Structural Deviation;
- Outliers;
- Publication.

Commit:

```text
feat(structlens-pymol): visualize comparative residue states
```

---

## Task P7 — Multi-target support

Test:

- target dropdown;
- active target;
- multiple displayed targets;
- target-specific mutations/selections.

Commit:

```text
feat(structlens-pymol): support multi-target analyses
```

---

## Task P8 — Publication export

Implement ray/render-driven JPEG/TIFF export.

Test dimensions.

Commit:

```text
feat(structlens-pymol): export publication-quality PyMOL scenes
```

---

## Task P9 — PyMOL command API

Register a stable command callable by launch handoff.

Conceptual API:

```text
structlens_open <path-to-bundle>
```

The exact PyMOL registration mechanism must follow supported PyMOL APIs.

This command:

1. validates the path;
2. opens the bundle;
3. activates the plugin;
4. focuses the default scene.

Commit:

```text
feat(structlens-pymol): add external bundle-open command
```

---

## Task P10 — Installable ZIP

Build:

```text
StructLens-PyMOL-v0.1.0.zip
```

Perform clean PyMOL Plugin Manager install test.

Test:

```text
fresh PyMOL profile
→ install ZIP
→ restart if needed
→ Plugin → StructLens-PyMOL
→ Open Analysis
→ fixture bundle
→ structures visible
```

Commit:

```text
build(structlens-pymol): add reproducible release package
```

---

## Task P11 — Documentation

README must include:

- purpose;
- screenshots;
- GitHub Release installation;
- supported PyMOL versions;
- `.structlens-pymol` workflow;
- highlight modes;
- multi-target workflow;
- image export;
- troubleshooting;
- uninstall/update.

Commit:

```text
docs(structlens-pymol): add installation and usage guide
```

---

## Task P12 — Release workflow

Create CI/release automation for:

```text
structlens-pymol-v*
```

Verify it does not trigger for other plugin tags.

Commit:

```text
ci(structlens-pymol): add dedicated release workflow
```

---

# 55. StructLens-PyMOL v0.1.0 Release Gate

Do not create the release until all checks pass.

Required:

- [ ] package unit tests pass;
- [ ] bundle security tests pass;
- [ ] schema compatibility tests pass;
- [ ] Plugin Manager ZIP install succeeds;
- [ ] plugin opens;
- [ ] pairwise bundle opens;
- [ ] Reference-vs-Many bundle opens;
- [ ] key residues show;
- [ ] mutated key residues are distinguishable;
- [ ] outliers show;
- [ ] active target switching works;
- [ ] reset does not destroy unrelated PyMOL objects;
- [ ] publication JPEG export works;
- [ ] publication TIFF export works;
- [ ] root repository README lists the plugin;
- [ ] plugin README installation instructions are verified;
- [ ] CHANGELOG is updated;
- [ ] tag is `structlens-pymol-v0.1.0`;
- [ ] release asset is `StructLens-PyMOL-v0.1.0.zip`;
- [ ] checksum asset is generated;
- [ ] release notes include tested PyMOL versions.

Only then publish:

```text
StructLens-PyMOL v0.1.0
```

---

# 56. End-to-End Scientific Acceptance Workflow

Use a fixture with:

```text
1 reference
3 targets
known substitutions
at least 1 key residue
at least 1 mutated key residue
at least 1 structural outlier
```

Run:

```text
StructLens
→ Project
→ Reference vs Many
→ Sequences
→ Structures
→ Residues
→ Charts
→ Export for PyMOL
```

Verify:

```text
.structlens-pymol bundle created
```

Then:

```text
PyMOL
→ Plugin
→ StructLens-PyMOL
→ Open Analysis
```

Verify:

- 4 structures load;
- correct transforms apply;
- target selector contains 3 targets;
- mutation counts match StructLens;
- key-residue selections match;
- mutated key residue is visibly distinguishable;
- outlier selection matches;
- residue identifiers preserve chains and insertion codes;
- no scientific metric is recomputed differently in the plugin.

Finally:

```text
StructLens → Open in PyMOL
```

Verify the same analysis opens automatically.

---

# 57. Scientific Regression Requirements

Add regression fixtures for:

- pairwise sequence/structure analysis;
- reference-vs-many;
- all-vs-all matrix;
- multiple structural alignment;
- mutation matrix;
- structural variability profile;
- key-residue comparison;
- `.structlens-pymol` bundle;
- PyMOL object/selection manifest.

Never approve a numerical change solely because the UI still works.

---

# 58. GUI Usability Requirements

The entire user-facing application remains English-only.

Every mode explains:

1. what it does;
2. when to use it;
3. what it needs;
4. what it produces;
5. important limitations.

Specific contextual help is mandatory for:

- Pairwise;
- Reference vs Many;
- All vs All;
- Multiple Structure Alignment;
- Auto;
- Sequence-guided;
- Structure-guided;
- Manual mapping;
- sequence identity;
- coverage;
- TM-score;
- strict RMSD;
- refined RMSD;
- Cα displacement;
- backbone RMSD;
- side-chain RMSD;
- structural variability;
- BLOSUM62;
- Grantham;
- key residues;
- PyMOL bundle;
- 300/600 dpi export.

The user should be able to complete a routine analysis without external documentation.

---

# 59. Explicit Non-Goals for This Revision

Do not implement without a separate approved plan:

- continuous live StructLens ↔ PyMOL synchronization;
- socket/HTTP/WebSocket bridge;
- full MD trajectory analysis;
- FEP integration;
- functional/pathogenic mutation prediction;
- automatic energetic interpretation;
- web/cloud backend;
- true general-purpose multiple sequence alignment dependency;
- DockLens integration;
- network interaction analysis.

These can be later extensions.

---

# 60. Recommended Subagent Allocation

After contracts are frozen:

### Luna A — Bundled backend

- US-align resources;
- platform discovery;
- provenance;
- tests.

### Luna B — Multi-structure core

- Reference vs Many;
- All vs All;
- multiple structural alignment.

### Luna C — Charts

- chart data models;
- six chart views;
- publication export.

### Luna D — Interchange

- `.structlens-pymol` schema;
- writer;
- validator;
- security.

### Luna E — PyMOL plugin core

Work in `amgoncalvesusp/pymol-plugins`:

- bundle reader;
- PyMOL adapter;
- semantic selections.

### Luna F — PyMOL plugin GUI

Only after plugin core contracts stabilize.

### Luna G — Validation

- golden data;
- end-to-end fixtures;
- cross-product consistency checks.

Main Luna owns:

- interface contracts;
- shared selection model;
- integration;
- release gates.

---

# 61. Review Gates

## Gate A — Main core

Must pass before GUI work:

- bundled US-align;
- new domain model;
- Reference vs Many;
- All vs All;
- MSTA;
- schema contract.

## Gate B — Main GUI

Must pass:

- new navigation;
- sequence/structure separation;
- charts;
- PyMOL page;
- export page.

## Gate C — Plugin core

Must pass:

- bundle read;
- safe extraction;
- structure load;
- transforms;
- semantic selections.

## Gate D — Plugin GUI

Must pass:

- target switching;
- highlight modes;
- presets;
- reset;
- publication export.

## Gate E — Release

Must pass:

- end-to-end StructLens → bundle → PyMOL;
- automatic Open in PyMOL;
- clean install from GitHub Release ZIP;
- scientific regression.

---

# 62. Completion Verification

Main StructLens:

```bash
pytest -q
ruff check .
mypy src/structlens
```

Plugin repository:

```bash
pytest -q structlens-pymol/tests
```

Build release ZIP:

```bash
python structlens-pymol/scripts/build_release_zip.py --version 0.1.0
```

Expected:

```text
structlens-pymol/dist/StructLens-PyMOL-v0.1.0.zip
```

Perform a real clean PyMOL installation test before tagging.

Use:

```text
superpowers:verification-before-completion
```

before declaring either repository ready.

---

# 63. Initial Directive for Luna Extra High

Use this as the execution prompt:

```text
Implement the approved StructLens v0.2 architectural revision according to
StructLens_v0.2_Architecture_and_Implementation_Plan.md.

Read the entire document before touching code. Use Superpowers. Start in isolated
git worktrees and use test-driven development.

The main StructLens application and StructLens-PyMOL are separate products.
StructLens owns scientific analysis. StructLens-PyMOL visualizes a versioned
.structlens-pymol data bundle and must not silently recompute or reinterpret the
scientific results.

Bundle US-align with StructLens so normal users never need to configure PATH.
Implement Pairwise, Reference vs Many, All vs All, and Multiple Structure
Alignment modes. Promote sequence analysis to a first-class workflow. Implement
the six required comparative charts and publication-ready XLSX/JPEG/TIFF exports.

Implement the PyMOL plugin in the existing GitHub repository:
https://github.com/amgoncalvesusp/pymol-plugins
under the new top-level directory:
structlens-pymol/

Do not create a separate GitHub repository for the plugin.

StructLens-PyMOL must have independent semantic versioning inside the monorepo.
Use tags named structlens-pymol-vMAJOR.MINOR.PATCH and create a dedicated GitHub
Release named StructLens-PyMOL vMAJOR.MINOR.PATCH with the installable asset
StructLens-PyMOL-vMAJOR.MINOR.PATCH.zip.

Before changing either GUI, discover and use the Impeccable skill if available,
then use it again for a final GUI critique/polish pass.

Do not implement live socket/HTTP/WebSocket synchronization in this milestone.
The approved integration is a validated .structlens-pymol bundle plus an optional
Open in PyMOL launch handoff.

Do not claim completion until the end-to-end test succeeds:
StructLens analysis → .structlens-pymol bundle → installable StructLens-PyMOL
plugin → PyMOL visualization with targets, mutations, key residues, mutated key
residues, outliers, and recorded structural transforms.
```
