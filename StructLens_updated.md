# StructLens Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: use `superpowers:subagent-driven-development` (recommended) or `superpowers:executing-plans` to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.
>
> **Primary implementation model:** Luna Extra High.
>
> **GUI requirement:** before creating or revising the graphical interface, discover whether the `impeccable` skill is installed in the execution environment. If it exists, it must be invoked for GUI architecture/design and again for GUI critique/polish. If it does not exist, record that fact in the implementation log and follow the GUI specification in this document. Never pretend that an unavailable skill was used.

**Goal:** Build StructLens, a standalone Python structural-bioinformatics engine plus a PyMOL plugin for reproducible sequence/structure alignment, explicit residue mapping, mutation detection, residue-by-residue structural comparison, and interactive highlighting of biologically relevant residues.

**Architecture:** The project is split into a PyMOL-independent scientific core, an application/service layer, external adapters, a CLI, and a thin PyMOL GUI/visualization integration. The authoritative scientific object is an explicit `ResidueCorrespondence` map; PyMOL is a visualization backend, never the scientific source of truth. Sequence and structure alignment may propose mappings, but all downstream calculations must use the persisted correspondence table.

**Tech Stack:** Python 3.11+, NumPy, Biopython, Qt bindings available inside the target PyMOL distribution, pytest, hypothesis where useful, Ruff, mypy, US-align as an external structural-alignment backend, `openpyxl` for publication-ready XLSX workbooks, Pillow and/or PyMOL rendering APIs for raster export, and standard-library `json`, `csv`, `hashlib`, and `subprocess`.

**Spec:** This file is the authoritative implementation spec and execution plan.

---

## 0. Required Superpowers Workflow

Before changing code, Luna must use the following workflow.

1. `superpowers:brainstorming` only to surface contradictions, missing scientific assumptions, or major architectural risks. The product decisions in this file remain authoritative unless a contradiction is demonstrated.
2. `superpowers:using-git-worktrees` before implementation. Work in an isolated worktree.
3. `superpowers:writing-plans` may refine task granularity, but must not weaken the scientific requirements in this document.
4. `superpowers:test-driven-development` for every feature and bug fix.
5. `superpowers:subagent-driven-development` for independent tasks where contracts are already frozen. Prefer this over broad parallel edits to the same modules.
6. `superpowers:systematic-debugging` for failing tests, PyMOL integration failures, numerical discrepancies, or unexplained GUI behavior.
7. `superpowers:requesting-code-review` after major milestones.
8. `superpowers:verification-before-completion` before claiming the project, milestone, feature, or release is complete.
9. `superpowers:finishing-a-development-branch` only after all verification gates pass.

Do not begin GUI implementation until the domain interfaces and application-service contracts are stable.

---


## 0.1 Additional Global Constraints

- User-facing language: English only.
- Approved StructLens icon must be packaged and used consistently.
- Every scientific table exposed by the product must support XLSX export.
- Every StructLens-generated molecular/figure view exposed by the product must support publication-quality JPEG/TIFF export when technically possible.
- Image export must provide 300 dpi and 600 dpi publication presets.
- Every workflow option and advanced metric must have contextual English help.
- Usability explanations are release requirements and must be tested/reviewed.

---

## 1. Product Definition

StructLens answers four scientific questions:

1. **Which residues in two or more protein structures are biologically equivalent?**
2. **Which equivalent positions are conserved, substituted, inserted, deleted, or non-standard?**
3. **How different are the aligned structures globally and at each mapped residue?**
4. **How can those differences be explored immediately and reproducibly inside PyMOL?**

StructLens is not merely a wrapper around `PyMOL align`, `super`, or `cealign`. It owns its own residue mapping, transformation, metrics, mutation model, project state, exports, and audit trail.


---

## 1.1 Product Identity, Language, and Usability Requirements

### Approved icon

StructLens has an approved application icon supplied with the project. The implementation must package and use that icon consistently in:

- PyMOL plugin/menu entry where supported;
- main StructLens panel/window;
- desktop/application bundle when one is produced;
- documentation;
- release assets;
- exported report metadata or cover sheet only when appropriate.

Store the packaged application icon under:

```text
src/structlens/plugin/assets/structlens_icon.png
```

If platform packaging needs derived formats such as `.ico` or `.icns`, generate those during packaging from the approved master asset rather than redesigning the icon.

Do not use PyMOL branding inside the StructLens icon and do not substitute the approved visual identity with a generic molecule icon.

### Language

StructLens v0.1 is **English-only**.

All of the following must be written in English:

- GUI labels;
- tooltips;
- workflow explanations;
- validation messages;
- errors;
- legends;
- column headers;
- sheet names in XLSX files;
- image annotations;
- CLI output;
- documentation;
- example projects;
- exported metadata.

Do not create a localization/i18n framework in v0.1 unless it is required by the chosen Qt stack. If the framework exists implicitly, ship only English strings.

### Intuitive, self-explanatory workflows

The user must not need to understand internal alignment algorithms before starting an analysis.

Every selectable workflow or analysis mode must include an immediately accessible explanation covering:

1. **What it does**
2. **When to use it**
3. **What data it needs**
4. **What it produces**
5. **Important limitations or caveats**

For compact controls, use a combination of:

- concise inline helper text;
- tooltips;
- information icons;
- expandable “What does this mean?” or “When should I use this?” help;
- contextual empty-state guidance.

The interface must prefer plain scientific language over implementation terminology.

Example for `AUTO`:

```text
Auto
Recommended for most comparisons.

StructLens first evaluates sequence similarity and coverage, then chooses
sequence-guided or structure-guided residue mapping. The selected method and
reason are always reported in Results.
```

Example for `Sequence`:

```text
Sequence
Best for homologous proteins with meaningful sequence similarity.

Residue equivalence is established from a global amino-acid sequence alignment
before structural superposition.
```

Example for `Structure`:

```text
Structure
Best when sequence identity is low or when fold similarity is more informative.

US-align is used to propose structurally equivalent residues. The resulting
correspondence remains explicit and inspectable.
```

Example for `Manual`:

```text
Manual
Best when specific biologically validated residue correspondences must be
enforced.

You define residue pairs explicitly. Locked pairs are never overwritten by
automatic mapping.
```

Workflow explanations are part of the product requirement, not optional documentation.

---

## 2. Scientific Invariants

These invariants are mandatory throughout the codebase.

### 2.1 Residue identity is not residue number

Never use only a PDB residue number such as `166` as a unique identity.

A residue must be identified by a structured locator containing at least:

```python
@dataclass(frozen=True, slots=True)
class ResidueId:
    structure_id: str
    model_id: str
    chain_id: str
    auth_seq_id: str
    insertion_code: str | None
    residue_name: str
```

For mmCIF, retain both author and label numbering when available:

```python
@dataclass(frozen=True, slots=True)
class ResidueNumbering:
    auth_seq_id: str
    label_seq_id: str | None
    insertion_code: str | None
```

### 2.2 Residue mapping is authoritative

All downstream calculations consume `ResidueCorrespondence` objects.

```python
@dataclass(slots=True)
class ResidueCorrespondence:
    alignment_index: int
    reference: ResidueId | None
    target: ResidueId | None
    reference_one_letter: str | None
    target_one_letter: str | None
    status: str
    sequence_score: float | None = None
    ca_displacement_angstrom: float | None = None
    backbone_rmsd_angstrom: float | None = None
    sidechain_rmsd_angstrom: float | None = None
    all_heavy_atom_rmsd_angstrom: float | None = None
    is_outlier: bool = False
    is_key_residue: bool = False
    mapping_source: str = "unknown"
    mapping_locked: bool = False
```

Allowed `status` values must be represented by an enum:

- `CONSERVED`
- `SUBSTITUTION`
- `INSERTION`
- `DELETION`
- `NONSTANDARD`
- `UNMAPPED`

### 2.3 A single Cα pair does not produce a residue RMSD

For one mapped residue pair:

- `Cα displacement` = Euclidean distance between aligned Cα coordinates.
- `Backbone RMSD` = RMSD over matched N, Cα, C, O atoms.
- `Side-chain RMSD` = RMSD over matched side-chain heavy atoms, with symmetry-aware matching.
- `All-heavy-atom RMSD` = RMSD over all matched heavy atoms.

The GUI, CLI, exports, docs, and tests must use these terms exactly.

### 2.4 Strict and refined RMSD are separate outputs

Never silently discard outliers.

Every global analysis must preserve:

- strict RMSD;
- strict atom/residue count;
- refined RMSD, if refinement is enabled;
- refined atom/residue count;
- explicit list of excluded residue correspondences;
- refinement cutoff and algorithm.

### 2.5 PyMOL is never the scientific state

The core must not import `pymol`.

Only `src/structlens/integrations/pymol/` and `src/structlens/plugin/` may import PyMOL APIs.

### 2.6 Units are explicit

All distances and RMSDs are stored and displayed in Å. Internal variable names that contain distance or RMSD should end in `_angstrom` unless the type makes the unit unambiguous.

### 2.7 No functional overclaiming

BLOSUM62 scores, Grantham distances, structural displacement, local RMSD, and mutation class are descriptors. StructLens must not report that a mutation is pathogenic, stabilizing, destabilizing, catalytic, or function-altering unless such a future feature explicitly implements and validates that inference.

---

## 3. Supported Workflows

### Workflow A — PyMOL objects

1. User opens structures in PyMOL.
2. StructLens lists available molecular objects and chains.
3. User chooses one reference and one or more targets.
4. StructLens analyzes without reloading the coordinates unless necessary.

### Workflow B — Files

Initial formats:

- PDB
- mmCIF / CIF
- gzip-compressed PDB
- gzip-compressed mmCIF/CIF
- FASTA / FA for sequence-only mapping support

### Workflow C — Mixed input

Reference may come from a PyMOL object while targets come from files, or vice versa. Normalize both into the same domain model before analysis.

### Workflow D — Manual locked mapping

The user may explicitly map residue pairs and lock them. Locked residue pairs cannot be overwritten by automatic mapping.

Example:

```text
A:70  <-> B:68
A:73  <-> B:71
A:130 <-> B:128
A:166 <-> B:164
```

---

## 4. Alignment Modes

Implement these modes as a stable enum.

```python
class AlignmentMode(str, Enum):
    AUTO = "auto"
    SEQUENCE = "sequence"
    STRUCTURE = "structure"
    MANUAL = "manual"
```

### AUTO

AUTO is an explicit policy, not a hidden heuristic.

Version 0.1 policy:

1. Obtain sequence mapping.
2. If sequence identity across aligned canonical residues is at least `0.30` and coverage is at least `0.70`, use sequence-guided correspondence as the initial map.
3. Perform structural superposition from the resulting mapped Cα pairs.
4. Validate geometrically and mark structural outliers without deleting correspondences from the authoritative map.
5. If identity or coverage falls below the thresholds, use US-align as the initial structural mapping source.
6. Report in `AnalysisResult.alignment_decision` which branch was used and why.

Thresholds must be configurable in `AnalysisSettings`; these are defaults, not hard-coded hidden constants.

### SEQUENCE

Use global pairwise sequence alignment for pairwise analyses. Default substitution matrix: BLOSUM62. Keep gap penalties explicit and serialized in project state.

### STRUCTURE

Use the US-align adapter. Parse correspondence output into the StructLens domain model. Do not treat US-align stdout as persistent state.

### MANUAL

Accept only explicit mappings supplied by the application layer. The core validates locators and then superposes on the locked pairs that contain required atoms.

---

## 5. Mutation Model

Mutation detection runs only after residue mapping.

```python
class MutationKind(str, Enum):
    CONSERVED = "conserved"
    SUBSTITUTION = "substitution"
    INSERTION = "insertion"
    DELETION = "deletion"
    NONSTANDARD = "nonstandard"
```

```python
@dataclass(frozen=True, slots=True)
class MutationEvent:
    alignment_index: int
    kind: MutationKind
    reference: ResidueId | None
    target: ResidueId | None
    reference_aa: str | None
    target_aa: str | None
    reference_label: str
    target_label: str
    canonical_notation: str
    blosum62_score: int | None
    grantham_distance: int | None
    physicochemical_class: str | None
```

Canonical notation uses reference numbering whenever a canonical reference residue exists.

Examples:

- `S130T`
- `E166A`
- `ΔGly145` for a deletion
- `ins145_G` only when a clear reference-anchored insertion notation can be generated reproducibly; otherwise display a structured insertion label rather than inventing ambiguous human notation.

Mutation classification must distinguish conservative/non-conservative descriptors from functional claims.

---

## 6. Visualization Model

Separate four concepts.

```text
FILTER LAYER
    what residues are included

DATA LAYER
    what value/category drives styling

COLOR LAYER
    color mapping

REPRESENTATION LAYER
    sticks / spheres / cartoon / surface / labels
```

### Highlight filters

Initial filters:

- all mapped residues;
- mutations only;
- conserved only;
- key residues only;
- mutated key residues;
- outliers;
- insertions/deletions;
- residues above Cα displacement threshold;
- residues above backbone RMSD threshold.

### Color modes

Initial color modes:

- fixed reference vs target colors;
- mutation status;
- physicochemical mutation class;
- BLOSUM62 category;
- Grantham-distance category;
- Cα displacement continuous scale;
- backbone RMSD continuous scale;
- side-chain RMSD continuous scale;
- key-residue category;
- outlier status.

Continuous scales must always render a visible legend with units.

### Representations

- sticks;
- spheres;
- sticks + spheres;
- cartoon + selected sticks;
- surface patch for selected residues;
- labels;
- optional distance lines for selected Cα pairs.

### Presets

Ship these named presets:

- `Minimal`
- `Publication`
- `Mutation focus`
- `Structural deviation`
- `Active site`
- `Presentation`

Presets only set visualization state. They must never modify scientific analysis state.

---

## 7. GUI Specification

The GUI must be designed as a scientific tool, not a SaaS dashboard.

### Required GUI skill workflow

Before Task 14:

1. Search the execution environment for an installed `impeccable` skill.
2. If found, invoke it before producing GUI code.
3. Apply its guidance to interaction hierarchy, spacing, affordances, accessibility, state visibility, and visual polish.
4. After the GUI is functionally complete, invoke `impeccable` again for critique/polish.
5. Record both invocations in `docs/implementation-log.md`.
6. If unavailable, record `Impeccable skill unavailable in this environment` and continue with this specification.

### GUI information architecture

Use a single dockable or plugin panel with these sections:

1. **Project**
2. **Alignment**
3. **Mutations**
4. **Residues**
5. **Visualization**
6. **Results**

Avoid modal dialogs for routine navigation.

### Project

Must support:

- select reference object/file;
- select reference chain;
- select multiple target structures;
- choose each target chain;
- add files;
- add current PyMOL objects;
- remove target;
- save/open `.structlens.json` project.

### Alignment

Must support:

- Auto;
- Sequence;
- Structure;
- Manual;
- BLOSUM62 selection for the implemented sequence engine;
- gap penalties;
- refined-RMSD toggle;
- refinement cutoff;
- run/cancel analysis;
- visible chosen alignment method after AUTO resolves.

### Mutations

Must show:

- mutation count;
- searchable/filterable mutation list;
- substitutions;
- insertions;
- deletions;
- non-standard residues;
- filters for all/key residues/conservative/non-conservative.

### Residue Browser

Columns:

- reference residue;
- target residue;
- status;
- mutation notation;
- Cα displacement (Å);
- backbone RMSD (Å);
- side-chain RMSD (Å);
- outlier flag;
- key-residue group.

Selection must synchronize with PyMOL when feasible:

- clicking a row focuses the mapped residues in PyMOL;
- clicking a StructLens-owned selection in PyMOL should update the relevant residue selection in the plugin without polling aggressively.

### Visualization

Controls:

- highlight filter;
- color mode;
- representation;
- labels on/off;
- show reference/target toggles;
- local environment toggle;
- local radius in Å;
- preset selector;
- Apply;
- Reset StructLens visualization.

### Results

Show at minimum:

- sequence identity;
- alignment coverage;
- TM-score when produced by US-align;
- strict Cα RMSD;
- refined Cα RMSD when enabled;
- mapped residue count;
- refined residue count;
- outlier count;
- mutation count;
- insertion count;
- deletion count.


### Workflow guidance and contextual help

The GUI must be understandable without consulting the manual for routine tasks.

Every major page must begin with one short English sentence explaining its purpose.

Required examples:

- **Project:** “Choose the reference structure and the structures you want to compare.”
- **Alignment:** “Choose how StructLens should determine equivalent residues before superposition.”
- **Mutations:** “Review substitutions, insertions, deletions, and non-standard residues detected from the residue map.”
- **Residues:** “Inspect residue-by-residue correspondence and structural differences.”
- **Visualization:** “Choose what to highlight in PyMOL and how it should be represented.”
- **Results:** “Review global alignment quality, structural metrics, and mutation counts.”

Every non-obvious control must have a tooltip or inline helper text.

At minimum, provide contextual explanations for:

- Auto / Sequence / Structure / Manual alignment;
- sequence identity;
- alignment coverage;
- TM-score;
- strict RMSD;
- refined RMSD;
- outlier cutoff;
- Cα displacement;
- backbone RMSD;
- side-chain RMSD;
- local RMSD;
- BLOSUM62;
- Grantham distance;
- key residues;
- mutation severity/color modes;
- visualization presets;
- export resolution and image format.

Tooltips must explain the scientific meaning, not merely repeat the control label.

### GUI design constraints

- scientific density without clutter;
- no decorative gradients;
- minimal use of cards;
- large tables get priority over decorative summary panels;
- colors never carry the only meaning;
- all continuous scales have visible legends;
- all metrics show units;
- keyboard navigation must remain usable where Qt supports it;
- long-running work must not freeze the PyMOL UI;
- PyMOL-owned user objects and selections must not be destroyed by StructLens reset.

---

## 8. Repository Layout

Create this structure.

```text
structlens/
├── pyproject.toml
├── README.md
├── LICENSE
├── CHANGELOG.md
├── docs/
│   ├── architecture.md
│   ├── scientific-methods.md
│   ├── implementation-log.md
│   └── validation.md
├── src/
│   └── structlens/
│       ├── __init__.py
│       ├── core/
│       │   ├── models/
│       │   │   ├── structure.py
│       │   │   ├── residue.py
│       │   │   ├── correspondence.py
│       │   │   ├── mutation.py
│       │   │   ├── settings.py
│       │   │   └── results.py
│       │   ├── parsing/
│       │   │   ├── pdb.py
│       │   │   ├── mmcif.py
│       │   │   ├── fasta.py
│       │   │   └── normalize.py
│       │   ├── mapping/
│       │   │   ├── sequence_mapper.py
│       │   │   ├── structural_mapper.py
│       │   │   ├── manual_mapper.py
│       │   │   └── validator.py
│       │   ├── alignment/
│       │   │   ├── protocols.py
│       │   │   ├── sequence.py
│       │   │   ├── superposition.py
│       │   │   └── refinement.py
│       │   ├── geometry/
│       │   │   ├── kabsch.py
│       │   │   ├── rmsd.py
│       │   │   ├── displacement.py
│       │   │   ├── symmetry.py
│       │   │   └── dihedrals.py
│       │   ├── mutations/
│       │   │   ├── detector.py
│       │   │   ├── blosum.py
│       │   │   ├── grantham.py
│       │   │   └── physicochemistry.py
│       │   ├── metrics/
│       │   │   ├── sequence_metrics.py
│       │   │   ├── structural_metrics.py
│       │   │   └── local_metrics.py
│       │   └── validation/
│       │       └── scientific_checks.py
│       ├── application/
│       │   ├── dto.py
│       │   ├── project_state.py
│       │   ├── analysis_service.py
│       │   ├── mutation_service.py
│       │   ├── visualization_service.py
│       │   └── export_service.py
│       ├── integrations/
│       │   ├── usalign/
│       │   │   ├── adapter.py
│       │   │   ├── parser.py
│       │   │   └── executable.py
│       │   └── pymol/
│       │       ├── adapter.py
│       │       ├── selections.py
│       │       ├── transforms.py
│       │       └── state_snapshot.py
│       ├── cli/
│       │   ├── main.py
│       │   └── formatting.py
│       └── plugin/
│           ├── assets/
│           │   └── structlens_icon.png
│           ├── entrypoint.py
│           ├── controller.py
│           ├── workers.py
│           ├── gui/
│           │   ├── main_panel.py
│           │   ├── project_page.py
│           │   ├── alignment_page.py
│           │   ├── mutations_page.py
│           │   ├── residues_page.py
│           │   ├── visualization_page.py
│           │   ├── results_page.py
│           │   └── widgets.py
│           └── visualization/
│               ├── renderer.py
│               ├── presets.py
│               ├── legends.py
│               └── namespace.py
└── tests/
    ├── unit/
    ├── integration/
    ├── regression/
    ├── scientific/
    └── fixtures/
```

Keep files focused. Do not collapse multiple subsystems into one large module.

---

## 9. Frozen Interfaces Before Parallel Work

The main agent owns these interfaces. Freeze them before dispatching subagents.

```python
class SequenceAlignmentEngine(Protocol):
    def align(
        self,
        reference: ProteinChain,
        target: ProteinChain,
        settings: SequenceAlignmentSettings,
    ) -> SequenceAlignmentResult: ...
```

```python
class StructuralAlignmentEngine(Protocol):
    def align(
        self,
        reference: ProteinChain,
        target: ProteinChain,
        settings: StructuralAlignmentSettings,
    ) -> StructuralAlignmentResult: ...
```

```python
class ResidueMapper(Protocol):
    def build_correspondence(
        self,
        reference: ProteinChain,
        target: ProteinChain,
        alignment: SequenceAlignmentResult | StructuralAlignmentResult,
    ) -> list[ResidueCorrespondence]: ...
```

```python
class StructureSuperposer(Protocol):
    def fit(
        self,
        correspondences: Sequence[ResidueCorrespondence],
        reference: ProteinChain,
        target: ProteinChain,
    ) -> SuperpositionResult: ...
```

```python
class MutationDetector(Protocol):
    def detect(
        self,
        correspondences: Sequence[ResidueCorrespondence],
    ) -> list[MutationEvent]: ...
```

```python
class VisualizationBackend(Protocol):
    def apply(self, state: VisualizationState, analysis: AnalysisResult) -> None: ...
    def reset(self) -> None: ...
```

Subagents must not rename these concepts or introduce competing domain representations without approval from the main agent.

---

## 10. Parallelization Strategy

Use Luna Extra High subagents only after contracts above are frozen.

Recommended independent lanes:

### Lane A — Parsing and sequence mapping

Owns:

- PDB/mmCIF/FASTA parsing;
- normalized residue IDs;
- sequence alignment;
- sequence-based correspondence.

### Lane B — Geometry

Owns:

- Kabsch superposition;
- strict RMSD;
- refinement algorithm;
- Cα displacement;
- backbone/side-chain/all-heavy RMSD;
- symmetry-aware atom matching.

### Lane C — Mutation engine

Owns:

- mutation detection;
- BLOSUM62 values;
- Grantham values;
- physicochemical classifications.

### Lane D — US-align integration

Owns:

- executable discovery;
- invocation;
- output parsing;
- structural correspondence;
- transform/TM-score integration.

### Lane E — Scientific fixtures and independent validation

Owns:

- artificial transformation fixtures;
- known mutation fixtures;
- insertion-code fixtures;
- mmCIF numbering fixtures;
- symmetric-residue tests;
- cross-validation calculations.

### Main agent

Retains:

- architecture;
- interface contracts;
- application services;
- project serialization;
- GUI controller;
- final PyMOL integration;
- release verification.

Do not dispatch two agents to edit the same file concurrently.

---

# Implementation Tasks

## Task 1: Repository, Quality Gates, and Domain Skeleton

**Files:**
- Create `pyproject.toml`
- Create `src/structlens/__init__.py`
- Create domain model files under `src/structlens/core/models/`
- Create `tests/unit/core/models/test_models.py`
- Create `docs/architecture.md`
- Create `docs/implementation-log.md`

**Produces:**
- importable `structlens` package;
- frozen core identifiers and enums;
- lint/type/test configuration.

- [ ] **Step 1: Write model tests first**

Tests must verify:

```python
def test_residue_id_distinguishes_insertion_codes():
    a = ResidueId("x", "1", "A", "100", None, "GLY")
    b = ResidueId("x", "1", "A", "100", "A", "GLY")
    assert a != b
```

```python
def test_correspondence_status_uses_enum():
    correspondence = ResidueCorrespondence(...)
    assert isinstance(correspondence.status, CorrespondenceStatus)
```

- [ ] **Step 2: Run model tests and verify failure**

Run:

```bash
pytest tests/unit/core/models/test_models.py -v
```

Expected: failure because the package/models do not exist.

- [ ] **Step 3: Implement minimal models and enums**

Implement the exact concepts defined in Sections 2, 4, and 5.

- [ ] **Step 4: Configure quality gates**

`pyproject.toml` must include:

- Python `>=3.11`;
- pytest;
- Ruff;
- mypy;
- NumPy;
- Biopython.

- [ ] **Step 5: Verify**

```bash
pytest -q
ruff check .
mypy src/structlens
```

- [ ] **Step 6: Commit**

```bash
git add pyproject.toml src tests docs
git commit -m "feat: establish StructLens domain model"
```

---

## Task 2: PDB, mmCIF, FASTA Parsing and Normalization

**Files:**
- Create parsing modules under `src/structlens/core/parsing/`
- Create `tests/unit/core/parsing/`
- Add fixtures to `tests/fixtures/parsing/`

**Consumes:** `ResidueId`, `ResidueNumbering`, `ProteinStructure`, `ProteinChain`.

**Produces:**

```python
def load_structure(path: Path) -> ProteinStructure
def load_fasta(path: Path) -> list[ProteinSequence]
```

- [ ] Write failing tests for PDB residue numbering, insertion codes, alternate locations, missing atoms, mmCIF author/label numbering, and gzip input.
- [ ] Verify tests fail.
- [ ] Implement PDB normalization.
- [ ] Implement mmCIF normalization.
- [ ] Implement FASTA parsing.
- [ ] Preserve original atom names and coordinates.
- [ ] Do not discard non-standard residues; mark them.
- [ ] Verify tests.
- [ ] Run Ruff and mypy.
- [ ] Commit:

```bash
git commit -am "feat: parse and normalize protein inputs"
```

**Acceptance:** the same biological chain loaded from PDB and equivalent mmCIF fixtures produces consistent canonical residue order while preserving source numbering metadata.

---

## Task 3: Sequence Alignment Engine

**Files:**
- `src/structlens/core/alignment/protocols.py`
- `src/structlens/core/alignment/sequence.py`
- `src/structlens/core/mapping/sequence_mapper.py`
- `src/structlens/core/metrics/sequence_metrics.py`
- `tests/unit/core/alignment/test_sequence_alignment.py`

**Produces:** `SequenceAlignmentResult`.

Default settings:

```python
SequenceAlignmentSettings(
    substitution_matrix="BLOSUM62",
    gap_open=-10.0,
    gap_extend=-0.5,
)
```

- [ ] Write tests for identical sequences.
- [ ] Write tests for substitutions.
- [ ] Write tests for insertion and deletion.
- [ ] Write a numbering-offset test proving that residue number equality is not required.
- [ ] Verify failure.
- [ ] Implement alignment.
- [ ] Convert aligned positions into explicit `ResidueCorrespondence` objects.
- [ ] Compute identity, similarity, and coverage from the mapped alignment.
- [ ] Verify tests and commit:

```bash
git commit -am "feat: add sequence-guided residue mapping"
```

---

## Task 4: Kabsch Superposition and Strict Structural Metrics

**Files:**
- `src/structlens/core/alignment/superposition.py`
- `src/structlens/core/geometry/kabsch.py`
- `src/structlens/core/geometry/rmsd.py`
- `src/structlens/core/geometry/displacement.py`
- `tests/unit/core/geometry/`

**Produces:** `SuperpositionResult`.

```python
@dataclass(frozen=True, slots=True)
class SuperpositionResult:
    rotation: np.ndarray
    translation: np.ndarray
    strict_rmsd_angstrom: float
    atom_count: int
    residue_count: int
```

- [ ] Create a fixture with known rotation + translation and no internal deformation.
- [ ] Test that fitted RMSD is approximately zero.
- [ ] Test that translation/rotation do not alter pairwise internal distances.
- [ ] Test Cα displacement after superposition.
- [ ] Implement Kabsch with explicit reflection handling.
- [ ] Verify numerical tolerance with `np.testing.assert_allclose`.
- [ ] Commit:

```bash
git commit -am "feat: add reproducible structural superposition"
```

---

## Task 5: Refinement Without Silent Data Loss

**Files:**
- `src/structlens/core/alignment/refinement.py`
- `tests/unit/core/alignment/test_refinement.py`

**Produces:** `RefinementResult`.

```python
@dataclass(frozen=True, slots=True)
class RefinementResult:
    refined_superposition: SuperpositionResult
    included_alignment_indices: tuple[int, ...]
    excluded_alignment_indices: tuple[int, ...]
    cutoff_angstrom: float
    cycles: int
```

- [ ] Test a structure with one artificial outlier.
- [ ] Verify strict RMSD remains unchanged in stored analysis.
- [ ] Verify refined result excludes the expected alignment index.
- [ ] Verify excluded correspondences remain present and are marked `is_outlier=True`.
- [ ] Implement deterministic iterative refinement.
- [ ] Commit:

```bash
git commit -am "feat: add transparent RMSD refinement"
```

---

## Task 6: Symmetry-Aware Residue RMSD

**Files:**
- `src/structlens/core/geometry/symmetry.py`
- Extend `rmsd.py`
- `tests/unit/core/geometry/test_residue_rmsd.py`

**Initial symmetric equivalences must cover at least:**

- ASP OD1/OD2
- GLU OE1/OE2
- ARG NH1/NH2
- PHE CD1/CD2 and CE1/CE2
- TYR CD1/CD2 and CE1/CE2

- [ ] Create fixtures where only symmetric atom naming is swapped.
- [ ] Prove naive name matching yields an inflated value in the fixture.
- [ ] Implement valid permutations.
- [ ] Return the minimum chemically valid RMSD.
- [ ] Implement backbone, side-chain, and all-heavy-atom metrics.
- [ ] Define behavior for missing atoms: metric is `None` when the minimum required atom set cannot be formed; never silently substitute atoms.
- [ ] Commit:

```bash
git commit -am "feat: calculate symmetry-aware residue RMSD"
```

---

## Task 7: Mutation Detection and Descriptors

**Files:**
- mutation modules under `src/structlens/core/mutations/`
- `tests/unit/core/mutations/`

- [ ] Test conserved residue.
- [ ] Test substitution `S -> T`.
- [ ] Test insertion.
- [ ] Test deletion.
- [ ] Test non-standard residue.
- [ ] Implement mutation detection.
- [ ] Add an embedded, test-verified BLOSUM62 lookup.
- [ ] Add a documented Grantham-distance lookup for canonical amino-acid substitutions.
- [ ] Add explicit physicochemical categories.
- [ ] Ensure none of the output copy implies functional impact.
- [ ] Commit:

```bash
git commit -am "feat: detect and characterize residue mutations"
```

---

## Task 8: US-align Adapter

**Files:**
- `src/structlens/integrations/usalign/executable.py`
- `src/structlens/integrations/usalign/adapter.py`
- `src/structlens/integrations/usalign/parser.py`
- `tests/unit/integrations/usalign/`
- `tests/integration/test_usalign.py`

**Produces:** implementation of `StructuralAlignmentEngine`.

- [ ] Write parser tests against committed, small text fixtures representing supported US-align output.
- [ ] Test missing executable behavior with a clear actionable error.
- [ ] Test executable discovery from explicit configuration and PATH.
- [ ] Implement `subprocess.run` with argument lists, never shell interpolation.
- [ ] Parse TM-score, aligned residue pairs, transform when available, and metadata.
- [ ] Convert structural pairs to StructLens `ResidueCorrespondence`.
- [ ] Record US-align executable version in analysis provenance.
- [ ] Do not download binaries automatically during analysis.
- [ ] Commit:

```bash
git commit -am "feat: integrate US-align structural mapping"
```

---

## Task 9: AUTO Alignment Policy and Mapping Validation

**Files:**
- `src/structlens/core/mapping/validator.py`
- `src/structlens/application/analysis_service.py`
- `tests/unit/application/test_auto_alignment.py`

- [ ] Test AUTO chooses sequence branch for identity >= 0.30 and coverage >= 0.70.
- [ ] Test AUTO chooses structural branch when either default threshold is missed.
- [ ] Test the decision string records identity, coverage, thresholds, and chosen engine.
- [ ] Test locked manual mappings survive validation unchanged.
- [ ] Implement geometric validation without deleting mapped positions.
- [ ] Commit:

```bash
git commit -am "feat: orchestrate automatic alignment strategy"
```

---

## Task 10: Key Residues and Residue Groups

**Files:**
- extend domain models
- `src/structlens/application/project_state.py`
- `tests/unit/application/test_key_residues.py`

Key-residue groups support:

- Catalytic
- Binding site
- Mutation hotspot
- Experimental mutations
- Custom user-defined names

- [ ] Test adding a key residue by reference locator.
- [ ] Test mapping it to target equivalent through correspondence.
- [ ] Test a key residue that is deleted in the target.
- [ ] Test serialization/deserialization.
- [ ] Commit:

```bash
git commit -am "feat: add biologically meaningful residue groups"
```

---

## Task 11: Analysis Result and Project Serialization

**Files:**
- `src/structlens/core/models/results.py`
- `src/structlens/application/project_state.py`
- `tests/unit/application/test_project_serialization.py`

Project file extension:

```text
.structlens.json
```

Persist:

- StructLens version;
- input structure identifiers;
- SHA-256 hashes for source files when available;
- chain selections;
- alignment mode;
- alignment settings;
- mapping source;
- full residue correspondence;
- locked mappings;
- strict/refined results;
- excluded outliers;
- mutations;
- key residues;
- US-align version when used;
- visualization state separately from scientific state.

- [ ] Round-trip test project -> JSON -> project.
- [ ] Verify all scientific results survive round-trip.
- [ ] Verify visualization changes do not mutate scientific hashes/results.
- [ ] Add explicit schema version.
- [ ] Commit:

```bash
git commit -am "feat: persist reproducible StructLens projects"
```

---

## Task 12: CLI

**Files:**
- `src/structlens/cli/main.py`
- `src/structlens/cli/formatting.py`
- `tests/integration/test_cli.py`

Initial commands:

```bash
structlens compare reference.pdb target.pdb
structlens compare reference.pdb target.pdb --mode auto
structlens compare reference.pdb target.pdb --mode sequence
structlens compare reference.pdb target.pdb --mode structure
structlens compare reference.pdb target.pdb --key-residues A:70,A:73,A:130,A:166
structlens compare reference.pdb target.pdb --export results.csv
```

Output must label strict and refined values separately.

- [ ] Write CLI tests first.
- [ ] Implement commands.
- [ ] Ensure non-zero exit codes for invalid input or missing required structural backend.
- [ ] Commit:

```bash
git commit -am "feat: expose StructLens analysis through CLI"
```

---

## Task 13: Publication-Ready Export Service

**Files:**
- `src/structlens/application/export_service.py`
- `src/structlens/application/export_models.py`
- `src/structlens/application/image_export_service.py`
- `tests/unit/application/test_export_service.py`
- `tests/unit/application/test_image_export_service.py`
- `tests/regression/exports/`

StructLens must allow the user to export **all analysis products that can reasonably be represented as a table or image**.

### Required export formats

Tabular/scientific data:

- XLSX
- CSV
- TSV
- JSON

Publication raster graphics:

- JPEG
- TIFF

### XLSX workbook requirements

Use `openpyxl`.

For a complete analysis export, generate one workbook with clearly named English sheets as applicable:

```text
Summary
Structures
Alignment
Residue Mapping
Mutations
Key Residues
Outliers
Per-Residue Metrics
Local Analysis
Settings
Provenance
```

For one-reference-vs-multiple-target analyses, additionally provide:

```text
Target Summary
```

Workbook requirements:

- freeze header rows;
- apply readable column widths;
- enable filters;
- use numeric cells for scientific values rather than preformatted strings;
- keep units in headers, e.g. `Cα displacement (Å)`;
- use English sheet names and headers;
- represent unavailable values as blank cells, not `NaN`, `None`, or textual placeholders;
- include StructLens version and analysis provenance;
- include source structure identifiers and hashes where available;
- include chosen alignment method and AUTO decision rationale;
- include strict/refined settings and excluded residue counts;
- preserve full residue identifiers, including chain and insertion code;
- do not use decorative spreadsheet styling that interferes with scientific readability.

### Image export requirements

Any view or figure generated by StructLens that is meaningful for a manuscript, presentation, or supplementary material must be exportable when technically possible.

At minimum support high-quality export of:

- current PyMOL comparison view configured by StructLens;
- mutation-focused view;
- key-residue view;
- outlier view;
- structural-displacement colored view;
- active-site view;
- legends;
- residue/mutation summary figure when such a figure is shown in the GUI;
- multiple-target summary visualization when such a figure is shown in the GUI.

Required formats:

```text
JPEG
TIFF
```

Required quality controls:

- resolution selector with at least 300 dpi and 600 dpi publication presets;
- width/height controls in pixels or physical dimensions;
- antialiasing enabled where supported;
- configurable background: white, transparent where technically compatible, or current PyMOL background;
- preserve visible legends when requested;
- export without GUI chrome unless the user explicitly chooses a screenshot-style export;
- do not upscale a low-resolution screenshot and label it 600 dpi;
- render at the requested pixel dimensions from the source scene whenever the rendering backend supports it.

For TIFF:

- prefer lossless or publication-appropriate compression;
- preserve full requested pixel dimensions;
- write DPI metadata when supported.

For JPEG:

- expose a high-quality default;
- avoid excessive compression artifacts;
- use JPEG only for opaque images.

If transparent background is requested, disable JPEG and explain why; allow TIFF where the backend supports alpha or export with an explicit solid background.

### Export naming

Default export names should be descriptive, for example:

```text
StructLens_CTXM15_vs_SHV1_residue_mapping.xlsx
StructLens_CTXM15_vs_SHV1_mutations.xlsx
StructLens_CTXM15_vs_SHV1_structural_deviation_600dpi.tiff
StructLens_CTXM15_vs_SHV1_active_site_300dpi.jpg
```

### Export UI

The export workflow must show:

- what will be exported;
- destination;
- format;
- units/resolution for images;
- whether legends/labels are included;
- concise explanation of recommended settings.

Recommended helper copy:

```text
300 dpi
Suitable for most journal figures at final publication size.

600 dpi
Recommended for line-rich molecular graphics, labels, and high-resolution
publication workflows.
```

### TDD steps

- [ ] Write failing XLSX tests before implementation.
- [ ] Snapshot-test workbook sheet names and column order.
- [ ] Verify scientific metric cells are numeric.
- [ ] Verify unavailable metrics produce blank cells.
- [ ] Verify insertion codes and target numbering survive export.
- [ ] Verify one-reference-vs-multiple-target workbook layout.
- [ ] Write image-export tests for requested dimensions, file format, and DPI metadata where supported.
- [ ] Verify JPEG export rejects transparency with an actionable English message.
- [ ] Verify TIFF/JPEG files are readable by Pillow after export.
- [ ] Verify an image requested at 600 dpi is rendered at the expected pixel dimensions for the requested publication size.
- [ ] Verify image export does not include plugin GUI chrome by default.
- [ ] Implement XLSX, CSV, TSV, JSON, JPEG, and TIFF export.
- [ ] Add regression fixtures for workbook metadata and publication image dimensions.
- [ ] Commit:

```bash
git commit -am "feat: export publication-ready data and figures"
```

**Acceptance:** every scientific table exposed by StructLens has an XLSX export path, and every StructLens-generated molecular/figure view exposed to the user has a JPEG/TIFF export path or an explicitly documented technical reason why that specific view cannot be raster-exported.

---

## Task 14: PyMOL Adapter and Safe Namespace

**Files:**
- `src/structlens/integrations/pymol/`
- `tests/unit/integrations/pymol/`
- integration tests runnable inside a PyMOL-capable environment

StructLens-owned names must use a predictable prefix, for example:

```text
SL_<project>_<target>_mutations
SL_<project>_<target>_key_residues
SL_<project>_<target>_outliers
SL_<project>_<target>_alignment
```

- [ ] Unit-test selection-string generation without importing a live PyMOL process where possible.
- [ ] Implement conversion from PyMOL objects/chains to normalized domain structures.
- [ ] Implement transformation application.
- [ ] Implement namespaced selections.
- [ ] Snapshot any PyMOL visualization state that StructLens must alter.
- [ ] `reset()` restores only StructLens-owned changes and must not delete unrelated user objects/selections.
- [ ] Commit:

```bash
git commit -am "feat: add safe PyMOL integration layer"
```

---


## Task 14A: Application Icon and Visual Identity Integration

**Files:**
- `src/structlens/plugin/assets/structlens_icon.png`
- packaging metadata as required by the target platform
- `README.md`
- plugin entrypoint/menu registration files

Use the approved StructLens icon supplied with the project.

- [ ] Copy the approved master icon into `src/structlens/plugin/assets/structlens_icon.png`.
- [ ] Verify the icon remains legible at 16, 24, 32, 48, 64, 128, and 256 px.
- [ ] Use the icon in the StructLens panel/window where Qt/PyMOL supports custom icons.
- [ ] Use it for plugin/menu entry where supported.
- [ ] Generate `.ico`/`.icns` only if a packaging target requires them.
- [ ] Do not redesign the icon during implementation.
- [ ] Add one documentation screenshot or README placement using the approved identity.
- [ ] Commit:

```bash
git commit -am "feat: integrate StructLens application icon"
```

---

## Task 15: GUI Architecture With Impeccable Gate

**Files:**
- plugin GUI modules under `src/structlens/plugin/gui/`
- `src/structlens/plugin/controller.py`
- `src/structlens/plugin/workers.py`
- GUI tests where practical
- `docs/implementation-log.md`

Before code:

- [ ] Discover `impeccable`.
- [ ] If available, invoke it and record the invocation.
- [ ] If unavailable, record exactly that it was unavailable.
- [ ] Produce a one-page GUI architecture note in `docs/architecture.md` mapping each UI page to application services and state.

Implementation:

- [ ] Build the six-section panel specified in Section 7.
- [ ] Add English explanatory subtitle/helper text to every major page.
- [ ] Add scientifically meaningful tooltips/help for every non-obvious workflow option and metric listed in Section 7.
- [ ] Ensure alignment mode cards/options explain “what it does” and “when to use it” before the user runs analysis.
- [ ] Add contextual empty states that tell the user the next valid action.
- [ ] Keep UI state separate from domain and application state.
- [ ] Use worker/QThread-compatible execution for long operations.
- [ ] Support progress, cancellation, and actionable errors.
- [ ] Disable actions whose prerequisites are not satisfied.
- [ ] Ensure running an analysis does not freeze the main PyMOL UI.
- [ ] Commit:

```bash
git commit -am "feat: add StructLens PyMOL scientific interface"
```

---

## Task 16: Visualization Renderer

**Files:**
- `src/structlens/application/visualization_service.py`
- `src/structlens/plugin/visualization/renderer.py`
- `src/structlens/plugin/visualization/presets.py`
- `src/structlens/plugin/visualization/legends.py`
- `tests/unit/plugin/visualization/`

- [ ] Implement independent filter/data/color/representation state.
- [ ] Add mutation-only highlighting.
- [ ] Add key-residue highlighting.
- [ ] Add mutated-key-residue highlighting.
- [ ] Add outlier highlighting.
- [ ] Add Cα displacement scale with visible Å legend.
- [ ] Add backbone RMSD scale with visible Å legend.
- [ ] Add representation selector.
- [ ] Implement six named presets.
- [ ] Verify preset application does not change analysis data.
- [ ] Commit:

```bash
git commit -am "feat: visualize mutations and structural differences"
```

---

## Task 17: Bidirectional Residue Browser Interaction

**Files:**
- `src/structlens/plugin/gui/residues_page.py`
- `src/structlens/plugin/controller.py`
- `src/structlens/integrations/pymol/selections.py`
- tests for controller events

- [ ] Clicking a residue table row selects/focuses corresponding residues in PyMOL.
- [ ] Mutations display canonical reference-based notation plus explicit target numbering.
- [ ] Selection changes generated by StructLens update the residue inspector.
- [ ] Avoid aggressive timers/polling.
- [ ] User-generated unrelated PyMOL selections do not corrupt the plugin state.
- [ ] Commit:

```bash
git commit -am "feat: synchronize residue inspection with PyMOL"
```

---

## Task 18: Mutation Inspector and Local Environment

**Files:**
- `src/structlens/core/metrics/local_metrics.py`
- mutation/residue GUI pages
- `tests/unit/core/metrics/test_local_metrics.py`

Version 0.1 local-environment definition:

- center = reference residue atoms;
- default radius = 5.0 Å;
- select neighboring reference residues that have at least one heavy atom within the radius;
- transfer that neighbor set through the authoritative correspondence map;
- compute local backbone RMSD only for mapped neighbors with complete required backbone atoms.

- [ ] Test neighborhood definition on a synthetic geometry fixture.
- [ ] Test deletion in the mapped local environment.
- [ ] Display local RMSD and local mapped residue count.
- [ ] Add mutation inspector action `Focus in PyMOL`.
- [ ] Commit:

```bash
git commit -am "feat: inspect local structural effects of mutations"
```

---

## Task 19: Multiple Targets Summary

**Files:**
- `src/structlens/application/analysis_service.py`
- results GUI
- export service
- tests

For one reference versus N targets, report:

- identity;
- coverage;
- TM-score when available;
- strict Cα RMSD;
- refined Cα RMSD;
- mutation count;
- key-site mutation count.

- [ ] Test three synthetic targets.
- [ ] Ensure each target has an independent correspondence set and transform.
- [ ] Ensure changing the selected target in the GUI changes the displayed residue table without recomputing unchanged analyses.
- [ ] Commit:

```bash
git commit -am "feat: compare one reference against multiple structures"
```

---

## Task 20: Scientific Regression and Golden Datasets

**Files:**
- `tests/scientific/`
- `tests/regression/`
- `docs/validation.md`

Required validation cases:

1. structure against itself: RMSD approximately 0, zero mutations;
2. rigid rotation/translation: fitted RMSD approximately 0;
3. single known substitution;
4. known insertion/deletion;
5. residue numbering offset;
6. insertion codes such as 100/100A/100B;
7. PDB vs equivalent mmCIF numbering;
8. swapped symmetric side-chain atom names;
9. missing backbone atom;
10. one extreme structural outlier;
11. manual locked mapping;
12. sequence-vs-US-align branch selection.

- [ ] Implement golden expected outputs in small human-readable fixtures.
- [ ] Verify results independently with direct NumPy calculations where feasible.
- [ ] Where PyMOL is available, cross-check equivalent atom selections under identical conditions.
- [ ] Where US-align is available, verify parsed TM-score and aligned pairs against raw output fixture.
- [ ] Document tolerances and rationale in `docs/validation.md`.
- [ ] Commit:

```bash
git commit -am "test: validate StructLens scientific calculations"
```

---

## Task 21: GUI Critique and Polish

Before this task:

- [ ] Invoke `impeccable` again if available.
- [ ] Record critique/polish invocation in `docs/implementation-log.md`.

Review:

- information hierarchy;
- excessive panels/cards;
- table density;
- error and empty states;
- unit visibility;
- legends;
- accessibility and non-color cues;
- long labels and small PyMOL window sizes;
- dark/light host behavior where available;
- focus order;
- disabled-state clarity;
- progress/cancellation feedback.

No architectural redesign is permitted here unless a usability defect demonstrates that the architecture blocks a required workflow.

Commit:

```bash
git commit -am "refactor: polish StructLens scientific GUI"
```

---

## Task 22: Documentation, Packaging, and PyMOL Installation

**Files:**
- `README.md`
- `docs/scientific-methods.md`
- packaging metadata
- plugin installation instructions
- `CHANGELOG.md`

README must contain:

- what StructLens does;
- supported formats;
- installation;
- US-align dependency;
- first analysis tutorial;
- mutation workflow;
- key-residue workflow;
- strict vs refined RMSD explanation;
- Cα displacement vs residue RMSD explanation;
- limitations.
- XLSX export workflow;
- JPEG/TIFF publication export workflow;
- 300/600 dpi guidance;
- explanation of each GUI workflow option;
- statement that StructLens v0.1 is English-only;
- approved StructLens icon/visual identity usage.

`scientific-methods.md` must document:

- sequence alignment settings;
- AUTO policy;
- structural alignment backend;
- Kabsch implementation;
- refinement;
- residue metrics;
- symmetry handling;
- mutation descriptors;
- local RMSD definition.

- [ ] Test installation in a clean environment.
- [ ] Test CLI.
- [ ] Test plugin loading in the supported PyMOL environment.
- [ ] Commit:

```bash
git commit -am "docs: document and package StructLens"
```

---

# 11. Deferred Features — Not Part of v0.1 Acceptance

Do not implement these during v0.1 unless all required tasks are complete and a separate approved plan is created:

- automatic online PDB/UniProt/AlphaFold download;
- Foldseek backend;
- hydrogen-bond network comparison;
- salt-bridge comparison;
- ligand-contact network comparison;
- energetic prediction;
- pathogenicity prediction;
- machine-learning mutation effects;
- MD trajectory comparison;
- FEP integration;
- DockLens direct integration;
- PDF/HTML report generation;
- N-by-N all-pairs structural matrix;
- cloud services.

This is deliberate YAGNI scope control.

---

# 12. PyMOL Visualization Safety Rules

1. Never use generic selection names such as `sele`, `site`, or `mut`.
2. Never call broad destructive commands on objects not owned by StructLens.
3. Preserve original source objects.
4. Transform display copies or use a reversible strategy when necessary.
5. Reset must restore only state changed by StructLens.
6. Every generated selection/object must be traceable to project + target + purpose.
7. The renderer must not recompute scientific metrics.
8. Visualization presets must not modify correspondence or analysis settings.

---

# 13. Error Model

Use typed application errors. At minimum:

```python
class StructLensError(Exception): ...
class InputFormatError(StructLensError): ...
class ChainNotFoundError(StructLensError): ...
class MappingError(StructLensError): ...
class InsufficientAtomsError(StructLensError): ...
class USAlignNotFoundError(StructLensError): ...
class USAlignExecutionError(StructLensError): ...
class ProjectSchemaError(StructLensError): ...
```

GUI presents concise user-facing messages and an expandable technical detail where appropriate.

CLI sends concise errors to stderr and returns a non-zero exit code.

Never convert a scientific failure into zero, empty RMSD, or a fabricated result.

---

# 14. Performance Rules

v0.1 must comfortably support one reference versus several dozen ordinary single-chain proteins without UI freezing.

Implementation rules:

- cache parsed structures by stable input hash;
- cache sequence alignments by sequence/settings hash;
- cache structural mappings by structure/settings hash;
- calculate side-chain RMSD lazily when large target sets make eager computation expensive;
- long US-align/batch analysis runs in worker threads/process-safe workers appropriate to the PyMOL Qt environment;
- progress must identify current target;
- cancel must stop queued work and terminate owned external US-align process when safe;
- do not introduce complex multiprocessing until profiling demonstrates need.

---

# 15. Scientific Acceptance Criteria

A release candidate is scientifically acceptable only if all conditions below pass.

### Mapping

- equivalent residues are not inferred from residue number equality;
- insertions/deletions are explicit;
- manual locked mapping cannot be silently overwritten;
- PDB insertion codes work;
- mmCIF author/label numbering is retained.

### Geometry

- rigid-body invariance passes;
- strict RMSD is reproducible;
- refinement reports excluded residues;
- Cα displacement is not mislabeled as RMSD;
- missing atoms do not silently substitute unrelated atoms;
- symmetric side chains are handled correctly.

### Mutations

- conserved/substitution/insertion/deletion/nonstandard are distinct;
- notation is reference-oriented and target numbering remains visible;
- BLOSUM/Grantham descriptors are reproducible;
- no unsupported functional conclusions appear.

### Reproducibility

- project JSON round-trips;
- source file hashes are recorded when possible;
- StructLens version is recorded;
- US-align version is recorded when used;
- settings and exclusions are recorded.

### Exports

- every scientific table visible in the product can be exported to XLSX;
- complete-analysis XLSX workbooks contain English sheet names, units, settings, and provenance;
- publication images can be exported as JPEG or TIFF;
- 300 dpi and 600 dpi presets produce internally consistent requested dimensions/metadata;
- raster export is rendered from the molecular scene at requested dimensions, not created by upscaling a low-resolution screenshot;
- export filenames are descriptive and deterministic enough for scientific workflows;
- export errors never silently create partial/corrupt files.

### Usability and language

- the complete v0.1 user-facing product is in English;
- every workflow choice explains what it does and when it should be used;
- every advanced scientific metric has contextual help;
- routine analysis can be completed without reading external documentation;
- empty states guide the user toward the next valid action;
- accessibility does not rely on color alone.

### PyMOL

- plugin loads;
- analysis runs without blocking the UI;
- row-to-structure focus works;
- visualization reset does not destroy unrelated user state;
- mutations/key residues/outliers can be highlighted independently.

---

# 16. Definition of Done per Feature

A feature is not done because a function returns a value or because one test passes.

Each feature must pass:

```text
unit test
    ↓
integration test when applicable
    ↓
regression/golden test when scientific
    ↓
scientific sanity check
    ↓
code review
```

Numerical features must document:

- input atoms;
- residue correspondences;
- transformation state;
- units;
- exclusions;
- tolerances.

---

# 17. Required Verification Commands

Before completion, run the commands appropriate to the environment and record results.

```bash
pytest -q
ruff check .
mypy src/structlens
```

If project packaging uses build:

```bash
python -m build
```

CLI smoke test:

```bash
structlens --help
```

Scientific self-comparison smoke test:

```bash
structlens compare tests/fixtures/proteins/reference.pdb \
  tests/fixtures/proteins/reference.pdb \
  --mode sequence
```

Expected scientific outcome:

- zero mutations;
- full or expected fixture coverage;
- RMSD within documented floating-point tolerance of zero.

When PyMOL is available, run the plugin smoke workflow and record it in `docs/implementation-log.md`.

When US-align is available, run one real adapter integration test and record detected version.

---

# 18. Completion Review Checklist

The main Luna agent must perform this review before calling `verification-before-completion`.

## Spec coverage

- [ ] PDB/mmCIF/FASTA
- [ ] PyMOL objects
- [ ] one reference vs multiple targets
- [ ] chain selection
- [ ] Auto/Sequence/Structure/Manual
- [ ] explicit residue mapping
- [ ] locked mapping
- [ ] mutation detection
- [ ] insertion/deletion
- [ ] strict RMSD
- [ ] refined RMSD
- [ ] Cα displacement
- [ ] backbone RMSD
- [ ] side-chain RMSD
- [ ] symmetry-aware matching
- [ ] key residues
- [ ] mutation highlights
- [ ] key-residue highlights
- [ ] outlier highlights
- [ ] displacement/RMSD color modes
- [ ] visualization presets
- [ ] residue browser
- [ ] PyMOL focus
- [ ] CSV/TSV/JSON export
- [ ] XLSX export for every scientific table
- [ ] JPEG export for StructLens-generated views
- [ ] TIFF export for StructLens-generated views
- [ ] 300/600 dpi publication presets
- [ ] approved StructLens icon packaged
- [ ] English-only UI and exports
- [ ] workflow explanations/tooltips for all non-obvious options
- [ ] `.structlens.json`
- [ ] CLI
- [ ] scientific validation
- [ ] GUI design pass
- [ ] GUI polish pass
- [ ] documentation

## Placeholder scan

Search for:

```bash
grep -RInE 'TODO|TBD|FIXME|pass[[:space:]]*$|NotImplementedError' \
  src tests docs
```

Every match must be justified or removed before release. Deferred features belong in documentation, not placeholder code.

## Architecture scan

Verify:

```bash
grep -RIn 'from pymol\|import pymol' src/structlens/core src/structlens/application
```

Expected: no matches.

## Scientific naming scan

Search GUI/docs/source for single-residue values mislabeled as RMSD. `Cα displacement` terminology must be consistent.

---

# 19. Milestones and Review Gates

### Milestone A — Scientific core

Tasks 1–7.

Gate:
- domain model frozen;
- sequence mapping validated;
- geometry validated;
- mutations validated.

Request code review.

### Milestone B — Structural backend and orchestration

Tasks 8–13.

Gate:
- US-align adapter works;
- AUTO strategy is transparent;
- project serialization and CLI work.

Request code review.

### Milestone C — PyMOL product

Tasks 14–19.

Gate:
- plugin works without corrupting user state;
- GUI remains responsive;
- mutation/residue visualization is usable.

Request code review.

### Milestone D — Validation and release

Tasks 20–22.

Gate:
- scientific suite passes;
- Impeccable second pass completed if available;
- documentation complete;
- packaging verified.

Use `superpowers:verification-before-completion`.

---

# 20. Final Product Behavior

The intended end-user flow is:

```text
Load structures
      ↓
Choose reference
      ↓
Choose one or more targets
      ↓
Choose Auto / Sequence / Structure / Manual
      ↓
Optionally mark key residues
      ↓
Analyze
      ↓
StructLens builds explicit residue correspondence
      ↓
Superposes structures
      ↓
Calculates strict + optional refined global metrics
      ↓
Detects substitutions / insertions / deletions
      ↓
Calculates residue-level structural metrics
      ↓
Displays synchronized residue/mutation table
      ↓
Highlights mutations, key residues, outliers, or structural displacement in PyMOL
      ↓
Exports a reproducible project, XLSX result workbooks, and publication-quality JPEG/TIFF figures
```

A user should be able to answer, with traceable evidence:

> Which residue in the target corresponds to my reference residue?

> Is it conserved or mutated?

> If mutated, what substitution occurred?

> How different is the local structure?

> Was this residue excluded from a refined RMSD?

> Can I see the difference immediately in PyMOL?

If StructLens cannot answer one of these reliably, it must say that the value is unavailable rather than fabricate or infer unsupported data.

---

# 21. Execution Handoff for Luna Extra High

Recommended execution mode:

**Subagent-Driven Development.**

The main Luna Extra High agent should:

1. create an isolated worktree;
2. complete Task 1 and freeze domain contracts;
3. dispatch independent Luna Extra High subagents for lanes A–E;
4. review each subagent result before integration;
5. integrate through application services;
6. build PyMOL adapter;
7. apply the Impeccable gate before GUI coding;
8. build GUI and visualization;
9. run scientific regression;
10. apply the Impeccable polish pass;
11. run full verification;
12. request final code review;
13. finish the development branch only after all gates pass.

Do not ask the user to manually resolve implementation details already defined in this plan. Escalate only genuine scientific ambiguity that can change the meaning of the analysis.

---

# 22. First Instruction to Give Luna

Use the following as the opening execution directive:

```text
Implement StructLens according to StructLens.md.

Before touching code, read the entire plan. Use Superpowers. Start with
using-git-worktrees and test-driven-development. Use subagent-driven-development
after the domain interfaces are frozen.

The scientific correspondence map is authoritative; PyMOL is a visualization
backend only. Never calculate residue equivalence from residue numbers alone.
Never silently remove residues from RMSD reporting. Never label a single Cα
distance as RMSD.

Use Luna Extra High for the main agent and independent subagents. Before GUI
implementation, discover and invoke the Impeccable skill if it exists in this
environment; use it again for a GUI critique/polish pass. If it is unavailable,
record that fact and follow the GUI specification in StructLens.md.

The user-facing v0.1 product is English-only. Every workflow option must explain
what it does, when it should be used, and any important limitation. Treat
contextual help as a required part of the GUI, not optional documentation.

Package the approved StructLens icon supplied with the project. Every scientific
table must be exportable as XLSX, and every StructLens-generated molecular or
figure view must have a publication-quality JPEG/TIFF export path when technically
possible, including 300 dpi and 600 dpi presets.

Work task-by-task with failing tests first, verification, review, and frequent
commits. Do not implement deferred v0.1 features without a separate approved
plan.
```
