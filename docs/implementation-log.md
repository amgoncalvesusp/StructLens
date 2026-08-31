# Implementation log

- 2026-08-22: Initialized the local git repository because the supplied workspace had no repository metadata.
- 2026-08-22: Applied Superpowers TDD and subagent-driven development for the scientific lanes.
- 2026-08-22: Impeccable skill discovered and invoked for GUI architecture/design; see `PRODUCT.md`, `DESIGN.md`, and `docs/architecture.md`.
- 2026-08-22: Impeccable post-build detector/critique pass ran against the GUI model and visualization renderer; no findings were returned.
- 2026-08-23: Replaced the placeholder Qt tab shell with the Evidence Bench workflow: source/chain loading, asynchronous comparison with cancel state, manual-pair input, evidence tables, project persistence, exports, visualization presets/legends, and PyMOL command-proxy rendering.
- 2026-08-23: Added offscreen Qt smoke coverage and a fake-command PyMOL adapter test; the optional GUI remains lazy-imported outside host environments.
- 2026-08-23: Implemented v0.2 multi-structure contracts (`ComparisonMode`, reference-vs-many, all-vs-all matrices, multiple-structure positions), project round-tripping, and shared selection primitives.
- 2026-08-23: Added bundled-backend resolution with custom → bundled → diagnostic PATH precedence, compressed-input normalization, provenance fields, and upstream US-align notices.
- 2026-08-23: Added renderer-independent scientific datasets for six chart families, real XLSX/JPEG/TIFF publication exports, and a validated atomic `.structlens-pymol` writer with ZIP security checks.
- 2026-08-23: Reorganized the visible desktop workflow around Project, Sequences, Structures, Residues, Charts, and Results while keeping legacy host section identifiers compatible; PyMOL export is explicit and actionable.
- 2026-08-23: Added the independent `structlens-pymol` package in `amgoncalvesusp/pymol-plugins`, including safe bundle reader, namespaced visualization/controller, command API, deterministic release ZIP, tests, and tag-scoped release workflow.
- 2026-08-30: Began StructLens v0.4 structural-reliability work in the isolated
  `codex/structlens-structural-reliability` branch at origin commit `a9d8bf1`.
  Applied the Superpowers brainstorming/writing-plans workflow and recorded the
  design plus resumable TDD plan under `docs/superpowers/` before production edits.
- 2026-08-30 baseline: Python 3.11.9; `python -m pytest -q` reported 110 passed and
  1 skipped; branch coverage was approximately 53%; `python -m ruff check .`
  passed; mypy was not installed in the global interpreter. The sole skipped module
  hard-imports PySide6 although PyQt5 is installed and the release workflow uses
  PyQt5. No production code had been changed at this checkpoint.
- 2026-08-30 audit: confirmed disconnected v0.3 GUI orchestration, incorrect
  C-alpha-only site “backbone” RMSD, unavailable-to-zero envelope semantics,
  unsafe difference-map missingness, false MSE self-mutation, absent resource caps,
  no true cavity detector, incomplete package validation, and several companion
  PyMOL v0.3 rendering/schema defects. These are ordered before new v0.4 claims in
  the implementation plan.
- 2026-08-30 Task 0 Step 3 (distribution-license strategy): wrote five failing
  metadata/release-policy tests first, confirmed each failed for its intended
  reason, then closed them. Two real compliance defects were found and fixed.
  (a) `release.yml` froze PyQt5 (GPL-3) into both official installers while
  passing `--exclude-module PySide6`, so every distributed GUI binary carried a
  GPL Qt binding. The `gui` extra now resolves to PySide6 (LGPL-3), PyQt5 moved
  to an opt-in `gui-pyqt5` extra, and both build jobs freeze PySide6 with
  `--exclude-module PyQt5`. (b) `licenses/MUSCLE-LICENSE.txt` was a two-line stub
  asserting "the full upstream license text is copied into release bundles by
  CI", but `scripts/fetch_backends.py` downloaded only the GPL-3 binary and never
  its licence or corresponding source. The notice now carries the verbatim FSF
  GPL-3 text (707 lines) plus the corresponding-source route; the fetcher
  retrieves `LICENSE` and `muscle-5.3-source.tar.gz` from the same pinned v5.3
  tag; `scripts/verify_release_backends.py` gained `check_muscle_compliance` and
  exits 1 when either file is absent (verified against a synthetic root); and
  `pyproject.toml` package-data ships the source archive. Added
  `docs/distribution-licensing.md` covering the component table, the LGPL
  relinking route, and the four GPL-3 obligations.
- 2026-08-30 Task 0 Step 4 (declared dependencies): SciPy `>=1.10` is declared as
  a runtime dependency and Hypothesis `>=6.90` as a dev dependency, both asserted
  by `tests/unit/test_packaging_metadata.py`. No GPL pocket library was added.
- 2026-08-30 Task 0 gate output: `python -m pytest -q` reported 136 passed with
  no skips; `python -m ruff check .` passed; `python -m mypy src` reported
  success across 109 source files; `git diff --check` was clean. Python 3.11.9.
- 2026-08-30 Task 1 (correct v0.3 scientific semantics): five RED tests were added
  to `tests/unit/application/test_v03_services.py` and each was confirmed to fail
  for its intended reason before any production edit. Four defects fixed in
  `site_service.py`. (a) Both fields named `*_backbone_rmsd_angstrom` measured the
  C-alpha alone, so a target whose carbonyl oxygens were displaced 2.0 A reported a
  perfect 0.0; the new `_backbone_coords` helper collects matched N/CA/C/O in a
  fixed order and returns None for an incomplete backbone instead of pairing
  mismatched atoms. (b) Site-fitted RMSD subtracted centroids only, which removes
  translation but leaves rotation; it now uses `core/geometry/kabsch.py` via
  `_site_fitted_rmsd`, so a rigidly rotated site fits to 0.0 as the name implies.
  (c) Global-frame RMSD was computed even with no authoritative transform,
  comparing coordinates in unrelated frames; it is now None unless
  `target_transform` is supplied. (d) `_envelope_volume` returned 0.0 for fewer
  than four atoms and let SciPy `QhullError` escape on collinear or coplanar input;
  both are now reported as unavailable. A fifth fix in `difference_map_service.py`:
  masked target distances were stored as NaN, which the strict finiteness check in
  `core/difference_maps/__init__.py` rejected, so every partially mapped structure
  pair raised `ValueError: target_distances_angstrom must be finite`. Masked cells
  now hold a finite placeholder and `valid_mask` remains the sole authority.
- 2026-08-30 Task 1 gate output: `python -m pytest -q` reported 141 passed with no
  skips; `python -m pytest tests/unit/core/geometry tests/unit/application -q`
  reported 37 passed; Ruff passed over the touched paths; mypy reported success;
  `git diff --check` was clean.
- 2026-08-30 Task 1 remainder: two RED tests were added to
  `tests/unit/application/test_analysis_service.py` and confirmed failing before
  any production edit. (a) `AnalysisResult.mutation_count` counted every event
  whose kind was not "conserved", so a non-standard residue such as MSE — which
  has no canonical one-letter code and is therefore neither conserved nor a
  substitution — was counted as a change. A structure containing MSE compared
  against itself reported `mutation_count == 1`. The property now excludes
  "nonstandard" while those positions remain in `mutations` as evidence, so the
  Residues view still shows them. (b) `_structural_transform` returned a
  transform only when US-align drove the mapping, so every sequence-guided and
  manual comparison stored `transform=None`. That mattered more after the
  earlier Task 1 change gated global-frame site RMSD on an authoritative
  transform: without it, sequence comparisons could never report a global-frame
  site measurement at all. `_calculate_geometry` already computes a Kabsch fit
  through `superpose()` and reports every RMSD and per-residue displacement
  under it, so that fit is now stored as the authoritative transform for the
  non-US-align paths. `plugin/gui/qt_panel.py` already forwards
  `result.transform` into `calculate_site_metrics`, so the Sites view inherits
  the fix with no GUI change.
- 2026-08-30 Task 1 end-to-end check: an MSE-containing structure compared with
  itself now reports `mutation_count == 0`, still flags the MSE position as
  non-standard, and stores an identity transform with `det(rotation) == 1.0`.
- 2026-08-30 Task 1 gate output: `python -m pytest -q` reported 143 passed with
  no skips; the plan's focused selection reported 35 passed and
  `tests/unit/core/geometry tests/unit/application` reported 39 passed;
  `python -m ruff check .` passed; `python -m mypy src` reported success across
  109 source files; `git diff --check` was clean.
- 2026-08-30 Task 2 (typed availability and provenance): RED tests first covered
  immutable availability states, severity-coded diagnostics, recursive scientific
  parameter freezing, finite values, required physical units, deterministic
  artifact identities, and separation of wall-clock audit events. GREEN added
  `Availability`, `Diagnostic`, `AuditEvent`, and `MethodProvenance`, including
  canonical JSON, SHA-256 input validation, dimensional length/area/volume/angle
  checks, and a v0.3 compatibility view. Nanometres are rejected until explicit
  coordinate conversion exists rather than being silently mislabeled.
- 2026-08-30 Task 2 integration: typed provenance is stored independently from
  the legacy flat map and now survives `AnalysisResult`, `TargetAnalysis`,
  reference-vs-many conversion, project serialization, GUI handoff, and PyMOL
  bundle export. Bundle validation recomputes artifact IDs, rejects malformed or
  conflicting typed payloads, and binds target-scoped provenance keys to the
  manifest target IDs while allowing documented partial coverage.
- 2026-08-30 Task 2 review gate: two read-only review passes found and closed
  loss through multi-analysis, dimensionally incompatible units, unverified or
  conflicting bundle provenance, and target-key substitution. Final Terra and
  code reviews reported no remaining CRITICAL/HIGH blocker. Fresh gate output:
  `python -m pytest -q` reported 171 passed; `python -m ruff check .` passed;
  `python -m mypy src` reported success across 111 source files; and
  `git diff --check` was clean. Focused branch coverage was 91% for the new
  availability/diagnostic module and 84% for the new provenance module.
  `pip-audit` resolved the declared runtime, GUI, and chart dependencies and
  reported no known vulnerabilities; a filename-only secret-pattern scan found
  no matches. The global interpreter's unrelated `open-webui`/`onnxruntime`
  `pip check` mismatch is outside StructLens and is not part of the release env.
- 2026-08-31 Task 3 (parser enrichment and component identity): RED tests first
  covered bounded source capture, gzip expansion/truncation/trailing-data behavior,
  raw versus logical hashes, immutable selections and retained components, parser
  limits, multiple-model selection, altloc identity, PDB/mmCIF metadata, blank
  chains, legacy wrapper compatibility, and replacement of the original path
  after capture. GREEN added `SourceSnapshot`, typed parser contracts, strict
  component/structure invariants, conservative versioned classification, and a
  single snapshot-consuming normalization boundary.
- 2026-08-31 Task 3 scientific semantics: the canonical protein view contains
  polymer residues only, while every component in the selected model is retained
  as evidence and records whether it was selected for analysis. PDB does not
  fabricate label/entity/label-sequence identifiers; mmCIF preserves per-residue
  author, label, entity, formal-charge, and source-atom identity. Unknown HETATM
  and carbohydrates without CCD/connectivity evidence remain `other`; MSE remains
  a modified polymer; and `AltlocPolicy.ALL` fails closed until conformer-aware
  analyses exist. Legacy loaders still retain all models and non-water HETATM.
- 2026-08-31 Task 3 review fixes: preflight now precedes both modern and legacy
  parsing; duplicate/non-preservable PDB serials and duplicate mmCIF atom-site IDs
  fail at the boundary; mmCIF row limits run before `MMCIF2Dict` materialization;
  valid zero-padded PDB/mmCIF model and atom IDs share canonical lookup identities
  while retaining raw mmCIF source IDs; absent mmCIF label chains remain absent;
  blank chains remain distinct from missing chains; and analyzed author/label/entity
  metadata is deduplicated without losing split label segments.
- 2026-08-31 Task 3 gate: final independent code and scientific reviews reported
  no remaining CRITICAL/HIGH blocker. Fresh output: `python -m pytest -q` reported
  236 passed; Ruff passed; mypy reported success across 117 source files;
  `git diff --check` passed. Combined branch coverage for the new/expanded Task 3
  model and parsing modules was 81%. `pip-audit` resolved the declared runtime,
  GUI, and chart dependency set with no known vulnerabilities, and the changed-file
  secret-pattern scan found no secret-like assignments or private keys.
- 2026-08-31 Task 4 (coordinate QC): RED tests first covered raw malformed PDB/mmCIF
  scalars, parser short-circuiting, semantic atom duplicates, altloc occupancy,
  missing backbone atoms, element/radius uncertainty, chain geometry, indexed
  overlaps, rigid-body invariance, sparse-input behavior, typed reports, and service
  provenance. The new raw audit keeps `AtomRecord` strict and retains every source
  row in lenient immutable snapshots; blocking errors raise `CoordinateQualityError`
  with the exact typed report before `PDBParser`/`MMCIFParser` can run. Direct PDB
  audit paths enforce byte limits before reading and atom/model/chain/residue limits
  before materializing residue snapshots. Repeated atom serials in distinct models
  remain valid while duplicate chemical identities in one model fail closed.
- 2026-08-31 Task 4 scientific semantics: the atomic gate reports non-numeric or
  non-finite coordinates, occupancy outside the declared inclusive interval,
  B-factors below the declared minimum, missing explicit elements, unknown radii,
  duplicate identities, altloc sums above 1.01, and missing polymer backbone atoms.
  Geometry uses permissive coordinate-sanity intervals (`C(i)-N(i+1) <= 2.0 Å` and
  `2.5–4.5 Å` for consecutive C-alpha atoms), not refinement validation. Heavy-atom
  overlaps use `cKDTree`, full VdW overlap `>= 0.40 Å`, deterministic pair ordering,
  and one versioned Bondi-style polymer table with an explicit Se extension.
  Hydrogens, same-residue pairs, adjacent residues, and probable disulfides are
  excluded; ligand/water/ion/metal contacts and source connectivity are outside this
  first screen. The provenance states that source connection records are not yet
  retained, so no MolProbity or crystallographic clash-score claim is made.
- 2026-08-31 Task 4 review/gate: delegated Luna implementation covered atomic QC,
  raw parser gating, geometry, clashes, and radii; Terra and Sol reviews defined the
  fail-closed boundary and scientific thresholds. All reported HIGH risks were
  addressed: no partial normalized structure, no post-Biopython duplicate loss, no
  element-name guessing, restricted clash scope, and explicit geometry limitations.
  A fresh gate reported 300 passing tests; Task 4 focused branch coverage was 95%
  (required 90%); Ruff passed; mypy passed across 125 source files; `git diff
  --check` passed; and `pip-audit .` reported no known vulnerabilities. Two final
  delegated re-review attempts were unavailable because the subagent quota was
  exhausted, so the primary agent completed the final diff/security/scientific
  audit locally. Next executable plan item is Task 5, canonical report and v0.3
  orchestration.
