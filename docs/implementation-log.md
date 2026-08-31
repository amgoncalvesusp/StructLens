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
- 2026-08-31 Task 5 (canonical report and v0.3 orchestration): closed the
  report-layer contracts around immutable snapshots, deterministic report IDs,
  typed section availability, site-definition persistence, and report-driven
  visualization compatibility. `AnalysisReport` now excludes local path and
  basename presentation metadata from canonical scientific identity, so renaming
  identical sources does not change report bytes or hashes. The orchestration path
  now preserves content-addressed `site_definitions`, records them in typed
  provenance, distinguishes unresolved sites as `not_detected` instead of a false
  zero-valued success, and propagates failed interaction availability into
  Evidence Card section quality instead of silently marking interactions present.
  `VisualizationService` now accepts both legacy `AnalysisResult` objects and the
  canonical `AnalysisReport`, returning an empty selection when a report has no
  analysis section rather than raising an attribute error.
- 2026-08-31 Task 5 gate: `python -m pytest -q` reported 350 passed. Focused
  verification used `python -m coverage run --branch -m pytest -p no:cov
  tests/unit/core/reports/test_models.py tests/unit/application/test_report_service.py`
  followed by `python -m coverage report -m --include="src/structlens/core/reports/*,src/structlens/application/report_service.py,src/structlens/application/visualization_service.py"`,
  yielding 45 passing focused tests and 93% total branch coverage across the new
  report modules (`report_service.py` 92%, `visualization_service.py` 100%,
  `core/reports/models.py` 93%). `python -m ruff check` passed on the touched
  report files and tests; `python -m mypy src` reported success across 128 source
  files; `git diff --check` was clean apart from Git's existing LF→CRLF warnings;
  and `python -m pip_audit .` reported no known vulnerabilities. `pytest-cov`
  triggered a Windows-local NumPy reimport collection error in this environment,
  so focused branch evidence was collected with `coverage run -p no:cov` instead.
- 2026-08-31 Task 5 post-gate refactor: split the report orchestration internals
  into `report_evidence.py`, `report_failures.py`, `report_geometry.py`,
  `report_input.py`, `report_provenance.py`, and `core/reports/serialization.py`
  so the report boundary and immutable contracts stay below the file-size limit
  and remain easier to review. Added regression checks for centroid-displacement
  missingness without a transform, nonphysical interaction thresholds, rigid-body
  transform validation, schema/comparison metadata, and a headless GUI
  `AnalysisReportController` boundary. Fresh focused verification reported
  `86 passed` for `tests/unit/core/reports`,
  `tests/unit/application/test_report_service.py`,
  `tests/unit/application/test_v03_services.py`,
  `tests/unit/core/test_v03_contracts.py`, and
  `tests/unit/plugin/test_gui_model.py`; `python -m ruff check` passed on the
  touched files; `python -m mypy src` reported success across 134 source files;
  and `git diff --check` remained clean apart from the existing LF→CRLF warnings.
  A fresh full-suite run now stops at collection in
  `tests/unit/core/pockets/test_models.py` and `test_radii.py` with
  `ModuleNotFoundError: No module named 'structlens.core.pockets'`, which is the
  expected RED starting point for Task 6 rather than a Task 5 regression.
- 2026-08-31 Task 6 (pocket geometry primitives): RED started from the new
  `tests/unit/core/pockets/` package with `ModuleNotFoundError` for
  `structlens.core.pockets`. GREEN added `core/pockets/__init__.py`,
  `models.py`, `radii.py`, and `geometry.py` with immutable
  `PocketGeometrySettings`, canonical-hash `AlphaSphere` and `PocketCandidate`,
  a conservative pocket-only Bondi-style radii table with an explicit selenium
  extension, and a typed `tetrahedron_circumsphere()` result that distinguishes
  malformed input (`ValueError`) from non-finite or degenerate simplices
  (`Diagnostic` with no NaN payload). Added/expanded tests cover deterministic
  identity under input reordering, residue and atom canonical ordering, invalid
  settings/contracts, rigid-body and permutation invariance, explicit tolerance
  handling for nearly degenerate tetrahedra, and rare linear-algebra failure
  paths. `tests/unit/core/pockets/__init__.py` was added to namespace the new
  tests and prevent basename collisions with existing `quality` and `reports`
  test modules during root collection.
- 2026-08-31 Task 6 gate: `python -m pytest tests/unit/core/pockets/test_models.py
  tests/unit/core/pockets/test_geometry.py tests/unit/core/pockets/test_radii.py
  --cov=src/structlens/core/pockets --cov-branch -q` reported 30 passed with
  96% total coverage across the new package (`geometry.py` 92%,
  `models.py` 99%, radii and package exports 100%). `python -m ruff check
  src/structlens/core/pockets tests/unit/core/pockets` passed after import
  normalization, `python -m mypy src/structlens/core/pockets` reported success,
  and `git diff --check` remained clean apart from the repository's existing
  LF→CRLF warnings. A fresh full-suite run on the shared worktree reported
  `391 passed, 2 failed`; both failures are outside Task 6 and come from
  concurrent uncommitted report-layer edits in
  `src/structlens/core/reports/models.py` and related application/tests, where
  default `AnalysisReport` availability no longer matches the two legacy tests
  `tests/unit/application/test_report_service.py::test_pairwise_report_rejects_a_selection_with_multiple_protein_chains`
  and `tests/unit/plugin/test_gui_model.py::test_report_controller_delivers_one_canonical_artifact_to_fake_widget`.
  The new `core.pockets` package itself remained green in the same workspace.
- 2026-08-31 Task 5/6 hardening after Task 6: closed the residual integration
  gap between canonical reports and the new pocket-ready branch. The report path
  now materializes exact `SourceSnapshot` bytes into private short-lived files
  only for backend adapters that still require filesystem paths, then removes
  them before serialization so no ephemeral path leaks into report identity.
  Pairwise selection diagnostics are projected back into the per-input QC reports,
  QC-runner failures are contained as `numerical_failure` without entering
  downstream analysis, evidence cards mark structure as unavailable when no
  authoritative displacement exists, ligand residue-name aliases are exposed only
  when unambiguous, and `AnalysisResult.transform` now always stores the strict
  Kabsch fit used by StructLens metrics rather than an arbitrary backend-native
  matrix. `AnalysisReport` also rejects incoherent input-quality/vector
  availability states explicitly. Fresh gate output:
  `python -m pytest tests/unit/application/test_report_service.py tests/unit/plugin/test_gui_model.py tests/unit/core/reports/test_models.py tests/unit/application/test_report_evidence.py tests/unit/application/test_analysis_service.py -q`
  reported 65 passed; focused `coverage run --branch -m pytest -p no:cov ...`
  plus `coverage report` measured `report_service.py` at 93%,
  `report_materialization.py` at 100%, `report_evidence.py` at 95%,
  `report_input.py` at 96%, `report_failures.py` at 100%, and
  `core/reports/models.py` at 94%; `python -m ruff check ...` passed;
  `python -m mypy src` reported success across 139 source files;
  `python -m pytest -q` reported 393 passed; and `python -m pip_audit .`
  reported no known vulnerabilities.
- 2026-08-31 Task 5 final structural-scope correction: a second scientific
  review showed that passing the entire captured source file to US-align could
  silently reintroduce unselected models/chains even though parsing and QC used
  only the user's selection. Structural backend materialization therefore moved
  into `USAlignAdapter`, where it runs only after STRUCTURE/AUTO has actually
  chosen structural mapping. The adapter writes a private, ordered C-alpha PDB
  trace built from exactly the normalized selected `ProteinChain`; it never
  reads the mutable original path, rejects missing/duplicate C-alpha atoms or
  PDB field overflow instead of shifting residue indices, and deletes the trace
  after subprocess completion. A multi-model/two-chain regression replaces the
  original file after snapshot capture and proves that only model 2, chain B
  reaches structural analysis. The US-align subprocess regression verifies
  selected-trace content, safe argument-list invocation, and cleanup.
- 2026-08-31 post-scope gate: 68 focused report/analysis/US-align/GUI tests
  passed; the integration suite excluding the newly introduced Task 7 RED tests
  reported 381 passed; Ruff passed on the touched paths; and mypy succeeded
  across 138 source files. The only root-collection errors now belong to Task 7
  RED tests awaiting the blind pocket detector implementation.
- 2026-08-31 Task 7 (blind alpha-sphere detector): RED synthetic tests first
  covered sealed/open and separated cavities, input reordering, rigid transforms,
  coordinate jitter, underspecified/unknown-radius inputs, non-vertex VDW
  intrusion, Qhull failures, candidate caps, progress/cancellation, and resource
  limits. GREEN added bounded SciPy Delaunay tessellation, deterministic
  union-find clustering, transparent geometric ranking components, and an
  application service with explicit availability, diagnostics, counts, and
  provenance. A scientific review then found six HIGH defects before acceptance:
  atom-centre radii were stored instead of physical VDW clearance; all-atom and
  grid work were quadratic; exposure checked only sphere centres; provenance did
  not bind the input selection; chain selection could be widened; and invalid or
  cancelled direct results collapsed to `not_detected`. The corrected detector
  now stores and range-filters exact `distance(center, atom) - vdw_radius`
  clearance, rejects closer non-vertex surfaces through `cKDTree` queries, derives
  lining residues from surface clearance, rasterizes only atom-local grid boxes,
  flood-fills boundary solvent once, tests every free candidate-region cell, and
  enforces separate simplex, clearance-check, grid-cell, and raster-work budgets.
  It reads only normalized selected-chain/model residue records (therefore also
  the selected altloc policy), binds full selection identity/details into method
  provenance, and keeps invalid input, cancellation, numerical failure, valid no
  detection, and available candidates distinct.
- 2026-08-31 Task 7 gate: `python -m pytest tests/unit/core/pockets
  tests/unit/application/test_pocket_service.py --cov=src/structlens/core/pockets
  --cov-branch --cov-report=term-missing --cov-fail-under=90 -q` reported 98
  passed and 93.06% total branch-aware coverage. `python -m pytest -q` reported
  462 passed; `python -m ruff check .` passed; `python -m mypy src` succeeded
  across 142 source files; `python -m pip_audit .` reported no known
  vulnerabilities; and `git diff --check` was clean apart from existing
  LF-to-CRLF notices. Fresh independent Terra scientific/integration and
  code/security re-reviews found no remaining CRITICAL/HIGH issue. The next
  accepted item is Task 8, pocket free-volume measurement.
- 2026-08-31 Task 8 (pocket free-volume measurement): GREEN now includes
  `src/structlens/core/pockets/volume.py`,
  `src/structlens/core/pockets/volume_models.py`,
  `src/structlens/application/pocket_service.py`, and the corresponding focused
  tests. The service measures bounded coarse/fine free volume for one
  `PocketCandidate` using only the selected primary polymer scope plus
  explicitly retained ligand/ion/other components under the chosen
  `component_exclusion_policy`, and binds selection/candidate/sphere settings
  into deterministic method provenance. Late hardening added a regression for
  retained ligand altloc selection and two contract fixes in
  `volume_models.py`: `PocketVolumeResult` now rejects a short
  `grid_origin_xyz` with a clean `ValueError` instead of leaking `IndexError`,
  and it rejects non-`PocketVolumeSettings` `settings` objects instead of
  silently accepting invalid state.
- 2026-08-31 Task 8 gate: `python -m pytest tests/unit/core/pockets/test_volume.py
  tests/unit/application/test_pocket_volume_service.py --cov=src/structlens/core/pockets
  --cov-branch --cov-report=term-missing -q` reported 54 passed; within the
  touched scientific modules `volume.py` reached 90% branch-aware coverage and
  `volume_models.py` reached 96%. `python -m pytest -q` then reported
  `516 passed in 6.59s`. `python -m ruff check .` passed. `python -m mypy src`
  passed across 144 source files. `git diff --check` remained clean apart from
  the repository's existing LF-to-CRLF warnings on
  `src/structlens/application/pocket_service.py` and
  `src/structlens/core/pockets/__init__.py`.
- 2026-08-31 Task 8 dependency audit note: `python -m pip_audit` on the shared
  interpreter reported 232 known vulnerabilities in 39 installed packages.
  The output is not attributable to the new Task 8 code: many flagged packages
  (for example `aiohttp`, `anthropic`, `chromadb`, `open-webui`, `torch`) are
  not StructLens dependencies, while the one directly relevant project
  dependency in the report was `Pillow 12.1.1`, and StructLens already declares
  `Pillow>=12.3.0` in `pyproject.toml`. Treat that audit result as a local
  environment/release-preparation follow-up for Task 17/18 rather than a Task 8
  implementation blocker. The next accepted item is Task 9, ligand support and
  focused pockets.
