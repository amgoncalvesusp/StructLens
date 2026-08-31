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
