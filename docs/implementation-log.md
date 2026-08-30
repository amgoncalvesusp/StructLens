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
- 2026-08-24: Shipped the v0.3 scientific modules, schema 3.0 exports, GUI result routing, and the offline installer pipeline; see `docs/v0.3.md`.
- 2026-08-30: Restored the pytest/Ruff/mypy gate by migrating eleven `(str, Enum)` classes to `enum.StrEnum` and clearing two mypy errors.
- 2026-08-30: Added `tests/golden/` regression cover for MSA, conservation, interactions, sites, distance maps, vectors, and Evidence Cards, driven by two synthetic PDB fixtures and pinned in `tests/fixtures/golden/expected.json`. Each lane was mutation-checked.
- 2026-08-30: Fixed three defects the golden fixtures exposed: site atomic-envelope volume crashed on degenerate hulls and returned `0.0` below four atoms, `scipy` backed that hull without being declared, and the conservation/entropy formula was duplicated in `msa_service` instead of reusing `column_statistics`.
- 2026-08-30: Added `.github/workflows/ci.yml` running Ruff, mypy, and pytest on Ubuntu and Windows across Python 3.11 and 3.12, gated the installer jobs on it, and made the Qt plugin tests run offscreen instead of skipping. Coverage moved from 55% to 80% with no skipped tests.
