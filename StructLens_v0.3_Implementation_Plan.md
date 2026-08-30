# STRUCTLENS v0.3.0 — Implementation Plan

This file records the v0.3.0 plan supplied in the implementation request. The
complete normative text is the user message in the task thread; this local
ledger copy keeps the execution order and release gates durable across context
compaction.

## Binding scientific capabilities

1. Multiple sequence alignment, alignment conservation, insertion columns,
   sequence logo, and provenance.
2. Structural neighborhoods and reference-normalized interaction differences.
3. Key-residue, ligand-radius, and residue-radius site definitions with
   global-frame and site-fitted RMSD, SASA, and atomic-envelope volume.
4. Distance-difference matrices and transformed structural displacement vectors.
5. Residue Evidence Card with sequence, structure, interactions, site, and
   evidence-quality sections; never an impact/damage/function score.

## Required contracts

- Preserve `ResidueId` source numbering and authoritative correspondence.
- Represent missing/unavailable values as `None`, never as zero.
- MSA column → sequence index → source `ResidueId`; preserve reference-gap and
  reference-relative insertion columns.
- Entropy excludes gaps and ambiguous residues from canonical amino-acid
  counts; gap and ambiguous fractions remain separate.
- Interaction comparison maps target residues to reference positions first.
- Heavy-atom hydrogen bonds are putative (`heavy_atom_geometry`) without
  explicit hydrogens.
- Site atomic-envelope volume is a convex-hull envelope, not cavity volume.
- Distance differences are internal and require no superposition.
- Vectors apply exactly one target→reference transform and match stored Cα
  displacement magnitudes.
- PyMOL consumes bundle results and never recomputes science.

## Required work order

1. Freeze architecture, v0.3 version, shared domain types, and service
   protocols.
2. Bundle MUSCLE and verify redistribution metadata.
3. Parallel-safe wave: MSA/conservation; interactions; sites; difference maps.
4. Review and integrate the four scientific modules.
5. Implement Evidence Card, schema v3, XLSX/charts, GUI, bundle, and plugin.
6. Build Windows Setup, Linux Debian, Linux AppImage, and plugin ZIP.
7. Run cross-platform parity, end-to-end fixture, documentation, and release
   automation checks.

## GUI navigation

`Project`, `Sequences`, `Structures`, `Residues`, `Sites`, `Charts`, `PyMOL`,
`Results`, `Export`.

Before GUI changes, use the available `impeccable` skill and repeat its final
critique after implementation.

## Release gates

- pytest, Ruff, and mypy pass.
- Golden tests cover MSA, conservation, interactions, sites, distance maps,
  vectors, and Evidence Cards.
- Windows/Linux parity and clean installs pass.
- US-align, MUSCLE, and FreeSASA work offline from installers.
- XLSX, JPEG, and genuine 300/600 dpi TIFF exports pass.
- `.structlens-pymol` v0.3 exports and imports in the plugin.
- Plugin interactions, sites, vectors, Evidence Inspector, reset, and image
  exports pass.
- Documentation, licenses, citations, checksums, and release notes are
  complete.

## Branches and artifacts

- Main branch: `feat/structlens-v0.3.0`, release tag `v0.3.0`.
- Plugin branch: `feat/structlens-pymol-v0.3.0`, tag
  `structlens-pymol-v0.3.0`.
- Main artifacts: `StructLens-v0.3.0-Windows-x86_64-Setup.exe`,
  `structlens_0.3.0_amd64.deb`,
  `StructLens-v0.3.0-Linux-x86_64.AppImage`, with SHA-256 files.
- Plugin artifacts: `StructLens-PyMOL-v0.3.0.zip` and
  `StructLens-PyMOL-v0.3.0.sha256`.

---

# Progress ledger

Each session appends one section here. The final section always states where
the next session resumes.

## Session 2026-08-30 — release gates restored and pinned

Work order item 7 (cross-platform parity, end-to-end fixture, documentation,
and release automation checks) was the open item. Four gates were failing or
unenforced; all four now hold.

### Gate: pytest, Ruff, and mypy pass — DONE (`1624bdb`)

Ruff reported eleven `UP042` errors and mypy was not runnable. Fixed by:

- Migrating the eleven `(str, Enum)` classes to `enum.StrEnum`. On Python 3.11
  a `str`/`Enum` mixin gives `str()` and `format()` different results;
  `StrEnum` makes both return the value. Verified safe first — no call site
  interpolates these enums, every use goes through `.value`.
- Narrowing the distance-map cell conversion so the `None` branch is expressed
  once (`qt_panel.py`), and coercing the legacy site `kind` keyword to `str`.
- Declaring `PySide6` as an optional untyped import in `[[tool.mypy.overrides]]`.

### Gate: golden tests over the seven scientific lanes — DONE (`d60a9a1`)

New `tests/golden/`. `pipeline.py` drives the real services over two checked-in
synthetic PDB fixtures; `expected.json` pins the numbers; `test_golden_v03.py`
asserts them per lane. Regenerate deliberately:

```
STRUCTLENS_REGEN_GOLDEN=1 python -m pytest tests/golden
```

Every lane was mutation-checked — a deliberate defect was injected into each
and the baseline confirmed to fail. The first fixture attempt was a pure rigid
translation, which made the distance-difference matrix identically zero and
caught nothing; the target now also displaces residue `A:5` by 1.2 Å so the
matrix carries real signal.

Building the fixtures surfaced three genuine defects, all fixed:

- `_envelope_volume` raised `QhullError` on collinear or coplanar site atoms
  instead of reporting the envelope unavailable, and returned `0.0` for sites
  with fewer than four atoms. The plan requires `None`, never zero.
- `scipy` backs that convex hull but was never declared, so site atomic
  envelope volume was silently unavailable in every clean install. Now a
  declared dependency.
- The conservation and entropy formula was duplicated in `msa_service` instead
  of reusing `core.msa.conservation.column_statistics`; the copy also let a
  fully conserved column report `-0.0` bits.

### Gate: release automation and cross-platform parity — DONE (`e53633b`, `8f3639e`)

`release.yml` built and published installers without ever running the gates, so
a red build could ship. New reusable `.github/workflows/ci.yml` runs Ruff,
mypy, and pytest across Ubuntu and Windows on Python 3.11 and 3.12 — that
matrix is also the parity evidence. Both installer jobs now list it under
`needs:`.

`test_qt_panel.py` skipped whenever no Qt binding was importable, which was
every environment including CI, so eighteen plugin-panel tests never ran
anywhere. CI now installs the PySide6 extra, runs Qt offscreen, and fails if
any test reports as skipped or if a test run rewrote the golden baseline.

### Coverage — 80% (`fae2c5a`)

Was 55%. Six modules with real branches had no test at all (PyMOL transform
validation, state-snapshot restore, MUSCLE platform resolver, FreeSASA
unavailable path, backend version table, plugin entry point). 153 tests, no
skips. Scientific lanes sit at 83–100%.

### Open findings for the next session

1. **Reference-only interactions vanish silently.** In
   `core/interactions/comparison.py`, a reference record whose residue has no
   target counterpart returns `None` from `key()` and is dropped, while the
   mirror case is reported as `target_only_unmapped`. In the golden fixture the
   reference PHE6 hydrophobic contact disappears rather than being reported.
   The plan lists four categories — conserved, gained, lost, target-only-
   unmapped — so this is a deliberate-looking asymmetry that changes documented
   science. It needs a decision, not a unilateral fix. The golden currently
   pins the existing behaviour, so the decision is visible in the diff.
2. **Golden fixture geometry is synthetic.** Detected hydrogen bonds sit at
   1.4–2.6 Å, well below a physiological 2.8–3.2 Å, so the golden pins the
   arithmetic rather than validating the chemistry. Raising the H-bond
   threshold from 3.5 to 2.8 Å does not move the fixture. Swap in a curated
   real structure pair if the thresholds themselves need regression cover.
3. **`qt_panel.py` is 1442 statements in one file**, far past the 800-line
   project limit, and its six sibling `*_page.py` modules are three-line
   re-export shims. Splitting the panel along those existing page boundaries is
   the obvious refactor, and it is now safe to attempt because the panel has
   77% test cover.

### Resume here

Gates 1, 2, and the automation/parity gates hold and are enforced in CI. Not
yet verified, in plan order:

- **Clean installs pass, and US-align, MUSCLE, and FreeSASA work offline from
  the installers.** `scripts/verify_release_backends.py` runs inside
  `release.yml`, but no installed artifact has been exercised this session.
  This needs an actual Windows and Linux install run.
- **XLSX, JPEG, and genuine 300/600 dpi TIFF exports.** `test_image_export.py`
  and `test_v03_exports_and_schema.py` cover the writers; the DPI claim on a
  real exported file is unverified.
- **`.structlens-pymol` v0.3 round trip in the plugin**, and the plugin's
  interactions, sites, vectors, Evidence Inspector, reset, and image exports.
- **Documentation, licenses, citations, checksums, and release notes.**
  `docs/implementation-log.md` stopped at 2026-08-23 and has been brought
  current; `RELEASE_NOTES_v0.3.0.md` has not been checked against the shipped
  behaviour.

Start with the export-fidelity gate — it is verifiable locally without a
release build, unlike the installer gates.
