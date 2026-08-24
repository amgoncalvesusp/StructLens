# StructLens-PyMOL interchange

StructLens and StructLens-PyMOL are separate products. StructLens performs all scientific analysis; the plugin only visualizes imported results. The optional integration is a validated `.structlens-pymol` ZIP bundle and does not use sockets, HTTP, WebSockets, or a shared Python interpreter.

## Bundle contract

Schema `1.0` contains a deterministic `manifest.json`, provenance, original structure files, correspondence/mutation/metric JSON, transforms, and visualization presets. The writer preserves source numbering and records transformations explicitly. It never embeds Python, executables, or arbitrary payloads. `validate_pymol_bundle` performs ZIP traversal checks, rejects executable entries, verifies required JSON and referenced structures, and rejects unsupported future schema majors.

The bundle can always be exported without PyMOL. `Open in PyMOL` remains an optional launch handoff; if PyMOL or the plugin is unavailable, the validated bundle is the portable fallback.

## Plugin repository and release

The canonical plugin source is `https://github.com/amgoncalvesusp/pymol-plugins/` under `structlens-pymol/`. Plugin releases use tags `structlens-pymol-vMAJOR.MINOR.PATCH` and ship `StructLens-PyMOL-vX.Y.Z.zip` plus a SHA-256 checksum. The plugin must not silently recompute mapping, RMSD, TM-score, or outlier status.
