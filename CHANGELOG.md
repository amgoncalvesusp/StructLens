# Changelog

## 0.1.4

- Switched native setup builds to the compact PyQt5 runtime, keeping PySide6 available through the `gui-pyside6` extra for plugin deployments.
- Added binding selection through `STRUCTLENS_QT_BINDING` for deterministic plugin and standalone environments.

## 0.1.3

- Optimized native GUI packaging and standalone startup flow.

## 0.1.2

- Added the standalone Evidence Bench desktop GUI with file-based PDB/mmCIF loading, comparison, evidence tables, and exports.
- Added a multi-resolution StructLens icon to the desktop application, Windows setup, and Linux desktop entry.
- Added GUI-aware Windows and Linux setup builds while preserving the PyMOL plugin workflow.

## 0.1.1

- Added native Windows and Linux setup installers for double-click installation.

## 0.1.0

- Credited Adriano Marques Gonçalves (UNIARA) as the software author and added Windows/Linux release bundles.
- Added normalized PDB/mmCIF/FASTA parsing and explicit residue mapping.
- Added sequence alignment, Kabsch geometry, mutation descriptors, US-align adapter, application services, CLI, project JSON, and XLSX/CSV/JSON exports.
- Added optional namespaced PyMOL integration and six-section English GUI model.
