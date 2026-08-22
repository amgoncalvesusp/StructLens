# Product

<!-- impeccable:product-schema 1 -->

## Platform

adaptive

## Users

Inferred from the implementation plan: structural biologists and computational researchers comparing homologous protein structures in PyMOL and from files.

## Product Purpose

StructLens determines explicit residue correspondence, mutation descriptors, global structural metrics, and residue-level differences, then exposes those results reproducibly in a PyMOL workflow.

## Positioning

The authoritative scientific object is a persisted residue correspondence table. PyMOL is a visualization backend and never the source of scientific truth.

## Operating Context

Users load PDB/mmCIF/FASTA inputs or existing PyMOL objects, choose reference and target chains, select Auto/Sequence/Structure/Manual mapping, inspect results, highlight residues, and export reproducible tables and publication figures.

## Capabilities and Constraints

Inferred from the supplied plan: English-only v0.1; explicit residue numbering including insertion codes; strict and refined RMSD remain separate; Cα displacement is not called residue RMSD; unsupported functional claims are prohibited; PyMOL imports are isolated to integration/plugin modules.

## Brand Commitments

The supplied StructLens icon is the approved visual identity and must be packaged at `src/structlens/plugin/assets/structlens_icon.png` when the plugin surface is built.

## Evidence on Hand

The supplied implementation plan and approved icon are the available product assets. No external customer, benchmark, pricing, or biological efficacy claims may be fabricated.

## Product Principles

- Make correspondence explicit and inspectable.
- Preserve strict scientific measurements and exclusions.
- Explain every workflow choice in plain scientific English.
- Keep visualization reversible and namespaced.
- Fail clearly when evidence is unavailable.

## Accessibility & Inclusion

Inferred from the plan: colors must not carry the only meaning, units and legends must be visible, keyboard navigation should remain usable where Qt supports it, and long-running analysis must not freeze the host UI.
