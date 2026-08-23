# StructLens visual system

<!-- impeccable:design-schema 1 -->

StructLens is an operate-mode scientific instrument inside a dark or light PyMOL host. The surface uses restrained neutrals, one controlled blue accent for active controls, and semantic status colors paired with text labels and icons. Layout is table-first: compact section headings, dense residue tables, explicit units, and generous separation between workflow stages. Typography follows the host system UI for controls and a readable tabular face for residue data; monospace is reserved for identifiers and numeric provenance. Avoid decorative gradients, dashboard card mosaics, and color-only encodings. Every continuous scale has a legend and every long-running action exposes progress, cancellation, and a recoverable error state.

## Evidence Bench Qt surface

- The shell uses a graphite canvas (`#111827`), cooler sidebar/header layers (`#0d1726`, `#0c1421`), and a single blue action accent (`#2f7af8`). Semantic colors appear in text, table status, and PyMOL selections together; color is never the only encoding.
- A 228–270 px workflow rail keeps Project, Alignment, Mutations, Residues, Visualization, and Results visible. The main canvas is table-first, with a persistent header action, footer status line, indeterminate progress bar, and explicit Cancel control for asynchronous comparisons.
- The first viewport is the Project source/chain workflow. Results tables use non-editable rows, explicit Å/fraction units, and double-click focus into namespaced PyMOL selections. Visualization controls expose filter, color mode, representation, preset, visibility, radius, and a textual range legend.
- Qt is discovered lazily (PySide6, then PyQt5); PyMOL is injected as a command proxy. Closing/resetting the panel deletes only StructLens-owned selections.
