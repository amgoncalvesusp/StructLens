# StructLens v0.5.0

## Comparing variants

- Select a reference in **Project**, open the variant comparison dialog and add target files. Choose a model and chain for ambiguous targets. The queue produces a separate report for each successful pair. Inspect any completed pair or cancel the remaining work without losing completed results. The list lasts for the current session. Save individual pairs to reopen them later.
- **Results** counts substitutions, insertions and deletions. Click a count to open the matching residue rows. In **Residues**, search either structure by author chain and position (including insertion code), residue name or mutation notation. The evidence card follows the selection.
- **Sequences** shows a window of aligned residue pairs with explicit gaps. Selecting a column opens that position in the residue table.

## Sites and ligand contacts

- Choose a deposited ligand by its structure identifier, preview the reference residues in the radius and keep multiple named sites. Site measurements and per-residue evidence use the same resolved membership.
- **Sites & pockets** shows a contact table and a selectable 2D schematic. For a deposited ligand, the distance is the nearest measured heavy-atom pair within 4.5 Å. The distance alone does not identify a hydrogen bond or establish a binding mechanism.
- To review docking interactions, compare the single-model PDB receptor used in DockLens, then import the DockLens JSON result and choose its original pose file. StructLens verifies both file hashes, displays the chosen pose and preserves the DockLens interaction family and profile. These imported contacts are a session view; keep the original DockLens files to reopen them. DSV-like and LUNA-like labels identify DockLens geometry rules. No external LUNA calculation is performed by StructLens. [DockLens interaction definitions](https://github.com/amgoncalvesusp/docklens/tree/a7711f2e6dc3034428bad507537ae8f34ee5b5ac/normalizacao_interacoes_intermoleculares)

## Charts and PyMOL

- **Charts** distinguishes available datasets from those the report cannot provide. Click a profile point to inspect its residue. Missing observations break the plotted line.
- In the distance-difference map, a cell reports the target distance minus the reference distance in Å. The scale is centered on zero. Select a measured cell to inspect both residue positions. Missing cells stay masked in the figure and empty in the spreadsheet.
- Export a verified comparison as a native PyMOL command file or open its overview and selected-residue scenes. The scenes use the report’s transform, so PyMOL does not perform a second fit. The exported file embeds both complete source structures, including models and chains outside the displayed selection. PyMOL is a separate application.

## Installers

Choose the Windows x86-64 Setup, Linux amd64 Debian package or Linux x86-64 AppImage from the release assets. Each download has a matching SHA-256 file. The desktop installers include the chart renderer.
