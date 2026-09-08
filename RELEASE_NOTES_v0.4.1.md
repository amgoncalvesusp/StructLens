# StructLens v0.4.1

This patch release improves the StructLens desktop GUI workflow and keeps all
scientific calculations and report data in the canonical application services.

## GUI

- Uses one clear Compare action with pending-report state and safe busy/export
  handling.
- Improves ordered page navigation, responsive scrolling, and collapsible
  advanced controls at 100% and 125% display scale.
- Restores manual pair selection after navigation and preserves legacy report
  workflows.
- Shows the PyMOL bundle writer as unavailable and disabled in this release;
  table and image exports remain available, and the legacy handoff is preserved.

## Validation

Native Windows smoke validation covered all nine GUI pages at 800x600 and
125% display scale, including real compare, JSON/XLSX/CSV export, navigation,
pending state, busy/export safety, and manual pair restoration.
