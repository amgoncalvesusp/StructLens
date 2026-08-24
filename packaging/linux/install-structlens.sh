#!/usr/bin/env bash
set -euo pipefail

BUNDLE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python3}"
WHEEL="$(find "$BUNDLE_ROOT" -maxdepth 1 -type f -name 'structlens-*.whl' -print -quit)"

if [[ -z "$WHEEL" ]]; then
  echo "The StructLens wheel was not found beside this installer." >&2
  exit 1
fi

echo "Installing StructLens from $(basename "$WHEEL")..."
"$PYTHON_BIN" -m pip install --upgrade "$WHEEL"
echo "Installation complete."
echo "Run structlens-gui to open the desktop application."
echo "Use Export for PyMOL from the app to generate a validated .structlens-pymol bundle."
