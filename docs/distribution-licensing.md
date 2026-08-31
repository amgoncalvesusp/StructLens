# Distribution and licensing policy

StructLens source code is MIT licensed. Installers are not: each one bundles
several separately licensed components, and the MIT licence never governs an
installer on its own.

This is conservative engineering compliance, not legal advice.

## Components redistributed in official artifacts

| Component | Role | Licence | How it is distributed |
|---|---|---|---|
| StructLens | The application itself | MIT | Source and frozen bytecode |
| PySide6 (Qt for Python) | Default GUI binding | LGPL-3 | Frozen into GUI artifacts |
| PyQt5 | Opt-in GUI binding | GPL-3 | **Never** frozen into official artifacts |
| MUSCLE 5.3 | Multiple sequence alignment | GPL-3 | Separate executable plus corresponding source |
| US-align | Structural superposition | Permissive upstream notice | Compiled from pinned source |
| FreeSASA | Solvent-accessible surface area | MIT | Python package and C library |
| Biopython, NumPy, SciPy, OpenPyXL, Pillow | Scientific and I/O libraries | Own permissive licences | Installed as dependencies |

## Qt binding policy

Official frozen GUI artifacts use **PySide6**, which is LGPL-3. The `gui` extra
resolves to PySide6, and the release workflow passes `--exclude-module PyQt5`
so no GPL Qt binding can be pulled into a distributed binary.

PyQt5 remains available as the explicitly named `gui-pyqt5` extra for users who
prefer it. Choosing it makes the resulting combined work GPL, which is why it is
opt-in and never the default binary path.

### LGPL obligations for the frozen PySide6 build

Freezing an LGPL library into a single-file binary is permitted only if the
recipient can replace that library. For every official GUI artifact we therefore:

- state that PySide6 and Qt are LGPL-3 and identify their upstream sources;
- ship the PySide6 and Qt licence texts inside the installer;
- provide the **relinking route** — the artifact is built with PyInstaller from
  published source, and the documented build command plus pinned dependency set
  lets a recipient rebuild it against a modified PySide6/Qt;
- impose no term that forbids reverse engineering for debugging those
  modifications.

## MUSCLE and the GPL-3 corresponding-source obligation

MUSCLE 5.3 is invoked as a **separate executable**. It is never linked into,
copied into, or derived from StructLens source, so it does not make StructLens a
derivative work. Redistributing the binary still carries GPL-3 obligations in
full:

1. identify MUSCLE as a separate component with its copyright intact;
2. include the complete GPL-3 text — see `licenses/MUSCLE-LICENSE.txt`;
3. provide the exact **corresponding source** for the redistributed binary. The
   release workflow downloads the tagged upstream archive
   (`muscle-5.3-source.tar.gz`) next to the binary so the source travels inside
   every artifact that carries it;
4. impose no further restriction on the rights GPL-3 grants for that binary.

If any artifact cannot satisfy all four points reproducibly, MUSCLE is omitted
from that artifact and treated as a user-supplied optional backend. Alignment
then reports the backend as unavailable rather than silently degrading.

## Release gate

`tests/unit/test_packaging_metadata.py` enforces the mechanical parts of this
policy: the default GUI extra is the LGPL binding, PyQt5 stays opt-in, the
release workflow never freezes PyQt5, this document names every redistributed
component, and the MUSCLE notice carries the full GPL-3 text rather than a
promise of one.
