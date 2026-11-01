#!/usr/bin/env bash
set -euo pipefail

VERSION="${1:?usage: build_packages.sh VERSION PYINSTALLER_DIST}"
PYINSTALLER_DIST="${2:?usage: build_packages.sh VERSION PYINSTALLER_DIST}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
APPDIR="$ROOT/build/appimage/StructLens.AppDir"
DEBROOT="$ROOT/build/deb-root"
rm -rf "$APPDIR"
rm -rf "$DEBROOT"
mkdir -p "$APPDIR/opt/structlens" "$APPDIR/usr/bin" "$APPDIR/usr/share/applications" "$APPDIR/usr/share/icons/hicolor/256x256/apps"
mkdir -p "$DEBROOT/DEBIAN" "$DEBROOT/opt/structlens" "$DEBROOT/usr/bin" "$DEBROOT/usr/share/applications" "$DEBROOT/usr/share/icons/hicolor/256x256/apps"
cp -a "$PYINSTALLER_DIST/StructLens/." "$APPDIR/opt/structlens/"
cp -a "$PYINSTALLER_DIST/structlens" "$APPDIR/opt/structlens/structlens"
cp -a "$ROOT/licenses" "$APPDIR/opt/structlens/"
cp "$ROOT/THIRD_PARTY_NOTICES.txt" "$APPDIR/opt/structlens/"
cp "$ROOT/RELEASE_NOTES_v${VERSION}.md" "$APPDIR/opt/structlens/"
cp -a "$PYINSTALLER_DIST/StructLens/." "$DEBROOT/opt/structlens/"
cp -a "$PYINSTALLER_DIST/structlens" "$DEBROOT/opt/structlens/structlens"
cp -a "$ROOT/licenses" "$DEBROOT/opt/structlens/"
cp "$ROOT/THIRD_PARTY_NOTICES.txt" "$DEBROOT/opt/structlens/"
cp "$ROOT/RELEASE_NOTES_v${VERSION}.md" "$DEBROOT/opt/structlens/"
cp "$ROOT/src/structlens/plugin/assets/structlens_icon.png" "$APPDIR/usr/share/icons/hicolor/256x256/apps/structlens.png"
cp "$ROOT/src/structlens/plugin/assets/structlens_icon.png" "$APPDIR/structlens.png"
cp "$ROOT/src/structlens/plugin/assets/structlens_icon.png" "$DEBROOT/usr/share/icons/hicolor/256x256/apps/structlens.png"
cat > "$APPDIR/AppRun" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
exec "$HERE/opt/structlens/StructLens" "$@"
EOF
chmod +x "$APPDIR/AppRun"
cat > "$APPDIR/usr/bin/structlens" <<'EOF'
#!/usr/bin/env bash
exec /opt/structlens/structlens "$@"
EOF
chmod +x "$APPDIR/usr/bin/structlens"
cat > "$APPDIR/usr/share/applications/structlens.desktop" <<EOF
[Desktop Entry]
Type=Application
Name=StructLens
Comment=Integrated protein sequence and structure comparison
Exec=/opt/structlens/StructLens
Icon=structlens
Terminal=false
Categories=Science;Education;
EOF
cp "$APPDIR/usr/share/applications/structlens.desktop" "$DEBROOT/usr/share/applications/structlens.desktop"
cp "$APPDIR/usr/share/applications/structlens.desktop" "$APPDIR/StructLens.desktop"
cat > "$DEBROOT/DEBIAN/control" <<EOF
Package: structlens
Version: ${VERSION}
Section: science
Priority: optional
Architecture: amd64
Maintainer: Adriano Marques Gonçalves (UNIARA)
Depends: libc6, libegl1, libgl1, libdbus-1-3, libxkbcommon-x11-0, libxcb-cursor0
Description: Integrated protein sequence and structure comparison
 StructLens provides sequence/structure comparison and reproducible evidence exports.
EOF
cat > "$DEBROOT/usr/bin/structlens" <<'EOF'
#!/usr/bin/env bash
exec /opt/structlens/structlens "$@"
EOF
chmod +x "$DEBROOT/usr/bin/structlens"
mkdir -p "$ROOT/release"
dpkg-deb --build --root-owner-group "$DEBROOT" "$ROOT/release/structlens_${VERSION}_amd64.deb" >/dev/null
APPIMAGETOOL="$ROOT/build/appimagetool.AppImage"
curl --fail --location --retry 3 -o "$APPIMAGETOOL" "https://github.com/AppImage/appimagetool/releases/download/1.9.1/appimagetool-x86_64.AppImage"
echo 'ed4ce84f0d9caff66f50bcca6ff6f35aae54ce8135408b3fa33abfc3cb384eb0  appimagetool.AppImage' | (cd "$ROOT/build" && sha256sum --check --strict)
chmod +x "$APPIMAGETOOL"
RUNTIME="$ROOT/build/runtime-x86_64"
curl --fail --location --retry 3 -o "$RUNTIME" "https://github.com/AppImage/type2-runtime/releases/download/20251108/runtime-x86_64"
echo '2fca8b443c92510f1483a883f60061ad09b46b978b2631c807cd873a47ec260d  runtime-x86_64' | (cd "$ROOT/build" && sha256sum --check --strict)
ARCH=x86_64 APPIMAGE_EXTRACT_AND_RUN=1 "$APPIMAGETOOL" --runtime-file "$RUNTIME" "$APPDIR" "$ROOT/release/StructLens-v${VERSION}-Linux-x86_64.AppImage"
(cd "$ROOT/release" && sha256sum "structlens_${VERSION}_amd64.deb") > "$ROOT/release/structlens_${VERSION}_amd64.deb.sha256"
(cd "$ROOT/release" && sha256sum "StructLens-v${VERSION}-Linux-x86_64.AppImage") > "$ROOT/release/StructLens-v${VERSION}-Linux-x86_64.AppImage.sha256"
