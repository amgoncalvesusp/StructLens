#!/usr/bin/env bash
# Install missing system libraries on ephemeral GitHub Ubuntu runners.
set -euo pipefail

missing=()
for package in "$@"; do
  if [[ "$(dpkg-query -W -f='${Status}' "$package" 2>/dev/null || true)" != "install ok installed" ]]; then
    missing+=("$package")
  fi
done
if (( ${#missing[@]} == 0 )); then
  exit 0
fi

. /etc/os-release
source_file="${RUNNER_TEMP:?}/structlens-ubuntu.sources"
cat > "$source_file" <<EOF
Types: deb
URIs: https://archive.ubuntu.com/ubuntu
Suites: ${VERSION_CODENAME} ${VERSION_CODENAME}-updates ${VERSION_CODENAME}-security
Components: main universe
Signed-By: /usr/share/keyrings/ubuntu-archive-keyring.gpg
EOF
options=(-o "Dir::Etc::sourcelist=$source_file" -o Dir::Etc::sourceparts=-
         -o Acquire::Retries=1 -o Acquire::https::Timeout=20)
sudo apt-get "${options[@]}" -o APT::Update::Error-Mode=any update
sudo apt-get "${options[@]}" install --no-install-recommends --yes "${missing[@]}"
