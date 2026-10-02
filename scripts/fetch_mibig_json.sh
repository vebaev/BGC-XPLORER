#!/usr/bin/env bash
# Download the MIBiG 4.0 JSON entries used by the AI summary (closest known cluster:
# compound, organism, class, activity, evidence, reference). About 1 MB.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DB_DIR="${BGC_DB_ROOT:-${ROOT_DIR}/db}/mibig"
URL="https://dl.secondarymetabolites.org/mibig/mibig_json_4.0.tar.gz"
SHA256="986f1a69fb5eeb84fd7c0175ae921304c22fd5f8ed04c96d8c1c0c95a145837b"
ARCHIVE="${DB_DIR}/mibig_json_4.0.tar.gz"

mkdir -p "${DB_DIR}"
if [ ! -f "${ARCHIVE}" ] || ! echo "${SHA256}  ${ARCHIVE}" | sha256sum -c --status; then
  echo "Downloading MIBiG 4.0 JSON from ${URL}"
  curl -fsSL --retry 3 -o "${ARCHIVE}.part" "${URL}"
  echo "${SHA256}  ${ARCHIVE}.part" | sha256sum -c --status || {
    echo "MIBiG 4.0 JSON checksum mismatch" >&2
    rm -f "${ARCHIVE}.part"
    exit 1
  }
  mv "${ARCHIVE}.part" "${ARCHIVE}"
fi
rm -rf "${DB_DIR}/mibig_json_4.0"
tar -xzf "${ARCHIVE}" -C "${DB_DIR}"
echo "MIBiG 4.0 JSON ready: $(ls "${DB_DIR}/mibig_json_4.0" | wc -l) entries in ${DB_DIR}/mibig_json_4.0"
