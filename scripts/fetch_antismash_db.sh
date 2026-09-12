#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DB_DIR="${BGC_DB_ROOT:-${ROOT_DIR}/db}/antismash"
BIN_DIR="${ROOT_DIR}/.local/bin"
MPL_DIR="${ROOT_DIR}/.cache/matplotlib"

mkdir -p "${DB_DIR}"
mkdir -p "${BIN_DIR}"
mkdir -p "${MPL_DIR}"

echo "Preparing antiSMASH database download helper"
HELPER_TMP="${BIN_DIR}/download_antismash_databases.tmp.$$"
trap 'rm -f "${HELPER_TMP}"' EXIT INT TERM
curl -fsSL "https://dl.secondarymetabolites.org/releases/latest/download_antismash_databases_docker" \
  -o "${HELPER_TMP}"
mv "${HELPER_TMP}" "${BIN_DIR}/download_antismash_databases"
chmod +x "${BIN_DIR}/download_antismash_databases"

echo "Downloading antiSMASH databases into ${DB_DIR}"
MPLCONFIGDIR="${MPL_DIR}" "${BIN_DIR}/download_antismash_databases" "${DB_DIR}"
