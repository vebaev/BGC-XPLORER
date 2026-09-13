#!/usr/bin/env bash
set -euo pipefail

DB_ROOT="${BGC_DB_ROOT:-/db}"
DB_TYPE="${BAKTA_DB_TYPE:-light}"
PYTHON_BIN="${BGC_PYTHON:-/opt/conda/bin/python}"
BAKTA_DB_BIN="${BAKTA_DB_BIN:-/opt/conda/envs/bakta/bin/bakta_db}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ "${DB_TYPE}" != "light" && "${DB_TYPE}" != "full" ]]; then
  echo "BAKTA_DB_TYPE must be 'light' or 'full', got '${DB_TYPE}'." >&2
  exit 2
fi

mkdir -p "${DB_ROOT}/bakta" "${DB_ROOT}/.locks"
exec 8>"${DB_ROOT}/.locks/bakta-${DB_TYPE}.lock"
flock -w "${DATABASE_LOCK_TIMEOUT:-3600}" 8

CHECK=("${PYTHON_BIN}" "${SCRIPT_DIR}/check_bakta_database.py" --db-root "${DB_ROOT}" --type "${DB_TYPE}")
if "${CHECK[@]}"; then
  exit 0
fi

echo "Downloading the latest Bakta-compatible ${DB_TYPE} database into ${DB_ROOT}/bakta..."
"${BAKTA_DB_BIN}" download --output "${DB_ROOT}/bakta" --type "${DB_TYPE}"
"${CHECK[@]}"

