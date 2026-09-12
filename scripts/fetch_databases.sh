/bin/bash: warning: setlocale: LC_ALL: cannot change locale (C.UTF-8)
#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DB_DIR="${BGC_DB_ROOT:-${ROOT_DIR}/db}"

mkdir -p "${DB_DIR}/.locks"
exec 9>"${DB_DIR}/.locks/database-bootstrap.lock"
if ! flock -w "${DATABASE_LOCK_TIMEOUT:-3600}" 9; then
  echo "Timed out waiting for the database bootstrap lock under ${DB_DIR}/.locks" >&2
  exit 1
fi

mkdir -p "${DB_DIR}/antismash"
mkdir -p "${DB_DIR}/deepbgc"
ARTS_REFERENCE="${ARTS_REFERENCE:-actinobacteria}"
mkdir -p "${DB_DIR}/arts"
mkdir -p "${DB_DIR}/eggnog"

PYTHON_BIN="${BGC_PYTHON:-$(command -v python3)}"
# Determine missing resources only after obtaining the lock. Another container
# may have completed the database installation while this process was waiting.
MISSING="$(${PYTHON_BIN} "${ROOT_DIR}/scripts/check_databases.py" \
  --manifest "${ROOT_DIR}/db/manifest.yaml" --db-root "${DB_DIR}" --missing-only)"

run_if_missing() {
  local database="$1"
  local downloader="$2"
  if printf '%s\n' "${MISSING}" | grep -qx "${database}"; then
    BGC_DB_ROOT="${DB_DIR}" bash "${ROOT_DIR}/scripts/${downloader}"
  else
    echo "${database}: already valid; download skipped."
  fi
}

run_if_missing antismash fetch_antismash_db.sh
run_if_missing deepbgc fetch_deepbgc_db.sh
run_if_missing eggnog fetch_eggnog_db.sh
if printf '%s\n' "${MISSING}" | grep -qx arts; then
  BGC_DB_ROOT="${DB_DIR}" ARTS_REFERENCE="${ARTS_REFERENCE}" \
    bash "${ROOT_DIR}/scripts/fetch_arts.sh"
else
  echo "arts: already valid; download skipped."
fi
