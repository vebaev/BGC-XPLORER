#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DB_DIR="${BGC_DB_ROOT:-${ROOT_DIR}/db}"

mkdir -p "${DB_DIR}/antismash"
mkdir -p "${DB_DIR}/deepbgc"
mkdir -p "${DB_DIR}/arts/actinobacteria"
mkdir -p "${DB_DIR}/eggnog"

PYTHON_BIN="${BGC_PYTHON:-$(command -v python3)}"
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
  echo "ARTS: automatic preparation is unavailable for the selected taxon." >&2
  echo "Populate ${DB_DIR}/arts/actinobacteria and start again." >&2
fi

cat <<'EOF'
ARTS reference files are not auto-downloaded here because the workflow expects a taxon-specific
precomputed reference set, for example:

  db/arts/actinobacteria

Place the required ARTS reference files in that directory, then validate with:

  python3 scripts/check_databases.py
EOF
