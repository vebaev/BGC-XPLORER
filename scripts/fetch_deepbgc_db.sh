#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DB_DIR="${BGC_DB_ROOT:-${ROOT_DIR}/db}/deepbgc"

mkdir -p "${DB_DIR}"

DEEPBGC_BIN="$(command -v deepbgc 2>/dev/null || true)"
if [ -z "${DEEPBGC_BIN}" ] && [ -x /opt/conda/envs/deepbgc/bin/deepbgc ]; then
  DEEPBGC_BIN=/opt/conda/envs/deepbgc/bin/deepbgc
fi

if [ -n "${DEEPBGC_BIN}" ]; then
  echo "Running local DeepBGC downloader"
  (
    cd "${DB_DIR}"
    DEEPBGC_DOWNLOADS_DIR="${DB_DIR}" "${DEEPBGC_BIN}" download
  )
  exit 0
fi

if docker image inspect aktino/deepbgc:latest >/dev/null 2>&1; then
  echo "Running DeepBGC downloader through aktino/deepbgc:latest"
  docker run --rm \
    --user "$(id -u):$(id -g)" \
    -e DEEPBGC_DOWNLOADS_DIR=/db/deepbgc \
    -v "${DB_DIR}:/db/deepbgc" \
    aktino/deepbgc:latest \
    deepbgc download
  exit 0
fi

cat <<'EOF'
DeepBGC downloader was not executed because neither local `deepbgc` nor
Docker image `aktino/deepbgc:latest` is available.

Required next step:
  bash scripts/build_tool_images.sh
  bash scripts/fetch_deepbgc_db.sh
EOF
exit 1
