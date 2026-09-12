#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DB_DIR="${BGC_DB_ROOT:-${ROOT_DIR}/db}/eggnog"

mkdir -p "${DB_DIR}"

EGGNOG_DOWNLOADER="$(command -v download_eggnog_data.py 2>/dev/null || true)"
if [ -z "${EGGNOG_DOWNLOADER}" ] && [ -x /opt/conda/envs/eggnog/bin/download_eggnog_data.py ]; then
  EGGNOG_DOWNLOADER=/opt/conda/envs/eggnog/bin/download_eggnog_data.py
fi

if [ -n "${EGGNOG_DOWNLOADER}" ]; then
  echo "Running local eggNOG-mapper downloader"
  "${EGGNOG_DOWNLOADER}" --data_dir "${DB_DIR}" -y
  exit 0
fi

if docker image inspect aktino/eggnog:latest >/dev/null 2>&1; then
  echo "Running eggNOG-mapper downloader through aktino/eggnog:latest"
  docker run --rm \
    --user "$(id -u):$(id -g)" \
    -v "${DB_DIR}:/db/eggnog" \
    aktino/eggnog:latest \
    download_eggnog_data.py --data_dir /db/eggnog -y
  exit 0
fi

cat <<'EOF'
eggNOG-mapper downloader was not executed because neither local `download_eggnog_data.py`
nor Docker image `aktino/eggnog:latest` is available.

Required next step:
  bash scripts/build_tool_images.sh
  bash scripts/fetch_eggnog_db.sh
EOF
exit 1
