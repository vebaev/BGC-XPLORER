#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TOOLS_DIR="${ROOT_DIR}/tools"
ARTS_DIR="${TOOLS_DIR}/arts"
REF_DIR="${ROOT_DIR}/db/arts"

mkdir -p "${TOOLS_DIR}"
mkdir -p "${REF_DIR}"

if [ ! -d "${ARTS_DIR}/.git" ]; then
  git clone https://github.com/ZiemertLab/ARTS "${ARTS_DIR}"
else
  git -C "${ARTS_DIR}" pull --ff-only
fi

if compgen -G "${ARTS_DIR}/reference/*.zip" >/dev/null; then
  unzip -o "${ARTS_DIR}/reference/*.zip" -d "${REF_DIR}"
fi

cat <<EOF
ARTS repository prepared at:
  ${ARTS_DIR}

Reference directory:
  ${REF_DIR}

Validate with:
  python3 scripts/check_databases.py
EOF
