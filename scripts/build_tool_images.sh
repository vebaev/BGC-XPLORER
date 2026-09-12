#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

docker build -t aktino/gecco:latest "${ROOT_DIR}/containers/gecco"
docker build -t aktino/deepbgc:latest "${ROOT_DIR}/containers/deepbgc"
docker build -t aktino/eggnog:latest "${ROOT_DIR}/containers/eggnog"

cat <<'EOF'
Built local workflow images:
  aktino/gecco:latest
  aktino/deepbgc:latest
  aktino/eggnog:latest

ARTS is handled separately because the official repository includes its own Docker/runtime assets and reference set handling.
EOF
