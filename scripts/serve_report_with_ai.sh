#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

REPORT="results/Soil_1/report/Soil_1.html"
if [[ ! -f "${REPORT}" ]]; then
  printf "Report is missing: %s\n" "${REPORT}" >&2
  printf "Generate the report first, then run this script again.\n" >&2
  exit 1
fi

scripts/start_ai_cluster_server.sh --start

cleanup() {
  scripts/start_ai_cluster_server.sh --stop >/dev/null 2>&1 || true
}
trap cleanup EXIT INT TERM

sleep 1
printf "AI service started on http://127.0.0.1:8787\n"
printf "Report will be available at http://127.0.0.1:8000/Soil_1.html\n"

python3 -m http.server 8000 --directory results/Soil_1/report
