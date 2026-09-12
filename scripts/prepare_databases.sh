#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DB_DIR="${BGC_DB_ROOT:-${ROOT_DIR}/db}"
ARTS_REFERENCE="${ARTS_REFERENCE:-actinobacteria}"

mkdir -p "${DB_DIR}/antismash"
mkdir -p "${DB_DIR}/deepbgc"
mkdir -p "${DB_DIR}/arts/${ARTS_REFERENCE}"
mkdir -p "${DB_DIR}/eggnog"

cat <<EOF
Database scaffold prepared under:
  ${DB_DIR}/antismash
  ${DB_DIR}/deepbgc
  ${DB_DIR}/arts/${ARTS_REFERENCE}
  ${DB_DIR}/eggnog

Next steps:
1. antiSMASH
   Populate db/antismash with antiSMASH databases.
   Typical standalone/local flow uses the antiSMASH database downloader.

2. DeepBGC
   Download DeepBGC models and Pfam resources into db/deepbgc.
   Typical flow uses: deepbgc download

3. ARTS
   The default Actinobacteria reference is installed automatically by:
   bash scripts/fetch_databases.sh

4. eggNOG-mapper
   Download eggNOG annotation and DIAMOND databases into db/eggnog.

5. Validate
   python3 scripts/check_databases.py
EOF
