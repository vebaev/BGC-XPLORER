#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DB_DIR="${ROOT_DIR}/db"

mkdir -p "${DB_DIR}/antismash"
mkdir -p "${DB_DIR}/deepbgc"
mkdir -p "${DB_DIR}/arts/actinobacteria"
mkdir -p "${DB_DIR}/eggnog"

cat <<'EOF'
Database scaffold prepared under:
  db/antismash
  db/deepbgc
  db/arts/actinobacteria
  db/eggnog

Next steps:
1. antiSMASH
   Populate db/antismash with antiSMASH databases.
   Typical standalone/local flow uses the antiSMASH database downloader.

2. DeepBGC
   Download DeepBGC models and Pfam resources into db/deepbgc.
   Typical flow uses: deepbgc download

3. ARTS
   Copy or unpack the selected ARTS reference set into:
   db/arts/actinobacteria

4. eggNOG-mapper
   Download eggNOG annotation and DIAMOND databases into db/eggnog.

5. Validate
   python3 scripts/check_databases.py
EOF
