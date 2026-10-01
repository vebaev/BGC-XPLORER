#!/usr/bin/env bash
# Run the BGC-XPLORER workflow for the benchmark genomes.
#
# The reported benchmark ran commit 81d608fe6571d1559a2e5ad609c44dd36b5fadc5
# mounted over the v1.1.1 image; tool versions are identical because the
# Dockerfile did not change between them. To repeat it with a later release,
# set IMAGE to that release and MOUNT_CODE=0.
#
# Databases: DB_ROOT holds antiSMASH, DeepBGC and ARTS; Bakta, eggNOG and dbCAN
# can live elsewhere and are mounted over it. PATH and BAKTA_DB mirror what
# scripts/entrypoint.sh exports.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
R="$(cd "$HERE/.." && pwd)"
B="${WORKDIR:-$HERE/work}"
NAME=${NAME:-bgc-benchmark}
IMAGE=${IMAGE:-ghcr.io/vebaev/bgc-xplorer:1.1.1}
COMMIT=${COMMIT:-81d608fe6571d1559a2e5ad609c44dd36b5fadc5}
DB_ROOT=${DB_ROOT:-$R/db}
BAKTA_DB_DIR=${BAKTA_DB_DIR:-$DB_ROOT/bakta}
EGGNOG_DB_DIR=${EGGNOG_DB_DIR:-$DB_ROOT/eggnog}
DBCAN_DB_DIR=${DBCAN_DB_DIR:-$DB_ROOT/dbcan}
BAKTA_DB_TYPE=${BAKTA_DB_TYPE:-light}
if [ "$BAKTA_DB_TYPE" = "full" ]; then BAKTA_SUBDIR=db; else BAKTA_SUBDIR=db-light; fi

code_mounts=()
if [ "${MOUNT_CODE:-1}" = "1" ]; then
  code_mounts=(-v "$R/Snakefile:/app/Snakefile:ro" -v "$R/rules:/app/rules:ro"
               -v "$R/scripts:/app/scripts:ro" -v "$R/envs:/app/envs:ro"
               -e "BGC_XPLORER_VERSION=${BGC_XPLORER_VERSION:-$COMMIT}" -e "BGC_XPLORER_COMMIT=$COMMIT")
fi

docker run --rm --name "$NAME" --entrypoint "" \
  -v "$B/config:/work/config" -v "$B/data:/work/data" -v "$B/results:/work/results" \
  -v "$DB_ROOT:/db" \
  -v "$BAKTA_DB_DIR:/db/bakta" \
  -v "$EGGNOG_DB_DIR:/db/eggnog" \
  -v "$DBCAN_DB_DIR:/db/dbcan" \
  "${code_mounts[@]}" \
  -e PATH=/opt/conda/envs/bakta/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin -e BAKTA_DB=/db/bakta/$BAKTA_SUBDIR \
  -e BGC_THREADS=${BGC_THREADS:-16} -e BAKTA_DB_TYPE=$BAKTA_DB_TYPE -e ARTS_REFERENCE=actinobacteria \
  "$IMAGE" \
  /opt/conda/bin/snakemake -s /app/Snakefile --configfile /work/config/config.yaml \
    --directory /work --rerun-incomplete --printshellcmds "$@"
