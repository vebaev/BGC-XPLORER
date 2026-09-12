/bin/bash: warning: setlocale: LC_ALL: cannot change locale (C.UTF-8)
#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DB_ROOT="${BGC_DB_ROOT:-${ROOT_DIR}/db}"
REFERENCE="${ARTS_REFERENCE:-actinobacteria}"
ARTS_COMMIT="${ARTS_COMMIT:-8922f296b2a532ba51f4d5daa6a807838c21be24}"
ARTS_ROOT="${DB_ROOT}/arts"
READY_FILE="${ARTS_ROOT}/.bgc-xplorer-${REFERENCE}.ready"

if [ "${REFERENCE}" != "actinobacteria" ]; then
  echo "Unsupported ARTS reference: ${REFERENCE}. Supported value: actinobacteria" >&2
  exit 2
fi

required_files=(
  "${ARTS_ROOT}/${REFERENCE}/coremodels.hmm"
  "${ARTS_ROOT}/${REFERENCE}/model_metadata.json"
  "${ARTS_ROOT}/${REFERENCE}/genematrix.txt"
  "${ARTS_ROOT}/knownresistance.hmm"
  "${ARTS_ROOT}/dufmodels.hmm"
  "${READY_FILE}"
)

is_ready() {
  local path
  for path in "${required_files[@]}"; do
    [ -s "${path}" ] || return 1
  done
}

if is_ready; then
  echo "ARTS ${REFERENCE}: already valid; download skipped."
  exit 0
fi

mkdir -p "${ARTS_ROOT}" "${DB_ROOT}/.locks"
LOCK_DIR="${DB_ROOT}/.locks/arts-${REFERENCE}.lock"
waited=0
until mkdir "${LOCK_DIR}" 2>/dev/null; do
  if is_ready; then
    echo "ARTS ${REFERENCE}: already valid; download skipped."
    exit 0
  fi
  if [ "${waited}" -ge 600 ]; then
    echo "Timed out waiting for ARTS database preparation lock: ${LOCK_DIR}" >&2
    exit 1
  fi
  sleep 2
  waited=$((waited + 2))
done

STAGE_DIR="$(mktemp -d "${DB_ROOT}/.arts-staging.XXXXXX")"
SOURCE_TEMP=""
cleanup() {
  rm -rf "${STAGE_DIR}"
  if [ -n "${SOURCE_TEMP}" ]; then
    rm -rf "${SOURCE_TEMP}"
  fi
  rmdir "${LOCK_DIR}" 2>/dev/null || true
}
trap cleanup EXIT INT TERM

if is_ready; then
  echo "ARTS ${REFERENCE}: already valid; download skipped."
  exit 0
fi

if [ -n "${ARTS_SOURCE_DIR:-}" ]; then
  SOURCE_DIR="${ARTS_SOURCE_DIR}"
elif [ -f "/opt/arts/reference/${REFERENCE}.zip" ] && [ -f "/opt/arts/reference/hmm_models.zip" ]; then
  SOURCE_DIR="/opt/arts/reference"
elif [ -f "${ROOT_DIR}/tools/arts/reference/${REFERENCE}.zip" ] && [ -f "${ROOT_DIR}/tools/arts/reference/hmm_models.zip" ]; then
  SOURCE_DIR="${ROOT_DIR}/tools/arts/reference"
else
  SOURCE_TEMP="$(mktemp -d "${DB_ROOT}/.arts-source.XXXXXX")"
  echo "Downloading pinned ARTS source ${ARTS_COMMIT}..."
  curl --fail --location --retry 3 --output "${SOURCE_TEMP}/arts.tar.gz" \
    "https://github.com/ZiemertLab/ARTS/archive/${ARTS_COMMIT}.tar.gz"
  mkdir -p "${SOURCE_TEMP}/arts"
  tar xzf "${SOURCE_TEMP}/arts.tar.gz" --strip-components=1 -C "${SOURCE_TEMP}/arts"
  SOURCE_DIR="${SOURCE_TEMP}/arts/reference"
fi

for archive in "${SOURCE_DIR}/${REFERENCE}.zip" "${SOURCE_DIR}/hmm_models.zip"; do
  if [ ! -s "${archive}" ]; then
    echo "Missing ARTS archive: ${archive}" >&2
    exit 1
  fi
done

echo "Installing ARTS ${REFERENCE} reference under ${ARTS_ROOT}..."
python_bin="${BGC_PYTHON:-$(command -v python3)}"
"${python_bin}" - "${SOURCE_DIR}/${REFERENCE}.zip" "${SOURCE_DIR}/hmm_models.zip" "${STAGE_DIR}" <<'PY'
import os
import sys
import zipfile

reference_archive, models_archive, destination = sys.argv[1:]
for archive_path in (reference_archive, models_archive):
    with zipfile.ZipFile(archive_path) as archive:
        for member in archive.infolist():
            target = os.path.realpath(os.path.join(destination, member.filename))
            if not target.startswith(os.path.realpath(destination) + os.sep):
                raise SystemExit("Unsafe path in ARTS archive: {0}".format(member.filename))
        archive.extractall(destination)
PY

for path in \
  "${STAGE_DIR}/${REFERENCE}/coremodels.hmm" \
  "${STAGE_DIR}/${REFERENCE}/model_metadata.json" \
  "${STAGE_DIR}/${REFERENCE}/genematrix.txt" \
  "${STAGE_DIR}/knownresistance.hmm" \
  "${STAGE_DIR}/dufmodels.hmm"; do
  if [ ! -s "${path}" ]; then
    echo "ARTS archive is incomplete; missing ${path#${STAGE_DIR}/}" >&2
    exit 1
  fi
done

rm -rf "${ARTS_ROOT:?}/${REFERENCE}"
mv "${STAGE_DIR}/${REFERENCE}" "${ARTS_ROOT}/${REFERENCE}"
mv "${STAGE_DIR}/knownresistance.hmm" "${ARTS_ROOT}/knownresistance.hmm"
mv "${STAGE_DIR}/dufmodels.hmm" "${ARTS_ROOT}/dufmodels.hmm"
if [ -f "${STAGE_DIR}/barnap_bact_rRna.hmm" ]; then
  mv "${STAGE_DIR}/barnap_bact_rRna.hmm" "${ARTS_ROOT}/barnap_bact_rRna.hmm"
fi
printf 'reference=%s\narts_commit=%s\n' "${REFERENCE}" "${ARTS_COMMIT}" > "${READY_FILE}"

echo "ARTS ${REFERENCE}: installed successfully."
