#!/bin/bash
set -e

WORK_DIR="/work"
APP_DIR="/app"
PYTHON="/opt/conda/bin/python"
LOCAL_AI_ENV="$WORK_DIR/config/local_ai.env"

mkdir -p "$WORK_DIR/data/fasta" "$WORK_DIR/data/bakta" "$WORK_DIR/results"

if [ ! -f "$WORK_DIR/config/config.yaml" ]; then
    cp -r "$APP_DIR/config/." "$WORK_DIR/config/"
fi

echo "Checking bundled tool runtimes..."
$PYTHON "$APP_DIR/scripts/startup_checks.py"

BAKTA_DB_TYPE="${BAKTA_DB_TYPE:-light}"
BAKTA_CHECK=("$PYTHON" "$APP_DIR/scripts/check_bakta_database.py" --db-root /db --type "$BAKTA_DB_TYPE")
if ! "${BAKTA_CHECK[@]}"; then
    if [ "${AUTO_PREPARE_DATABASES:-true}" = "true" ]; then
        BGC_DB_ROOT=/db BGC_PYTHON="$PYTHON" BAKTA_DB_TYPE="$BAKTA_DB_TYPE" \
            bash "$APP_DIR/scripts/fetch_bakta_db.sh"
    else
        echo "The selected Bakta database is incomplete and AUTO_PREPARE_DATABASES is disabled." >&2
        exit 1
    fi
fi
export BAKTA_DB_TYPE
if [ "$BAKTA_DB_TYPE" = "light" ]; then
    export BAKTA_DB="/db/bakta/db-light"
else
    export BAKTA_DB="/db/bakta/db"
fi

DB_CHECK=("$PYTHON" "$APP_DIR/scripts/check_databases.py" --manifest "$APP_DIR/db/manifest.yaml" --db-root /db)
if "${DB_CHECK[@]}"; then
    echo "Required databases are already available; database preparation skipped."
elif [ "${AUTO_PREPARE_DATABASES:-true}" = "true" ]; then
    echo "Required databases are incomplete; downloading supported resources..."
    BGC_DB_ROOT=/db BGC_PYTHON="$PYTHON" ARTS_REFERENCE="${ARTS_REFERENCE:-actinobacteria}" \
        bash "$APP_DIR/scripts/fetch_databases.sh"
    "${DB_CHECK[@]}"
else
    echo "Required databases are incomplete and AUTO_PREPARE_DATABASES is disabled." >&2
    exit 1
fi

if [ -z "${NVIDIA_API_KEY:-}" ] && [ -f "$LOCAL_AI_ENV" ]; then
    echo "Loading local AI credentials from $LOCAL_AI_ENV"
    set -a
    # shellcheck source=/dev/null
    . "$LOCAL_AI_ENV"
    set +a
fi

echo "Starting AI cluster server on port 8484..."
export AI_TIMEOUT="${AI_TIMEOUT:-240}"
export AI_MAX_TOKENS="${AI_MAX_TOKENS:-1800}"
$PYTHON "$APP_DIR/scripts/ai_cluster_server.py" \
    --host 0.0.0.0 \
    --port 8484 \
    --results-dir "$WORK_DIR/results" &

AI_PID=$!
echo "AI server PID: $AI_PID"

echo "Starting NiceGUI on port 8778..."
exec $PYTHON "$APP_DIR/scripts/web_ui.py" \
    --host 0.0.0.0 \
    --port 8778 \
    --work-dir "$WORK_DIR" \
    --app-dir "$APP_DIR"
