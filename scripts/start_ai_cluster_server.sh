#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

LOCAL_ENV="${ROOT_DIR}/config/local_ai.env"
RUNTIME_DIR="${ROOT_DIR}/.runtime"
PID_FILE="${RUNTIME_DIR}/ai_cluster_server.pid"
LOG_FILE="${RUNTIME_DIR}/ai_cluster_server.log"
HOST="${AI_HOST:-127.0.0.1}"
PORT="${AI_PORT:-8787}"

mkdir -p "${RUNTIME_DIR}"

load_env() {
  if [[ -f "${LOCAL_ENV}" ]]; then
    set -a
    # shellcheck source=/dev/null
    source "${LOCAL_ENV}"
    set +a
  fi
}

usage() {
  cat <<EOF
Usage: $(basename "$0") [--start|--stop|--restart] [--model MODEL_ID]

  --start    Start the AI cluster server in the background
  --stop     Stop the background AI cluster server
  --restart  Restart the background AI cluster server
  --model    Override the NVIDIA model for this start/restart

If no option is provided, --start is used.
EOF
}

pid_from_file() {
  if [[ -f "${PID_FILE}" ]]; then
    tr -d '[:space:]' < "${PID_FILE}"
  fi
}

is_running() {
  local pid="${1:-}"
  [[ -n "${pid}" ]] || return 1
  kill -0 "${pid}" >/dev/null 2>&1
}

cleanup_stale_pid() {
  local pid
  pid="$(pid_from_file || true)"
  if [[ -n "${pid}" ]] && ! is_running "${pid}"; then
    rm -f "${PID_FILE}"
  fi
}

require_key() {
  load_env
  if [[ -z "${NVIDIA_API_KEY:-}" ]]; then
    printf "NVIDIA API key: "
    read -r -s NVIDIA_API_KEY
    printf "\n"
    export NVIDIA_API_KEY
  fi
  export NVIDIA_MODEL="${CLI_MODEL:-${NVIDIA_MODEL:-deepseek-ai/deepseek-v4-pro}}"
}

start_server() {
  cleanup_stale_pid
  local pid
  pid="$(pid_from_file || true)"
  if is_running "${pid}"; then
    printf "AI cluster server is already running (PID %s)\n" "${pid}"
    if [[ -n "${CLI_MODEL:-}" ]]; then
      printf "Requested model override was ignored because the running service was kept alive.\n"
    fi
    printf "Log: %s\n" "${LOG_FILE}"
    return 0
  fi

  require_key

  nohup python3 scripts/ai_cluster_server.py \
    --host "${HOST}" \
    --port "${PORT}" \
    --results-dir results \
    >> "${LOG_FILE}" 2>&1 < /dev/null &

  pid=$!
  echo "${pid}" > "${PID_FILE}"

  sleep 1
  if is_running "${pid}"; then
    printf "AI cluster server started in background.\n"
    printf "PID: %s\n" "${pid}"
    printf "Model: %s\n" "${NVIDIA_MODEL}"
    printf "Endpoint: http://%s:%s/analyze_cluster\n" "${HOST}" "${PORT}"
    printf "Log: %s\n" "${LOG_FILE}"
    return 0
  fi

  printf "AI cluster server failed to start. Last log lines:\n" >&2
  tail -n 20 "${LOG_FILE}" >&2 || true
  rm -f "${PID_FILE}"
  return 1
}

stop_server() {
  cleanup_stale_pid
  local pid
  pid="$(pid_from_file || true)"
  if [[ -z "${pid}" ]]; then
    printf "AI cluster server is not running.\n"
    return 0
  fi

  if ! is_running "${pid}"; then
    rm -f "${PID_FILE}"
    printf "AI cluster server is not running.\n"
    return 0
  fi

  kill "${pid}" >/dev/null 2>&1 || true
  for _ in {1..20}; do
    if ! is_running "${pid}"; then
      rm -f "${PID_FILE}"
      printf "AI cluster server stopped.\n"
      return 0
    fi
    sleep 0.25
  done

  kill -9 "${pid}" >/dev/null 2>&1 || true
  rm -f "${PID_FILE}"
  printf "AI cluster server force-stopped.\n"
}

restart_server() {
  stop_server
  start_server
}

ACTION="--start"
CLI_MODEL=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --start|--stop|--restart)
      ACTION="$1"
      shift
      ;;
    --model)
      if [[ $# -lt 2 ]]; then
        printf "Missing value for --model\n\n" >&2
        usage >&2
        exit 1
      fi
      CLI_MODEL="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      printf "Unknown option: %s\n\n" "$1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

case "${ACTION}" in
  --start)
    start_server
    ;;
  --stop)
    stop_server
    ;;
  --restart)
    restart_server
    ;;
  -h|--help)
    usage
    ;;
  *)
    printf "Unknown option: %s\n\n" "${ACTION}" >&2
    usage >&2
    exit 1
    ;;
esac
