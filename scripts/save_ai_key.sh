#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

mkdir -p config

printf "NVIDIA API key: "
read -r -s NVIDIA_API_KEY
printf "\n"

MODEL="${NVIDIA_MODEL:-deepseek-ai/deepseek-v4-pro}"

umask 077
{
  printf "NVIDIA_API_KEY=%q\n" "${NVIDIA_API_KEY}"
  printf "NVIDIA_MODEL=%q\n" "${MODEL}"
} > config/local_ai.env

chmod 600 config/local_ai.env
printf "Saved local AI credentials to config/local_ai.env (ignored by git).\n"
