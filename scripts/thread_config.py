"""Validate the shared Docker thread limit used by the workflow and tools."""

import os


def configured_threads(environment=None, default=4):
    environment = os.environ if environment is None else environment
    raw = environment.get("BGC_THREADS")
    if raw is None:
        return default
    if not raw or raw != raw.strip() or not raw.isdigit() or int(raw) < 1:
        raise ValueError("BGC_THREADS must be a positive integer")
    return int(raw)
