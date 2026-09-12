#!/usr/bin/env python3

import os
import shutil
from pathlib import Path

import yaml


TOOLS = [
    {
        "name": "antismash",
        "label": "antiSMASH",
        "config_key": "executable",
        "path_checks": ["results/{sample}/antismash/index.html", "results/{sample}/antismash/regions.js"],
        "path_mode": "all",
        "fallback_cmd": "antismash",
    },
    {
        "name": "gecco",
        "label": "GECCO",
        "config_key": "executable",
        "path_checks": ["results/{sample}/gecco/{sample}.clusters.tsv"],
        "path_mode": "any",
        "fallback_cmd": "gecco",
    },
    {
        "name": "deepbgc",
        "label": "DeepBGC",
        "config_key": "executable",
        "path_checks": ["results/{sample}/deepbgc/deepbgc.bgc.tsv"],
        "path_mode": "any",
        "fallback_cmd": "deepbgc",
    },
    {
        "name": "arts",
        "label": "ARTS",
        "config_key": "python",
        "path_checks": ["results/{sample}/arts/tables/knownhits.tsv", "results/{sample}/arts/tables/coretable.tsv"],
        "path_mode": "all",
        "fallback_cmd": "python3",
        "script_key": "script",
    },
    {
        "name": "eggnog",
        "label": "eggNOG-mapper",
        "config_key": "executable",
        "path_checks": ["results/{sample}/eggnog/{sample}.emapper.annotations"],
        "path_mode": "any",
        "fallback_cmd": "emapper.py",
    },
    {
        "name": "dbcan",
        "label": "dbCAN",
        "config_key": "binary",
        "path_checks": [
            "results/{sample}/dbcan/substrate_prediction.tsv",
            "results/{sample}/dbcan/cgc_standard_out_summary.tsv",
            "results/{sample}/dbcan/overview.tsv",
        ],
        "path_mode": "any",
        "fallback_cmd": "run_dbcan",
    },
]


def load_config():
    with open("config/config.yaml", "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def check_existing_outputs(sample, patterns, mode):
    resolved = [pattern.format(sample=sample) for pattern in patterns]
    existing = [path for path in resolved if Path(path).exists()]
    if mode == "all":
        return len(existing) == len(resolved), existing
    return bool(existing), existing


def resolve_command(configured, fallback_cmd):
    configured = (configured or "").strip()
    if configured:
        if Path(configured).exists() and os.access(configured, os.X_OK):
            return configured, "configured"
        return configured, "configured_missing"
    fallback = shutil.which(fallback_cmd)
    if fallback:
        return fallback, "path"
    return "", "missing"


def main():
    config = load_config()
    sample_rows = Path("config/samples.tsv").read_text(encoding="utf-8").strip().splitlines()
    sample = sample_rows[1].split("\t")[0] if len(sample_rows) > 1 else ""
    reuse_existing = bool(config.get("execution", {}).get("reuse_existing_outputs", True))

    print("Tool runtime status")
    print("===================")
    print("sample: {0}".format(sample or "<none>"))
    print("reuse_existing_outputs: {0}".format(str(reuse_existing).lower()))

    for tool in TOOLS:
        tool_cfg = config.get("tools", {}).get(tool["name"], {})
        runtime_value = tool_cfg.get(tool["config_key"], "")
        command, source = resolve_command(runtime_value, tool["fallback_cmd"])
        has_outputs, existing_outputs = check_existing_outputs(sample, tool["path_checks"], tool["path_mode"]) if sample else (False, [])

        status_bits = []
        if command and source in {"configured", "path"}:
            status_bits.append("runtime=ready")
        elif source == "configured_missing":
            status_bits.append("runtime=configured_missing")
        else:
            status_bits.append("runtime=missing")

        status_bits.append("cached_outputs={0}".format("ready" if has_outputs else "missing"))

        if command and source in {"configured", "path"}:
            recommended = "can_run_locally"
        elif reuse_existing and has_outputs:
            recommended = "reuse_existing_outputs"
        else:
            recommended = "needs_runtime_or_docker"

        print("{label}: {bits} recommended={recommended}".format(
            label=tool["label"],
            bits=" ".join(status_bits),
            recommended=recommended,
        ))
        if command:
            print("  command: {0}".format(command))
            print("  command_source: {0}".format(source))
        elif runtime_value:
            print("  configured_command: {0}".format(runtime_value))
            print("  command_source: {0}".format(source))
        if "script_key" in tool:
            script_path = tool_cfg.get(tool["script_key"], "")
            print("  script: {0}".format(script_path or "<unset>"))
        if existing_outputs:
            print("  existing_outputs:")
            for path in existing_outputs:
                print("    - {0}".format(path))


if __name__ == "__main__":
    main()
