#!/usr/bin/env python3
from __future__ import print_function

import argparse


PRIVATE_FILES = {".env", "config/local_ai.env"}
PRIVATE_PREFIXES = (
    ".cache/", ".conda-pkgs/", ".local/", ".runtime/", ".snakemake/",
    ".snakemake-conda/", "backup/", "data/", "results/",
)


def is_private_path(path):
    normalized = path.replace("\\", "/")
    while normalized.startswith("./"):
        normalized = normalized[2:]
    if normalized in PRIVATE_FILES or normalized.startswith(PRIVATE_PREFIXES):
        return True
    parts = normalized.split("/")
    if ".git" in parts or "__pycache__" in parts:
        return True
    if normalized.startswith("db/") and normalized != "db/manifest.yaml":
        return True
    return normalized.endswith((".pyc", ".pem", ".key"))


def find_public_tree_violations(paths):
    return [path for path in paths if is_private_path(path)]


def main():
    parser = argparse.ArgumentParser(description="Reject private/generated files from the public source tree")
    parser.add_argument("paths", nargs="*")
    args = parser.parse_args()
    violations = find_public_tree_violations(args.paths)
    for path in violations:
        print("Public tree violation: {0}".format(path))
    if violations:
        raise SystemExit(1)
    print("Public source tree contains no blocked paths.")


if __name__ == "__main__":
    main()
