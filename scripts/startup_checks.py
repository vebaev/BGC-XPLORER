#!/usr/bin/env python3
from __future__ import print_function

import argparse
import os


DEFAULT_EXECUTABLES = [
    "/opt/conda/bin/python",
    "/opt/conda/bin/snakemake",
    "/opt/conda/envs/antismash/bin/antismash",
    "/opt/conda/envs/gecco/bin/gecco",
    "/opt/conda/envs/deepbgc/bin/deepbgc",
    "/opt/conda/envs/eggnog/bin/emapper.py",
    "/opt/conda/envs/arts/bin/python",
]


def missing_executables(paths):
    return [path for path in paths if not (os.path.isfile(path) and os.access(path, os.X_OK))]


def main():
    parser = argparse.ArgumentParser(description="Check BGC-XPLORER container tools")
    parser.add_argument("paths", nargs="*", default=DEFAULT_EXECUTABLES)
    args = parser.parse_args()
    missing = missing_executables(args.paths)
    if missing:
        for path in missing:
            print("Missing executable: {0}".format(path))
        raise SystemExit(1)
    print("All required tool runtimes are available; installation skipped.")


if __name__ == "__main__":
    main()
