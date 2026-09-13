/bin/bash: warning: setlocale: LC_ALL: cannot change locale (C.UTF-8)
#!/usr/bin/env python3
from __future__ import print_function

import argparse
import os

try:
    from .thread_config import configured_threads
except ImportError:
    from thread_config import configured_threads


DEFAULT_EXECUTABLES = [
    "/opt/conda/bin/python",
    "/opt/conda/bin/snakemake",
    "/opt/conda/envs/antismash/bin/antismash",
    "/opt/conda/envs/gecco/bin/gecco",
    "/opt/conda/envs/deepbgc/bin/deepbgc",
    "/opt/conda/envs/eggnog/bin/emapper.py",
    "/opt/conda/envs/arts/bin/python",
    "/opt/conda/envs/bakta/bin/bakta",
    "/opt/conda/envs/bakta/bin/bakta_db",
]


def missing_executables(paths):
    return [path for path in paths if not (os.path.isfile(path) and os.access(path, os.X_OK))]


def main():
    parser = argparse.ArgumentParser(description="Check BGC-XPLORER container tools")
    parser.add_argument("paths", nargs="*", default=DEFAULT_EXECUTABLES)
    args = parser.parse_args()
    missing = missing_executables(args.paths)
    try:
        thread_count = configured_threads()
    except ValueError as error:
        print(str(error))
        raise SystemExit(2)
    if missing:
        for path in missing:
            print("Missing executable: {0}".format(path))
        raise SystemExit(1)
    print("All required tool runtimes are available; installation skipped.")
    print("Workflow thread limit: {0}".format(thread_count))


if __name__ == "__main__":
    main()
