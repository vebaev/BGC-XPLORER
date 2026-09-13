#!/usr/bin/env python3
"""Validate the externally mounted Bakta database selected at container start."""

import argparse
import json
from pathlib import Path


def database_directory(database_root, db_type):
    return Path(database_root) / "bakta" / ("db-light" if db_type == "light" else "db")


def inspect_bakta_database(database_root, db_type):
    if db_type not in {"light", "full"}:
        raise ValueError("BAKTA_DB_TYPE must be 'light' or 'full'")
    path = database_directory(database_root, db_type)
    required = ["version.json", "bakta.db", "pscc.dmnd" if db_type == "light" else "psc.dmnd"]
    missing = [name for name in required if not (path / name).is_file()]
    amr_dir = path / "amrfinderplus-db"
    if not amr_dir.is_dir() or not any(item.is_file() for item in amr_dir.rglob("*")):
        missing.append("amrfinderplus-db")
    version = {}
    if (path / "version.json").is_file():
        try:
            version = json.loads((path / "version.json").read_text())
        except (OSError, ValueError):
            missing.append("valid version.json")
    recorded_type = version.get("type")
    expected_type = db_type == "light"
    if version and recorded_type not in {expected_type, db_type}:
        missing.append("matching database type")
    return {"ready": path.is_dir() and not missing, "path": str(path), "missing": missing, "version": version}


def main():
    parser = argparse.ArgumentParser(description="Validate a Bakta light or full database")
    parser.add_argument("--db-root", default="/db")
    parser.add_argument("--type", choices=("light", "full"), required=True)
    args = parser.parse_args()
    status = inspect_bakta_database(args.db_root, args.type)
    if not status["ready"]:
        print("Bakta {0} database is incomplete at {1}: {2}".format(args.type, status["path"], ", ".join(status["missing"])))
        raise SystemExit(1)
    version = status["version"]
    print("Bakta {0} database is ready at {1} (v{2}.{3}); download skipped.".format(
        args.type, status["path"], version.get("major", "?"), version.get("minor", "?")))


if __name__ == "__main__":
    main()

