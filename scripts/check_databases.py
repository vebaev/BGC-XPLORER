from __future__ import print_function

import os
import argparse

import yaml


def resolve_database_path(host_path, database_root):
    if not host_path:
        return ""
    relative = host_path
    if relative == "db":
        relative = ""
    elif relative.startswith("db/"):
        relative = relative[3:]
    return os.path.join(str(database_root), relative)


def inspect_databases(manifest, database_root):
    status = {}
    for tool_name in sorted(manifest):
        entry = manifest[tool_name] or {}
        required = bool(entry.get("required", False))
        path = resolve_database_path(entry.get("host_path", ""), database_root)
        required_files = entry.get("required_files", []) or []
        exists = bool(path) and os.path.isdir(path)
        file_count = count_real_files(path) if exists else 0
        missing_files = [name for name in required_files if not os.path.exists(os.path.join(path, name))]
        ready = (not required) or (exists and file_count > 0 and not missing_files)
        status[tool_name] = {
            "required": required,
            "ready": ready,
            "path": path,
            "exists": exists,
            "file_count": file_count,
            "missing_files": missing_files,
            "manifest_status": entry.get("status", "unknown"),
            "notes": entry.get("notes", ""),
        }
    return status


def missing_required_names(database_status):
    return [name for name in sorted(database_status) if database_status[name]["required"] and not database_status[name]["ready"]]


def count_real_files(path):
    total = 0
    for _, _, files in os.walk(path):
        for name in files:
            if not name.startswith("."):
                total += 1
    return total


def main():
    parser = argparse.ArgumentParser(description="Validate BGC-XPLORER reference databases")
    parser.add_argument("--manifest", default=os.path.join("db", "manifest.yaml"))
    parser.add_argument("--db-root", default="db")
    parser.add_argument("--missing-only", action="store_true")
    args = parser.parse_args()
    manifest_path = args.manifest
    if not os.path.exists(manifest_path):
        raise SystemExit("Missing manifest: {0}".format(manifest_path))

    with open(manifest_path, "r") as handle:
        manifest = yaml.safe_load(handle) or {}

    print("Database status")
    print("===============")

    missing_required = []

    database_status = inspect_databases(manifest, args.db_root)
    if args.missing_only:
        for name in missing_required_names(database_status):
            print(name)
        return
    for tool_name in sorted(database_status):
        entry = database_status[tool_name]
        required = entry["required"]
        ready = entry["ready"]

        if required and not ready:
            missing_required.append(tool_name)

        print("{tool}: required={required} ready={ready} manifest_status={status} exists={exists} files={files}".format(
            tool=tool_name,
            required=str(required).lower(),
            ready=str(ready).lower(),
            status=entry["manifest_status"],
            exists=str(entry["exists"]).lower(),
            files=entry["file_count"],
        ))
        if entry["path"]:
            print("  path: {0}".format(entry["path"]))
        if entry["notes"]:
            print("  notes: {0}".format(entry["notes"]))
        if entry["missing_files"]:
            print("  missing_required_files: {0}".format(", ".join(entry["missing_files"])))

    if missing_required:
        raise SystemExit(
            "Missing required database resources for: {0}".format(", ".join(missing_required))
        )


if __name__ == "__main__":
    main()
