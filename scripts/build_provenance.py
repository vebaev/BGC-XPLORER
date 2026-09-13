import hashlib
import json
import os
from datetime import datetime

import yaml

try:
    import fcntl
except ImportError:  # pragma: no cover - non-POSIX development hosts
    fcntl = None


CHUNK_SIZE = 1024 * 1024
SECRET_KEY_PARTS = ("api_key", "apikey", "password", "secret", "token")


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(CHUNK_SIZE)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


class DatabaseChecksumCache(object):
    """Cache immutable database hashes by path, size, and nanosecond mtime."""

    def __init__(self, database_root):
        self.database_root = os.path.abspath(database_root)
        self.cache_path = os.path.join(self.database_root, ".bgc-xplorer-checksums.json")
        self.lock_path = self.cache_path + ".lock"

    def checksum(self, path):
        stat = os.stat(path)
        key = os.path.abspath(path)
        signature = {
            "size": stat.st_size,
            "mtime_ns": getattr(stat, "st_mtime_ns", int(stat.st_mtime * 1000000000)),
        }
        os.makedirs(self.database_root, exist_ok=True)
        with open(self.lock_path, "a+") as lock_handle:
            if fcntl is not None:
                fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
            cache = self._read()
            cached = cache.get(key, {})
            if cached.get("size") == signature["size"] and cached.get("mtime_ns") == signature["mtime_ns"]:
                return cached.get("sha256")
            signature["sha256"] = sha256_file(path)
            cache[key] = signature
            self._write(cache)
            return signature["sha256"]

    def _read(self):
        try:
            with open(self.cache_path, "r") as handle:
                value = json.load(handle)
            return value if isinstance(value, dict) else {}
        except (IOError, ValueError):
            return {}

    def _write(self, cache):
        temporary = self.cache_path + ".tmp.{0}".format(os.getpid())
        with open(temporary, "w") as handle:
            json.dump(cache, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary, self.cache_path)


def file_record(path, checksum=None):
    exists = os.path.isfile(path)
    return {
        "path": str(path),
        "exists": exists,
        "size_bytes": os.path.getsize(path) if exists else None,
        "sha256": checksum(path) if exists and checksum else (sha256_file(path) if exists else None),
    }


def _database_path(host_path, database_root):
    relative = host_path or ""
    if relative == "db":
        relative = ""
    elif relative.startswith("db/"):
        relative = relative[3:]
    return os.path.join(database_root, relative)


def redact_secrets(value):
    if isinstance(value, dict):
        clean = {}
        for key, item in value.items():
            normalized = str(key).lower().replace("-", "_")
            clean[key] = "[REDACTED]" if any(part in normalized for part in SECRET_KEY_PARTS) else redact_secrets(item)
        return clean
    if isinstance(value, list):
        return [redact_secrets(item) for item in value]
    return value


def build_provenance(sample, input_paths, artifact_paths, config_path,
                     database_manifest_path, database_root, environment=None):
    environment = dict(environment or {})
    with open(config_path, "r") as handle:
        effective_config = yaml.safe_load(handle) or {}
    with open(database_manifest_path, "r") as handle:
        database_manifest = yaml.safe_load(handle) or {}

    try:
        tool_versions = json.loads(environment.get("BGC_TOOL_VERSIONS_JSON", "{}"))
    except ValueError:
        tool_versions = {"error": "invalid BGC_TOOL_VERSIONS_JSON"}

    checksum_cache = DatabaseChecksumCache(database_root)
    databases = {}
    for name in sorted(database_manifest):
        definition = database_manifest[name] or {}
        root = _database_path(definition.get("host_path", ""), database_root)
        files = []
        for relative in definition.get("required_files", []) or []:
            path = os.path.join(root, relative) if root else relative
            record = file_record(path, checksum=checksum_cache.checksum)
            record["relative_path"] = relative
            files.append(record)
        databases[name] = {
            "required": bool(definition.get("required", False)),
            "source": definition.get("source"),
            "status": definition.get("status"),
            "files": files,
        }

    arts_reference = environment.get("ARTS_REFERENCE") or (
        effective_config.get("tools", {}).get("arts", {}).get("reference_set", "actinobacteria")
    )
    return {
        "schema_version": 1,
        "generated_at_utc": datetime.utcnow().replace(microsecond=0).isoformat() + "Z",
        "sample": sample,
        "application": {
            "name": "BGC-XPLORER",
            "version": environment.get("BGC_XPLORER_VERSION", "development"),
            "git_commit": environment.get("BGC_XPLORER_COMMIT", "unknown"),
            "image_reference": environment.get("BGC_IMAGE_REFERENCE", "unknown"),
        },
        "ai": {
            "provider": "NVIDIA",
            "model": environment.get("NVIDIA_MODEL", "not-configured"),
        },
        "execution": {
            "threads": int(environment.get("BGC_THREADS", "4")),
        },
        "arts_reference": arts_reference,
        "tools": tool_versions,
        "configuration": {
            "path": config_path,
            "sha256": sha256_file(config_path),
            "effective": redact_secrets(effective_config),
        },
        "inputs": [file_record(path) for path in sorted(set(input_paths))],
        "artifacts": [file_record(path) for path in sorted(set(artifact_paths))],
        "databases": databases,
    }


def write_provenance(provenance, output_path):
    parent = os.path.dirname(output_path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    temporary = output_path + ".tmp.{0}".format(os.getpid())
    with open(temporary, "w") as handle:
        json.dump(provenance, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(temporary, output_path)


if "snakemake" in globals():
    provenance = build_provenance(
        sample=str(snakemake.wildcards.sample),
        input_paths=[str(path) for path in snakemake.input.raw],
        artifact_paths=[str(path) for path in snakemake.input.artifacts],
        config_path=str(snakemake.input.config),
        database_manifest_path=str(snakemake.input.database_manifest),
        database_root=str(snakemake.params.database_root),
        environment=os.environ,
    )
    write_provenance(provenance, str(snakemake.output[0]))
