from pathlib import Path

from common import write_json


required = {
    "gbff": snakemake.input.gbff,
    "fna": snakemake.input.fna,
    "faa": snakemake.input.faa,
    "gff3": snakemake.input.gff3,
    "json": snakemake.input.json,
    "tsv": snakemake.input.tsv,
}

status = {"sample": snakemake.params.sample, "files": {}, "ok": True}

for key, path in required.items():
    exists = Path(path).exists()
    size = Path(path).stat().st_size if exists else 0
    status["files"][key] = {"path": path, "exists": exists, "size": size}
    status["ok"] = status["ok"] and exists and size > 0

if not status["ok"]:
    missing = [name for name, meta in status["files"].items() if not meta["exists"] or meta["size"] == 0]
    raise ValueError(
        "Missing or empty Bakta inputs for {sample}: {missing}".format(
            sample=snakemake.params.sample,
            missing=", ".join(missing),
        )
    )

write_json(status, snakemake.output[0])
