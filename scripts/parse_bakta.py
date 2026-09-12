import json

import pandas as pd

from common import write_json, write_tsv


sample = snakemake.params.sample

features = pd.read_csv(
    snakemake.input.tsv,
    sep="\t",
    comment="#",
    header=None,
    names=[
        "contig",
        "type",
        "start",
        "end",
        "strand",
        "locus_tag",
        "gene",
        "product",
        "dbxrefs",
    ],
)
features["sample"] = sample
write_tsv(features, snakemake.output.features)

metadata = {
    "sample": sample,
    "feature_count": int(len(features)),
    "cds_count": int((features.get("type", pd.Series(dtype=str)) == "cds").sum()) if "type" in features.columns else None,
    "input_files": dict(snakemake.input.items()),
}

with open(snakemake.input.json, "r", encoding="utf-8") as handle:
    bakta_json = json.load(handle)

if isinstance(bakta_json, dict):
    metadata["bakta"] = {
        "genome": bakta_json.get("genome"),
        "stats": bakta_json.get("stats"),
        "version": bakta_json.get("version"),
    }

write_json(metadata, snakemake.output.metadata)
