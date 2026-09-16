from pathlib import Path
import json

import pandas as pd

from common import BGC_COLUMNS, coerce_frame_columns, empty_bgc_frame, first_existing_column, first_existing_path, write_tsv
from core_gene_evidence import antismash_core_records, join_records, parse_records


sample = snakemake.wildcards.sample
outdir = Path(snakemake.input.done).parent
regions_tsv = first_existing_path(outdir, ["*regions*.tsv", "*region*.tsv", "**/*regions*.tsv"])
regions_json = first_existing_path(outdir, [
    "*regions*.json",
    "*region*.json",
    "*.json",
    "**/*regions*.json",
    "**/*region*.json",
])

if regions_tsv is not None:
    df = pd.read_csv(str(regions_tsv), sep="\t")
    out = pd.DataFrame({
        "sample": sample,
        "tool": "antismash",
        "contig": first_existing_column(df, ["contig", "sequence", "record_id", "seq_id"]),
        "start": first_existing_column(df, ["start", "region_start"]),
        "end": first_existing_column(df, ["end", "region_end"]),
        "strand": first_existing_column(df, ["strand"], "."),
        "bgc_id": first_existing_column(df, ["region", "region_number", "bgc_id"], ""),
        "bgc_type": first_existing_column(df, ["type", "product", "product_prediction"], "unknown"),
        "product": first_existing_column(df, ["product", "product_prediction", "type"], ""),
        "score": first_existing_column(df, ["score"], ""),
        "confidence": first_existing_column(df, ["confidence"], ""),
        "core_gene_records": "",
        "source_file": str(regions_tsv),
    })
    out = coerce_frame_columns(out, BGC_COLUMNS)
    out["bgc_id"] = out["bgc_id"].astype(object)
    blank_ids = out["bgc_id"].astype(str).str.strip() == ""
    out.loc[blank_ids, "bgc_id"] = ["antismash_{0}".format(idx + 1) for idx in out.index[blank_ids]]
elif regions_json is not None:
    with open(str(regions_json), "r", encoding="utf-8") as handle:
        payload = json.load(handle)

    records = []
    areas = payload.get("records", []) if isinstance(payload, dict) else []
    for entry in areas:
        structured_core = antismash_core_records(entry.get("features", []))
        regions = entry.get("regions", [])
        if not regions:
            regions = entry.get("areas", [])
        for region_number, region in enumerate(regions, start=1):
            region_start = region.get("start", "")
            region_end = region.get("end", "")
            local_core = []
            for record in structured_core:
                parsed = parse_records(record)
                if not parsed:
                    continue
                gene = parsed[0]
                if (
                    gene["start"] is not None and gene["end"] is not None
                    and region_start != "" and region_end != ""
                    and gene["start"] <= int(region_end)
                    and gene["end"] >= int(region_start)
                ):
                    local_core.append(record)
            products = region.get("products", [])
            if not products and isinstance(region.get("protoclusters"), dict):
                products = [
                    cluster.get("product", "")
                    for cluster in region.get("protoclusters", {}).values()
                    if cluster.get("product", "")
                ]
            records.append({
                "sample": sample,
                "tool": "antismash",
                "contig": entry.get("id", ""),
                "start": region.get("start", ""),
                "end": region.get("end", ""),
                "strand": region.get("strand", "."),
                "bgc_id": region.get("idx", region.get("region_number", region_number)),
                "bgc_type": ",".join(products) if isinstance(products, list) else products or "unknown",
                "product": ",".join(products) if isinstance(products, list) else products,
                "score": "",
                "confidence": region.get("confidence", ""),
                "core_gene_records": join_records(local_core),
                "source_file": str(regions_json),
            })
    out = pd.DataFrame(records) if records else empty_bgc_frame()
    out = coerce_frame_columns(out, BGC_COLUMNS)
    out["bgc_id"] = out["bgc_id"].astype(object)
    blank_ids = out["bgc_id"].astype(str).str.strip() == ""
    out.loc[blank_ids, "bgc_id"] = ["antismash_{0}".format(idx + 1) for idx in out.index[blank_ids]]
else:
    out = empty_bgc_frame()

write_tsv(out, snakemake.output[0])
