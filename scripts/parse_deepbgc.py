from pathlib import Path

import pandas as pd

from common import BGC_COLUMNS, coerce_frame_columns, empty_bgc_frame, first_existing_column, first_existing_path, write_tsv
from core_gene_evidence import domain_core_record, join_records


sample = snakemake.wildcards.sample
outdir = Path(snakemake.input.done).parent
source = first_existing_path(outdir, ["**/*bgc*.tsv", "**/*cluster*.tsv", "*.bgc.tsv", "*.clusters.tsv"])

if source is not None:
    df = pd.read_csv(str(source), sep="\t")
    pfam_path = first_existing_path(outdir, ["*.pfam.tsv", "**/*.pfam.tsv"])
    pfam = pd.read_csv(str(pfam_path), sep="\t") if pfam_path is not None else pd.DataFrame()

    def core_records(cluster):
        if pfam.empty:
            return ""
        biological = {
            value.strip().split(".")[0]
            for value in str(cluster.get("bio_pfam_ids", "")).split(";")
            if value.strip() and value.strip().lower() != "nan"
        }
        if not biological:
            return ""
        local = pfam[
            (pfam["sequence_id"].astype(str) == str(cluster.get("sequence_id", "")))
            & (pd.to_numeric(pfam["gene_start"], errors="coerce") <= float(cluster.get("nucl_end", 0)))
            & (pd.to_numeric(pfam["gene_end"], errors="coerce") >= float(cluster.get("nucl_start", 0)))
        ]
        records = []
        for protein_id, protein in local.groupby("protein_id", sort=False):
            domains = [domain for domain in protein["pfam_id"].tolist() if str(domain).split(".")[0] in biological]
            record = domain_core_record(
                protein_id, protein["gene_start"].min(), protein["gene_end"].max(),
                "deepbgc", domains,
            )
            if record:
                records.append(record)
        return join_records(records)

    out = pd.DataFrame({
        "sample": sample,
        "tool": "deepbgc",
        "contig": first_existing_column(df, ["sequence_id", "contig_id", "sequence"]),
        "start": first_existing_column(df, ["nucl_start", "start"]),
        "end": first_existing_column(df, ["nucl_end", "end"]),
        "strand": first_existing_column(df, ["strand"], "."),
        "bgc_id": first_existing_column(df, ["bgc_candidate_id", "cluster_id", "id"]),
        "bgc_type": first_existing_column(df, ["product_class", "classification"], "unknown"),
        "product": first_existing_column(df, ["product_activity", "classification", "product_class"], ""),
        "score": first_existing_column(df, ["detector_score", "score"], ""),
        "confidence": first_existing_column(df, ["classifier_score", "probability"], ""),
        "core_gene_records": df.apply(core_records, axis=1),
        "source_file": str(source),
    })
    out = coerce_frame_columns(out, BGC_COLUMNS)
else:
    out = empty_bgc_frame()

write_tsv(out, snakemake.output[0])
