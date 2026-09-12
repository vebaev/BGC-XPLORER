from pathlib import Path

import pandas as pd

from common import BGC_COLUMNS, coerce_frame_columns, empty_bgc_frame, first_existing_column, first_existing_path, write_tsv


sample = snakemake.wildcards.sample
outdir = Path(snakemake.input.done).parent
source = first_existing_path(outdir, ["**/*bgc*.tsv", "**/*cluster*.tsv", "*.bgc.tsv", "*.clusters.tsv"])

if source is not None:
    df = pd.read_csv(str(source), sep="\t")
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
        "source_file": str(source),
    })
    out = coerce_frame_columns(out, BGC_COLUMNS)
else:
    out = empty_bgc_frame()

write_tsv(out, snakemake.output[0])
