from pathlib import Path

import pandas as pd

from common import BGC_COLUMNS, coerce_frame_columns, empty_bgc_frame, first_existing_column, first_existing_path, write_tsv


sample = snakemake.wildcards.sample
outdir = Path(snakemake.input.done).parent
source = first_existing_path(outdir, ["*.clusters.tsv", "**/*.clusters.tsv", "*clusters*.tsv", "**/*clusters*.tsv"])

if source is not None:
    df = pd.read_csv(str(source), sep="\t")
    renamed = pd.DataFrame({
        "sample": sample,
        "tool": "gecco",
        "contig": first_existing_column(df, ["sequence_id", "contig_id", "seq_id"]),
        "start": first_existing_column(df, ["start", "cluster_start"]),
        "end": first_existing_column(df, ["end", "cluster_end"]),
        "strand": first_existing_column(df, ["strand"], "."),
        "bgc_id": first_existing_column(df, ["cluster_id", "id"]),
        "bgc_type": first_existing_column(df, ["type", "biosyn_class"], "unknown"),
        "product": first_existing_column(df, ["type", "biosyn_class"], ""),
        "score": first_existing_column(df, ["average_p", "mean_probability", "score"], ""),
        "confidence": first_existing_column(df, ["max_p", "max_probability", "probability"], ""),
        "source_file": str(source),
    })
    out = coerce_frame_columns(renamed, BGC_COLUMNS)
else:
    out = empty_bgc_frame()

write_tsv(out, snakemake.output[0])
