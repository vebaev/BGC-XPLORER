from pathlib import Path

import pandas as pd

from common import write_tsv


OUTPUT_COLUMNS = [
    "sample",
    "cgcid",
    "cgc_id",
    "contig",
    "cluster_start",
    "cluster_end",
    "length_bp",
    "genes",
    "cazyme_genes",
    "tc_genes",
    "tf_genes",
    "stp_genes",
    "sulfatase_genes",
    "peptidase_genes",
    "signature_genes",
    "pul_id",
    "pul_substrate",
    "pul_bitscore",
    "signature_pairs",
    "dbcan_sub_substrate",
    "dbcan_sub_substrate_score",
    "predicted_substrate",
    "prediction_source",
    "source_file",
]


def empty_output():
    return pd.DataFrame(columns=OUTPUT_COLUMNS)


def find_file(outdir, names):
    for name in names:
        candidate = outdir / name
        if candidate.exists():
            return candidate
    return None


def clean_text(value):
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if text.lower() == "nan":
        return ""
    return text


def normalize_substrate_frame(path):
    if path is None or not path.exists():
        return pd.DataFrame(columns=[
            "cgcid",
            "pul_id",
            "pul_substrate",
            "pul_bitscore",
            "signature_pairs",
            "dbcan_sub_substrate",
            "dbcan_sub_substrate_score",
        ])
    df = pd.read_csv(path, sep="\t")
    rename_map = {
        "#cgcid": "cgcid",
        "PULID": "pul_id",
        "dbCAN-PUL substrate": "pul_substrate",
        "bitscore": "pul_bitscore",
        "signature pairs": "signature_pairs",
        "dbCAN-sub substrate": "dbcan_sub_substrate",
        "dbCAN-sub substrate score": "dbcan_sub_substrate_score",
    }
    df = df.rename(columns=rename_map)
    keep = list(rename_map.values())
    for column in keep:
        if column not in df.columns:
            df[column] = ""
    return df[keep]


def summarize_cgc_standard(path):
    if path is None or not path.exists():
        return pd.DataFrame(columns=[
            "cgcid",
            "cgc_id",
            "contig",
            "cluster_start",
            "cluster_end",
            "length_bp",
            "genes",
            "cazyme_genes",
            "tc_genes",
            "tf_genes",
            "stp_genes",
            "sulfatase_genes",
            "peptidase_genes",
            "signature_genes",
        ])
    df = pd.read_csv(path, sep="\t")
    required = ["CGC#", "Gene Type", "Contig ID", "Gene Start", "Gene Stop"]
    if any(column not in df.columns for column in required):
        return pd.DataFrame(columns=[
            "cgcid",
            "cgc_id",
            "contig",
            "cluster_start",
            "cluster_end",
            "length_bp",
            "genes",
            "cazyme_genes",
            "tc_genes",
            "tf_genes",
            "stp_genes",
            "sulfatase_genes",
            "peptidase_genes",
            "signature_genes",
        ])

    def summarize_group(group):
        start = int(pd.to_numeric(group["Gene Start"], errors="coerce").min())
        end = int(pd.to_numeric(group["Gene Stop"], errors="coerce").max())
        gene_type = group["Gene Type"].astype(str)
        cazyme = int((gene_type == "CAZyme").sum())
        tc = int((gene_type == "TC").sum())
        tf = int((gene_type == "TF").sum())
        stp = int((gene_type == "STP").sum())
        sulfatase = int((gene_type == "Sulfatase").sum())
        peptidase = int((gene_type == "Peptidase").sum())
        return pd.Series({
            "cgc_id": clean_text(group["CGC#"].iloc[0]),
            "contig": clean_text(group["Contig ID"].iloc[0]),
            "cluster_start": start,
            "cluster_end": end,
            "length_bp": end - start + 1,
            "genes": int(len(group)),
            "cazyme_genes": cazyme,
            "tc_genes": tc,
            "tf_genes": tf,
            "stp_genes": stp,
            "sulfatase_genes": sulfatase,
            "peptidase_genes": peptidase,
            "signature_genes": cazyme + tc + tf + stp + sulfatase + peptidase,
        })

    summary = df.groupby("CGC#", sort=False).apply(summarize_group).reset_index(drop=True)
    summary["cgcid"] = summary["contig"].astype(str) + "|" + summary["cgc_id"].astype(str)
    keep = [
        "cgcid",
        "cgc_id",
        "contig",
        "cluster_start",
        "cluster_end",
        "length_bp",
        "genes",
        "cazyme_genes",
        "tc_genes",
        "tf_genes",
        "stp_genes",
        "sulfatase_genes",
        "peptidase_genes",
        "signature_genes",
    ]
    return summary[keep]


def normalize_summary_frame(path):
    if path is None or not path.exists():
        return summarize_cgc_standard(None)
    df = pd.read_csv(path, sep="\t")
    rename_map = {
        "CGC#": "cgc_id",
        "Contig ID": "contig",
        "Cluster Start": "cluster_start",
        "Cluster End": "cluster_end",
        "Genes": "genes",
        "CAZymes": "cazyme_genes",
        "TC": "tc_genes",
        "TF": "tf_genes",
        "STP": "stp_genes",
        "Sulfatase": "sulfatase_genes",
        "Peptidase": "peptidase_genes",
        "Signatures": "signature_genes",
        "Length (bp)": "length_bp",
    }
    df = df.rename(columns=rename_map)
    keep = list(rename_map.values())
    for column in keep:
        if column not in df.columns:
            df[column] = ""
    df["cgc_id"] = df["cgc_id"].astype(str)
    if "cgcid" not in df.columns:
        df["cgcid"] = df["contig"].astype(str) + "|" + df["cgc_id"].astype(str)
    keep = ["cgcid"] + keep
    return df[keep]


def predicted_substrate(row):
    pul_substrate = clean_text(row.get("pul_substrate"))
    dbsub_substrate = clean_text(row.get("dbcan_sub_substrate"))
    if pul_substrate and dbsub_substrate:
        if pul_substrate == dbsub_substrate:
            return pul_substrate, "dbCAN-PUL + dbCAN-sub"
        return "{0} | {1}".format(pul_substrate, dbsub_substrate), "dbCAN-PUL + dbCAN-sub"
    if pul_substrate:
        return pul_substrate, "dbCAN-PUL only"
    if dbsub_substrate:
        return dbsub_substrate, "dbCAN-sub only"
    return "", "unresolved"


def build_table(sample, outdir):
    substrate_path = find_file(outdir, ["substrate_prediction.tsv"])
    summary_path = find_file(outdir, ["cgc_standard_out_summary.tsv"])
    standard_path = find_file(outdir, ["cgc_standard_out.tsv"])

    substrate = normalize_substrate_frame(substrate_path)
    if summary_path is not None:
        summary = normalize_summary_frame(summary_path)
    else:
        summary = summarize_cgc_standard(standard_path)

    if summary.empty and substrate.empty:
        return empty_output()

    merged = summary.merge(substrate, on="cgcid", how="outer")
    if "cgc_id" not in merged.columns:
        merged["cgc_id"] = merged["cgcid"].astype(str).str.split("|").str[-1]
    if "contig" not in merged.columns:
        merged["contig"] = merged["cgcid"].astype(str).str.split("|").str[0]
    merged["sample"] = sample
    predictions = merged.apply(predicted_substrate, axis=1, result_type="expand")
    merged["predicted_substrate"] = predictions[0]
    merged["prediction_source"] = predictions[1]
    merged["source_file"] = str(substrate_path) if substrate_path is not None else str(summary_path or standard_path or "")

    for column in OUTPUT_COLUMNS:
        if column not in merged.columns:
            merged[column] = ""

    merged = merged[OUTPUT_COLUMNS].sort_values(
        by=["contig", "cluster_start", "cluster_end", "cgcid"],
        na_position="last",
    ).reset_index(drop=True)
    return merged


def main(snakemake_obj):
    sample = snakemake_obj.wildcards.sample
    outdir = Path(snakemake_obj.input.done).parent
    table = build_table(sample, outdir)
    write_tsv(table, snakemake_obj.output[0])


if "snakemake" in globals():
    main(snakemake)
