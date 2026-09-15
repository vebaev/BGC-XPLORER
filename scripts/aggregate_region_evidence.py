"""Join observable tool, ARTS, and MIBiG evidence without ranking regions."""

import os

import pandas as pd

from common import add_consensus_label, load_table_if_exists, write_tsv


def interval_overlap(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def arts_counts_for_region(region, arts):
    counts = {"known_hit": 0, "duf_hit": 0, "other": 0}
    if arts.empty:
        return counts
    same_contig = arts[arts["contig"].astype(str) == str(region["contig"])]
    if same_contig.empty:
        return counts
    start, end = pd.to_numeric([region["start"], region["end"]], errors="coerce")
    if pd.isna(start) or pd.isna(end):
        return counts
    seen = set()
    for _, hit in same_contig.iterrows():
        hit_start, hit_end = pd.to_numeric([hit.get("start"), hit.get("end")], errors="coerce")
        if pd.isna(hit_start) or pd.isna(hit_end):
            continue
        if interval_overlap(start, end, hit_start, hit_end) <= 0:
            continue
        evidence = str(hit.get("evidence", "")).strip()
        key = (evidence, str(hit.get("feature", "")).strip(), int(hit_start), int(hit_end))
        if key in seen:
            continue
        seen.add(key)
        category = evidence if evidence in ("known_hit", "duf_hit") else "other"
        counts[category] += 1
    return counts


def aggregate_evidence(consensus, arts, mibig):
    if consensus.empty:
        return consensus.copy()
    result = consensus.copy()
    result["arts_known_hits"] = 0
    result["arts_duf_hits"] = 0
    result["arts_other_hits"] = 0
    for index, row in result.iterrows():
        counts = arts_counts_for_region(row, arts)
        result.at[index, "arts_known_hits"] = counts["known_hit"]
        result.at[index, "arts_duf_hits"] = counts["duf_hit"]
        result.at[index, "arts_other_hits"] = counts["other"]
    result["arts_hits"] = result[["arts_known_hits", "arts_duf_hits", "arts_other_hits"]].sum(axis=1)

    mibig_fields = [
        "best_mibig_id", "best_mibig_product", "best_mibig_class",
        "mibig_similarity", "match_score", "score_metric", "matched_genes",
        "core_gene_hits", "evidence_source",
    ]
    if not mibig.empty and "consensus_id" in mibig.columns:
        available = [field for field in mibig_fields if field in mibig.columns]
        result = result.merge(
            mibig[["consensus_id"] + available].drop_duplicates("consensus_id"),
            on="consensus_id", how="left", validate="one_to_one",
        )
    for field in mibig_fields:
        if field not in result.columns:
            result[field] = ""
    left = pd.to_numeric(result.get("distance_to_left_edge"), errors="coerce")
    right = pd.to_numeric(result.get("distance_to_right_edge"), errors="coerce")
    result["nearest_contig_edge_bp"] = pd.concat([left, right], axis=1).min(axis=1)
    result["nearest_contig_edge_bp"] = result["nearest_contig_edge_bp"].fillna("")
    result = result.sort_values(["contig", "start", "end"], kind="mergesort")
    return add_consensus_label(result)


if "snakemake" in globals():
    consensus = pd.read_csv(snakemake.input.consensus, sep="\t")
    arts = load_table_if_exists(snakemake.input.arts, [
        "sample", "tool", "contig", "start", "end", "feature", "score", "evidence", "source_file"
    ])
    mibig = pd.read_csv(snakemake.input.mibig, sep="\t") if os.path.exists(snakemake.input.mibig) else pd.DataFrame()
    write_tsv(aggregate_evidence(consensus, arts, mibig), snakemake.output.evidence)
