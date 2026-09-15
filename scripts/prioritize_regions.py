import pandas as pd

from mibig_hits import novelty_score_for_priority

from common import add_consensus_label, load_table_if_exists, write_tsv


def to_num(series):
    return pd.to_numeric(series, errors="coerce")


def interval_overlap(a_start, a_end, b_start, b_end):
    left = max(a_start, b_start)
    right = min(a_end, b_end)
    return max(0, right - left + 1)


def unique_preserve_order(values):
    seen = []
    for value in values:
        text = str(value).strip()
        if not text or text.lower() == "nan":
            continue
        if text not in seen:
            seen.append(text)
    return seen


def summarize_arts_hits(hit_rows):
    if hit_rows.empty:
        return ""
    evidence_labels = {
        "known_hit": "known resistance hit",
        "duf_hit": "DUF-associated hit",
        "core_gene": "core-gene signal",
        "duplication_signal": "duplication signal",
    }
    evidence_counts = []
    for evidence, group in hit_rows.groupby("evidence", dropna=False):
        label = evidence_labels.get(str(evidence).strip(), str(evidence).strip() or "ARTS signal")
        features = unique_preserve_order(group["feature"].tolist())[:3]
        feature_text = ", ".join(features) if features else "unspecified features"
        evidence_counts.append(
            "{count} {label} ({features})".format(
                count=len(group),
                label=label,
                features=feature_text,
            )
        )
    return "ARTS overlap: " + "; ".join(evidence_counts)


def build_why_prioritized(row):
    reasons = []
    support_count = int(row.get("num_supporting_tools", 0) or 0)
    arts_hits = int(row.get("arts_hits", 0) or 0)
    support_tools = unique_preserve_order(str(row.get("supporting_tools", "")).split(","))
    dereplication_status = str(row.get("dereplication_status", "")).strip()
    novelty_score = row.get("novelty_score", "")
    edge_truncated = bool(row.get("edge_truncated", False))
    products_text = str(row.get("products", "")).lower()

    if support_count >= 3:
        reasons.append("{0} tools agree".format(support_count))
    elif support_count == 2:
        reasons.append("2 tools agree")
    else:
        reasons.append("single-tool only")

    if arts_hits > 0:
        reasons.append("ARTS support ({0} hits)".format(arts_hits))
    else:
        reasons.append("no ARTS support")

    if "shared_core_gene" in str(row.get("overlap_relationship", "")):
        reasons.append("shared core-gene support")

    if dereplication_status == "novel_candidate":
        reasons.append("no close MIBiG match")
    elif dereplication_status == "divergent":
        reasons.append("divergent from known MIBiG cluster")
    elif dereplication_status == "related":
        reasons.append("related to known MIBiG cluster")
    elif dereplication_status == "known-like":
        reasons.append("known-like MIBiG match")
    elif novelty_score != "" and not pd.isna(novelty_score):
        reasons.append("MIBiG status unknown")

    if "antibacterial" in products_text:
        reasons.append("antibacterial signal")
    if "cytotoxic" in products_text:
        reasons.append("cytotoxic signal")
    if edge_truncated:
        reasons.append("contig-edge region")
    if support_tools:
        reasons.append("callers: {0}".format(", ".join(support_tools)))

    return " + ".join(reasons)


def build_why_not_prioritized(row):
    reasons = []
    support_count = int(row.get("num_supporting_tools", 0) or 0)
    arts_hits = int(row.get("arts_hits", 0) or 0)
    edge_truncated = bool(row.get("edge_truncated", False))
    dereplication_status = str(row.get("dereplication_status", "")).strip()
    boundary_confidence = str(row.get("boundary_confidence", "")).strip()
    priority_class = str(row.get("priority_class", "")).strip()

    if support_count <= 1:
        reasons.append("single-tool only")
    elif support_count == 2:
        reasons.append("only 2 tools support this locus")

    if arts_hits == 0:
        reasons.append("no ARTS support")
    elif support_count <= 1 and arts_hits > 0:
        reasons.append("ARTS support without cross-tool confirmation")

    if boundary_confidence in ["low", "medium_edge"]:
        reasons.append("uncertain region boundaries")

    if edge_truncated:
        reasons.append("near contig edge; possible partial BGC")

    if dereplication_status == "known-like":
        reasons.append("strong similarity to known MIBiG cluster")
    elif dereplication_status == "related":
        reasons.append("close to a known MIBiG family")

    if priority_class == "low" and not reasons:
        reasons.append("limited supporting evidence")

    return "; ".join(unique_preserve_order(reasons))


consensus = pd.read_csv(snakemake.input.consensus, sep="\t")
arts = load_table_if_exists(snakemake.input.arts, [
    "sample", "tool", "contig", "start", "end", "feature", "score", "evidence", "source_file"
])
mibig = pd.read_csv(snakemake.input.mibig, sep="\t") if snakemake.input.mibig else pd.DataFrame()
weights = snakemake.config["consensus"].get("score_weights", {})
priority_thresholds = snakemake.config["consensus"].get("priority_thresholds", {})

if consensus.empty:
    prioritized = pd.DataFrame(columns=[
        "sample", "consensus_id", "consensus_label", "priority_score", "priority_class", "confidence_category",
        "interest_category", "arts_hits", "notes", "why_prioritized", "why_not_prioritized",
        "recommended_followup"
    ])
else:
    prioritized = consensus.copy()
    prioritized["priority_score"] = prioritized["num_supporting_tools"].astype(float) * float(weights.get("tool_support", 1.0))
    prioritized["priority_class"] = ""
    prioritized["confidence_category"] = ""
    prioritized["interest_category"] = ""
    prioritized["arts_hits"] = 0
    prioritized["notes"] = ""
    prioritized["why_prioritized"] = ""
    prioritized["why_not_prioritized"] = ""
    prioritized["recommended_followup"] = ""
    prioritized["best_mibig_id"] = ""
    prioritized["best_mibig_product"] = ""
    prioritized["best_mibig_class"] = ""
    prioritized["mibig_similarity"] = ""
    prioritized["match_score"] = ""
    prioritized["score_metric"] = ""
    prioritized["matched_genes"] = ""
    prioritized["core_gene_hits"] = ""
    prioritized["dereplication_status"] = ""
    prioritized["novelty_score"] = ""
    prioritized["evidence_source"] = ""

    for idx, row in prioritized.iterrows():
        support_count = int(row.get("num_supporting_tools", 0) or 0)
        if support_count >= 3:
            prioritized.at[idx, "priority_score"] += float(weights.get("tool_triple_bonus", 3.0))
        elif support_count == 2:
            prioritized.at[idx, "priority_score"] += float(weights.get("tool_pair_bonus", 1.5))
        else:
            prioritized.at[idx, "priority_score"] -= float(weights.get("single_tool_penalty", 2.5))

        derep = mibig[mibig["consensus_id"] == row["consensus_id"]] if not mibig.empty else pd.DataFrame()
        if not derep.empty:
            derep_row = derep.iloc[0]
            prioritized.at[idx, "best_mibig_id"] = derep_row.get("best_mibig_id", "")
            prioritized.at[idx, "best_mibig_product"] = derep_row.get("best_mibig_product", "")
            prioritized.at[idx, "best_mibig_class"] = derep_row.get("best_mibig_class", "")
            prioritized.at[idx, "mibig_similarity"] = derep_row.get("mibig_similarity", "")
            prioritized.at[idx, "match_score"] = derep_row.get("match_score", "")
            prioritized.at[idx, "score_metric"] = derep_row.get("score_metric", "")
            prioritized.at[idx, "matched_genes"] = derep_row.get("matched_genes", "")
            prioritized.at[idx, "core_gene_hits"] = derep_row.get("core_gene_hits", "")
            prioritized.at[idx, "dereplication_status"] = derep_row.get("dereplication_status", "")
            prioritized.at[idx, "evidence_source"] = derep_row.get("evidence_source", "")
            novelty_value = novelty_score_for_priority(derep_row.get("novelty_score", ""))
            prioritized.at[idx, "novelty_score"] = (
                novelty_value if derep_row.get("evidence_source", "") == "comparippson_html" else ""
            )
            prioritized.at[idx, "priority_score"] += float(weights.get("novelty_score", 3.0)) * novelty_value

        if arts.empty or not row["contig"]:
            local = pd.DataFrame()
        else:
            local = arts[arts["contig"] == row["contig"]].copy()
        overlapping_hits = []
        hits = 0
        if not local.empty:
            local["start_num"] = to_num(local["start"])
            local["end_num"] = to_num(local["end"])
            for _, hit in local.iterrows():
                if pd.isna(hit["start_num"]) or pd.isna(hit["end_num"]) or pd.isna(row["start"]) or pd.isna(row["end"]):
                    continue
                if interval_overlap(row["start"], row["end"], hit["start_num"], hit["end_num"]) > 0:
                    hits += 1
                    overlapping_hits.append(hit)
            prioritized.at[idx, "arts_hits"] = hits
            capped_hits = min(hits, float(weights.get("arts_signal_cap", 4.0)))
            prioritized.at[idx, "priority_score"] += float(weights.get("arts_signal", 0.75)) * capped_hits

        if row["boundary_confidence"] == "high":
            prioritized.at[idx, "priority_score"] += float(weights.get("boundary_high", 1.5))
        elif row["boundary_confidence"] == "medium":
            prioritized.at[idx, "priority_score"] += float(weights.get("boundary_medium", 0.5))
        if "shared_core_gene" in str(row.get("overlap_relationship", "")):
            prioritized.at[idx, "priority_score"] += float(weights.get("shared_core_gene_bonus", 1.5))

        products_text = str(row.get("products", "")).lower()
        if "antibacterial" in products_text:
            prioritized.at[idx, "priority_score"] += float(weights.get("antibacterial_bonus", 1.5))
        if "cytotoxic" in products_text:
            prioritized.at[idx, "priority_score"] += float(weights.get("cytotoxic_bonus", 1.0))
        if row["edge_truncated"]:
            prioritized.at[idx, "priority_score"] -= float(weights.get("edge_penalty", 1.5))
            prioritized.at[idx, "recommended_followup"] = "reassemble or inspect flanking sequence before deprioritizing biologically"
        else:
            prioritized.at[idx, "recommended_followup"] = "manual gene-neighborhood review and targeted dereplication"

        if hits:
            hit_frame = pd.DataFrame(overlapping_hits)
            prioritized.at[idx, "notes"] = summarize_arts_hits(hit_frame)

        if prioritized.at[idx, "dereplication_status"] == "known-like":
            prioritized.at[idx, "recommended_followup"] = "review if this is a known-like recovery target before novelty-focused follow-up"
        elif prioritized.at[idx, "dereplication_status"] in ["divergent", "novel_candidate"]:
            prioritized.at[idx, "recommended_followup"] = "prioritize dereplication review, manual domain inspection, and targeted validation"

        score = float(prioritized.at[idx, "priority_score"])
        novelty_value = prioritized.at[idx, "novelty_score"]
        novelty = float(novelty_value) if novelty_value != "" and not pd.isna(novelty_value) else None
        if support_count >= 3 and row["boundary_confidence"] in ["high", "medium"] and not row["edge_truncated"]:
            prioritized.at[idx, "confidence_category"] = "high-confidence BGC"
        elif support_count >= 2:
            prioritized.at[idx, "confidence_category"] = "supported candidate"
        else:
            prioritized.at[idx, "confidence_category"] = "single-tool candidate"

        if (
            prioritized.at[idx, "dereplication_status"] in ["divergent", "novel_candidate"]
            or (novelty is not None and novelty >= 0.75)
            or hits >= 3
            or "cytotoxic" in products_text
        ):
            prioritized.at[idx, "interest_category"] = "high-interest / potentially novel"
        else:
            prioritized.at[idx, "interest_category"] = "standard follow-up"

        if score >= float(priority_thresholds.get("high", 12.0)):
            prioritized.at[idx, "priority_class"] = "high"
        elif score >= float(priority_thresholds.get("medium", 7.0)):
            prioritized.at[idx, "priority_class"] = "medium"
        else:
            prioritized.at[idx, "priority_class"] = "low"

        prioritized.at[idx, "why_prioritized"] = build_why_prioritized(prioritized.loc[idx])
        prioritized.at[idx, "why_not_prioritized"] = build_why_not_prioritized(prioritized.loc[idx])

    prioritized = prioritized.sort_values(["priority_score", "num_supporting_tools"], ascending=[False, False])

prioritized = add_consensus_label(prioritized)

write_tsv(prioritized, snakemake.output.prioritized)
write_tsv(prioritized, snakemake.output.prioritized_regions)
