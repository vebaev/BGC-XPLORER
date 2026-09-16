from html import escape

import pandas as pd

from common import df_to_html_table, format_consensus_label, html_page, load_table_if_exists, write_tsv


def to_num(series):
    return pd.to_numeric(series, errors="coerce")


def interval_overlap(a_start, a_end, b_start, b_end):
    left = max(a_start, b_start)
    right = min(a_end, b_end)
    return max(0, right - left + 1)


def unique_join(values):
    seen = []
    for value in values:
        text = str(value).strip()
        if not text or text.lower() == "nan":
            continue
        if text not in seen:
            seen.append(text)
    return "; ".join(seen)


def stat_card(label, value, note):
    return (
        "<div class='stat'>"
        "<div class='label'>{label}</div>"
        "<div class='value'>{value}</div>"
        "<p>{note}</p>"
        "</div>"
    ).format(label=escape(label), value=escape(str(value)), note=escape(note))


def tool_support_summary(tool_df, contig, start, end):
    if tool_df.empty:
        return "", "", ""
    local = tool_df[tool_df["contig"] == contig].copy()
    if local.empty:
        return "", "", ""
    local["start_num"] = to_num(local["start"])
    local["end_num"] = to_num(local["end"])
    local = local.dropna(subset=["start_num", "end_num"])
    local = local[
        local.apply(
            lambda row: interval_overlap(start, end, row["start_num"], row["end_num"]) > 0,
            axis=1,
        )
    ]
    if local.empty:
        return "", "", ""
    ids = unique_join(local["bgc_id"].tolist())
    classes = unique_join(local["bgc_type"].tolist())
    products = unique_join(local["product"].tolist())
    details = unique_join([classes, products])
    return ids, classes, details


sample = snakemake.wildcards.sample
merged = pd.read_csv(snakemake.input.evidence, sep="\t")
antismash = load_table_if_exists(snakemake.input.antismash, [
    "sample", "tool", "contig", "start", "end", "strand", "bgc_id", "bgc_type", "product", "score", "confidence", "source_file"
])
gecco = load_table_if_exists(snakemake.input.gecco, antismash.columns.tolist())
deepbgc = load_table_if_exists(snakemake.input.deepbgc, antismash.columns.tolist())

rows = []
for _, row in merged.iterrows():
    start = int(row["start"])
    end = int(row["end"])
    antismash_ids, antismash_types, antismash_details = tool_support_summary(antismash, row["contig"], start, end)
    gecco_ids, gecco_types, gecco_details = tool_support_summary(gecco, row["contig"], start, end)
    deepbgc_ids, deepbgc_types, deepbgc_details = tool_support_summary(deepbgc, row["contig"], start, end)
    rows.append({
        "sample": sample,
        "consensus_id": row["consensus_id"],
        "consensus_label": row.get("consensus_label", format_consensus_label(row["consensus_id"])),
        "contig": row["contig"],
        "start": start,
        "end": end,
        "support_tools": row["support_tools"],
        "support_count": row["support_count"],
        "antismash_ids": antismash_ids,
        "antismash_types": antismash_types,
        "antismash_details": antismash_details,
        "gecco_ids": gecco_ids,
        "gecco_types": gecco_types,
        "gecco_details": gecco_details,
        "deepbgc_ids": deepbgc_ids,
        "deepbgc_types": deepbgc_types,
        "deepbgc_details": deepbgc_details,
        "merged_bgc_types": row.get("bgc_types", ""),
        "merged_products": row.get("products", ""),
        "biological_interpretation": row.get("biological_interpretation", ""),
        "core_gene_support": row.get("core_gene_support", ""),
        "core_gene_evidence": row.get("core_gene_evidence", ""),
        "arts_hits": row.get("arts_hits", ""),
        "arts_known_hits": row.get("arts_known_hits", ""),
        "arts_duf_hits": row.get("arts_duf_hits", ""),
        "arts_other_hits": row.get("arts_other_hits", ""),
        "nearest_contig_edge_bp": row.get("nearest_contig_edge_bp", ""),
        "best_mibig_id": row.get("best_mibig_id", ""),
        "best_mibig_product": row.get("best_mibig_product", ""),
        "best_mibig_class": row.get("best_mibig_class", ""),
        "mibig_similarity": row.get("mibig_similarity", ""),
        "match_score": row.get("match_score", ""),
        "score_metric": row.get("score_metric", ""),
        "matched_genes": row.get("matched_genes", ""),
        "core_gene_hits": row.get("core_gene_hits", ""),
        "evidence_source": row.get("evidence_source", ""),
    })

catalog = pd.DataFrame(rows).sort_values(
    ["contig", "start", "end"],
    ascending=[True, True, True],
)

write_tsv(catalog, snakemake.output.tsv)

multi_tool = catalog[pd.to_numeric(catalog["support_count"], errors="coerce").fillna(0) > 1]
sections = [
    (
        "<section class='hero'>"
        "<div class='eyebrow'>Cluster Catalog</div>"
        "<h1>{sample}</h1>"
        "<p>One row per grouped candidate locus. Each tool-specific block preserves its original prediction; grouped boundaries are approximate.</p>"
        "</section>"
    ).format(sample=escape(sample)),
    "<section class='section'><h2>Catalog Summary</h2><div class='stats'>{cards}</div></section>".format(
        cards="".join([
            stat_card("Grouped loci", len(catalog), "approximate candidate loci in the sample"),
            stat_card("Multi-caller loci", len(multi_tool), "contain predictions from at least two callers"),
            stat_card("antiSMASH-backed", int((catalog["antismash_ids"].fillna("").str.len() > 0).sum()), "catalog rows with antiSMASH evidence"),
            stat_card("GECCO-backed", int((catalog["gecco_ids"].fillna("").str.len() > 0).sum()), "catalog rows with GECCO evidence"),
            stat_card("DeepBGC-backed", int((catalog["deepbgc_ids"].fillna("").str.len() > 0).sum()), "catalog rows with DeepBGC evidence"),
        ])
    ),
    (
        "<section class='section'>"
        "<h2>Full Cluster Table</h2>"
        "<p class='muted'>When more than one tool supports a region, the row keeps the combined interpretation so differences stay visible instead of being flattened away.</p>"
        "{table}"
        "</section>"
    ).format(table=df_to_html_table(catalog)),
]

html = html_page("Cluster Catalog - {sample}".format(sample=sample), sections)
with open(snakemake.output.html, "w", encoding="utf-8") as handle:
    handle.write(html)
