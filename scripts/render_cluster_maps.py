from html import escape
from pathlib import Path

import pandas as pd

from common import ensure_parent, format_consensus_label, write_tsv


GENE_COLUMNS = [
    "sample",
    "cluster_key",
    "cluster_label",
    "cluster_type",
    "consensus_id",
    "consensus_label",
    "contig",
    "cluster_start",
    "cluster_end",
    "locus_tag",
    "gene_start",
    "gene_end",
    "strand",
    "bakta_gene",
    "bakta_product",
    "eggnog_description",
    "COG_category",
    "preferred_name",
    "ec",
    "kegg_ko",
    "pfams",
    "dbcan_hmm",
    "dbcan_subfamily",
    "dbcan_diamond",
    "dbcan_recommendation",
    "dbcan_substrate",
    "arts_evidence",
    "predictor_core_evidence",
    "gene_category",
    "display_label",
    "tooltip_text",
]


MAP_COLUMNS = [
    "sample",
    "cluster_key",
    "cluster_label",
    "cluster_type",
    "ai_cluster_id",
    "consensus_id",
    "consensus_label",
    "contig",
    "start",
    "end",
    "length_bp",
    "gene_count",
    "category_summary",
    "map_path",
]


CATEGORY_COLORS = {
    "biosynthetic_core": "#8f3f23",
    "cazyme": "#2f855a",
    "transporter": "#2b6cb0",
    "regulator": "#6b46c1",
    "resistance": "#c53030",
    "tailoring_enzyme": "#b7791f",
    "mobile_element": "#718096",
    "hypothetical": "#a0aec0",
    "other": "#4a5568",
}


CATEGORY_LABELS = {
    "biosynthetic_core": "Biosynthetic core",
    "cazyme": "CAZyme / carbohydrate",
    "transporter": "Transporter",
    "regulator": "Regulator",
    "resistance": "Resistance evidence",
    "tailoring_enzyme": "Tailoring enzyme",
    "mobile_element": "Mobile element",
    "hypothetical": "Hypothetical",
    "other": "Other",
}


GENERIC_NAME_TERMS = [
    "hypothetical",
    "unknown",
    "uncharacterized",
    "putative protein",
    "domain-containing protein",
]


def clean(value):
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if text.lower() == "nan" or text == "-":
        return ""
    return text


def is_generic_name(value):
    text = clean(value).lower()
    if not text:
        return True
    if text in {"-", "na", "n/a", "none", "nan"}:
        return True
    return any(term in text for term in GENERIC_NAME_TERMS)


def shorten_name(value, limit=28):
    text = clean(value)
    if not text:
        return ""
    replacements = [
        (" family protein", ""),
        (" domain-containing protein", ""),
        (" protein", ""),
        ("Putative ", ""),
        ("putative ", ""),
    ]
    for old, new in replacements:
        text = text.replace(old, new)
    text = " ".join(text.split())
    if len(text) <= limit:
        return text
    cut = text[:limit].rsplit(" ", 1)[0].strip()
    return cut if cut else text[:limit].strip()


def cazyme_family(value):
    text = clean(value)
    if not text:
        return ""
    first = text.split("+")[0].split("|")[0].split(",")[0].split(";")[0].strip()
    family = first.split("(")[0].strip()
    return family if family and family != "-" else ""


def best_cazyme_family(row):
    for column in ["dbcan_recommendation", "dbcan_hmm", "dbcan_subfamily", "dbcan_diamond"]:
        family = cazyme_family(row.get(column))
        if family:
            return family
    return ""


def load_optional_table(path):
    if path and Path(path).exists():
        return pd.read_csv(path, sep="\t")
    return pd.DataFrame()


def region_key(sample, prefix, raw_id):
    text = clean(raw_id)
    safe = "".join(ch if ch.isalnum() else "_" for ch in text).strip("_")
    return "{0}_{1}_{2}".format(sample, prefix, safe or "region")


def load_dbcan_overview(outdir):
    path = Path(outdir) / "overview.tsv"
    if not path.exists():
        return pd.DataFrame(columns=[
            "locus_tag",
            "dbcan_ec",
            "dbcan_hmm",
            "dbcan_subfamily",
            "dbcan_diamond",
            "dbcan_tools",
            "dbcan_recommendation",
            "dbcan_substrate",
        ])
    df = pd.read_csv(path, sep="\t")
    rename = {
        "Gene ID": "locus_tag",
        "EC#": "dbcan_ec",
        "dbCAN_hmm": "dbcan_hmm",
        "dbCAN_sub": "dbcan_subfamily",
        "DIAMOND": "dbcan_diamond",
        "#ofTools": "dbcan_tools",
        "Recommend Results": "dbcan_recommendation",
        "Substrate": "dbcan_substrate",
    }
    df = df.rename(columns=rename)
    for column in rename.values():
        if column not in df.columns:
            df[column] = ""
    return df[list(rename.values())]


def load_dbcan_sub(outdir):
    path = Path(outdir) / "dbCANsub_hmm_results.tsv"
    if not path.exists():
        return pd.DataFrame(columns=["locus_tag", "dbcan_substrate_detail"])
    df = pd.read_csv(path, sep="\t")
    if "Target Name" not in df.columns:
        return pd.DataFrame(columns=["locus_tag", "dbcan_substrate_detail"])
    if "Substrate" not in df.columns:
        df["Substrate"] = ""
    substrate = (
        df.assign(Substrate=df["Substrate"].fillna("").astype(str))
        .groupby("Target Name", sort=False)["Substrate"]
        .apply(lambda values: "; ".join(sorted({clean(value) for value in values if clean(value)})))
        .reset_index()
        .rename(columns={"Target Name": "locus_tag", "Substrate": "dbcan_substrate_detail"})
    )
    return substrate


def overlaps(start_a, end_a, start_b, end_b):
    return int(start_a) <= int(end_b) and int(end_a) >= int(start_b)


def classify_gene(row):
    text = " ".join([
        clean(row.get("bakta_gene")),
        clean(row.get("bakta_product")),
        clean(row.get("eggnog_description")),
        clean(row.get("preferred_name")),
        clean(row.get("pfams")),
        clean(row.get("COG_category")),
    ]).lower()
    if clean(row.get("arts_evidence")):
        return "resistance"
    if clean(row.get("predictor_core_evidence")):
        return "biosynthetic_core"
    if any(clean(row.get(column)) for column in ["dbcan_hmm", "dbcan_subfamily", "dbcan_diamond", "dbcan_substrate"]):
        return "cazyme"
    if any(term in text for term in [
        "transposase", "integrase", "recombinase", "insertion sequence",
        "mobile", "phage", "tniq",
    ]):
        return "mobile_element"
    if any(term in text for term in [
        "transporter", "permease", "efflux", "abc transporter", "mfs",
        "exporter", "importer", "antiporter", "symporter",
    ]):
        return "transporter"
    if any(term in text for term in [
        "regulator", "transcription", "response regulator", "sigma",
        "histidine kinase", "two-component", "sensor kinase",
    ]) or clean(row.get("COG_category")) in {"K", "T", "KT"}:
        return "regulator"
    if any(term in text for term in [
        "oxidoreductase", "methyltransferase", "hydroxylase", "dehydrogenase",
        "reductase", "transferase", "hydrolase", "monooxygenase",
        "dioxygenase", "aminotransferase", "cytochrome p450", "p450",
    ]):
        return "tailoring_enzyme"
    if "hypothetical" in text or "unknown" in text or not text.strip():
        return "hypothetical"
    return "other"


def tooltip_for(row):
    parts = [
        "Locus: {0}".format(clean(row.get("locus_tag"))),
        "Location: {0}:{1}-{2} ({3})".format(
            clean(row.get("contig")),
            clean(row.get("gene_start")),
            clean(row.get("gene_end")),
            clean(row.get("strand")),
        ),
        "Bakta: {0}".format(clean(row.get("bakta_product")) or "n/a"),
    ]
    if clean(row.get("eggnog_description")):
        parts.append("eggNOG: {0}".format(clean(row.get("eggnog_description"))))
    if clean(row.get("preferred_name")):
        parts.append("Preferred name: {0}".format(clean(row.get("preferred_name"))))
    if clean(row.get("COG_category")):
        parts.append("COG: {0}".format(clean(row.get("COG_category"))))
    if clean(row.get("kegg_ko")):
        parts.append("KEGG: {0}".format(clean(row.get("kegg_ko"))))
    if clean(row.get("pfams")):
        parts.append("PFAMs: {0}".format(clean(row.get("pfams"))))
    if clean(row.get("predictor_core_evidence")):
        parts.append("Core-gene evidence: {0}".format(clean(row.get("predictor_core_evidence"))))
    dbcan_bits = [
        clean(row.get("dbcan_hmm")),
        clean(row.get("dbcan_subfamily")),
        clean(row.get("dbcan_diamond")),
    ]
    dbcan_bits = [bit for bit in dbcan_bits if bit]
    if dbcan_bits:
        parts.append("dbCAN: {0}".format(" | ".join(dbcan_bits)))
    if clean(row.get("dbcan_substrate")):
        parts.append("dbCAN substrate: {0}".format(clean(row.get("dbcan_substrate"))))
    if clean(row.get("arts_evidence")):
        parts.append("ARTS: {0}".format(clean(row.get("arts_evidence"))))
    parts.append("Category: {0}".format(CATEGORY_LABELS.get(clean(row.get("gene_category")), "Other")))
    return " | ".join(parts)


def display_label(row):
    gene = clean(row.get("bakta_gene"))
    preferred = clean(row.get("preferred_name"))
    bakta_product = clean(row.get("bakta_product"))
    eggnog_description = clean(row.get("eggnog_description"))
    candidates = [
        gene if not is_generic_name(gene) else "",
        preferred if not is_generic_name(preferred) else "",
        eggnog_description if not is_generic_name(eggnog_description) else "",
        bakta_product if not is_generic_name(bakta_product) else "",
    ]
    label = next((shorten_name(candidate) for candidate in candidates if clean(candidate)), "")
    if not label:
        label = "hypothetical" if is_generic_name(bakta_product) else "unannotated"
    family = best_cazyme_family(row)
    if family and family not in label:
        label = "{0} ({1})".format(label, family)
    return label


def build_cluster_genes(sample, clusters, bakta, eggnog, dbcan, dbcan_sub, arts):
    if clusters.empty or bakta.empty:
        return pd.DataFrame(columns=GENE_COLUMNS)

    bakta = bakta.copy()
    if "type" in bakta.columns:
        bakta = bakta[bakta["type"].astype(str).str.lower().isin(["cds", "gene", "orf"])]
    bakta["start"] = pd.to_numeric(bakta["start"], errors="coerce")
    bakta["end"] = pd.to_numeric(bakta["end"], errors="coerce")
    bakta = bakta.dropna(subset=["start", "end", "locus_tag"])

    if not eggnog.empty:
        keep = [
            "locus_tag", "eggnog_description", "COG_category", "preferred_name",
            "ec", "KEGG_ko", "pfams",
        ]
        for column in keep:
            if column not in eggnog.columns:
                eggnog[column] = ""
        bakta = bakta.merge(eggnog[keep], on="locus_tag", how="left")
    else:
        for column in ["eggnog_description", "COG_category", "preferred_name", "ec", "KEGG_ko", "pfams"]:
            bakta[column] = ""

    if not dbcan.empty:
        bakta = bakta.merge(dbcan, on="locus_tag", how="left")
    else:
        for column in ["dbcan_hmm", "dbcan_subfamily", "dbcan_diamond", "dbcan_recommendation", "dbcan_substrate"]:
            bakta[column] = ""

    if not dbcan_sub.empty:
        bakta = bakta.merge(dbcan_sub, on="locus_tag", how="left")
        bakta["dbcan_substrate"] = bakta["dbcan_substrate"].fillna("")
        detail = bakta["dbcan_substrate_detail"].fillna("")
        bakta.loc[bakta["dbcan_substrate"].astype(str).isin(["", "-"]), "dbcan_substrate"] = detail
    else:
        bakta["dbcan_substrate_detail"] = ""

    rows = []
    for _, cluster in clusters.iterrows():
        contig = clean(cluster.get("contig")) or clean(cluster.get("contig_id"))
        start = int(float(cluster.get("start")))
        end = int(float(cluster.get("end")))
        genes = bakta[
            (bakta["contig"].astype(str) == contig) &
            (bakta["start"].astype(int) <= end) &
            (bakta["end"].astype(int) >= start)
        ].copy()
        core_by_locus = {}
        for item in clean(cluster.get("core_gene_evidence")).split("; "):
            if " [" in item:
                locus, evidence = item.split(" [", 1)
                core_by_locus[locus.strip()] = evidence.rstrip("]")
        for _, gene in genes.sort_values(["start", "end"]).iterrows():
            arts_evidence = ""
            if not arts.empty:
                hits = arts[
                    (arts["contig"].astype(str) == contig) &
                    (pd.to_numeric(arts["start"], errors="coerce").fillna(-1).astype(int) <= int(gene["end"])) &
                    (pd.to_numeric(arts["end"], errors="coerce").fillna(-1).astype(int) >= int(gene["start"]))
                ]
                if not hits.empty:
                    arts_evidence = "; ".join(sorted({clean(value) for value in hits.get("evidence", []) if clean(value)}))
            row = {
                "sample": sample,
                "cluster_key": clean(cluster.get("cluster_key")) or clean(cluster.get("consensus_id")),
                "cluster_label": clean(cluster.get("cluster_label")) or clean(cluster.get("consensus_id")),
                "cluster_type": clean(cluster.get("cluster_type")) or "consensus",
                "consensus_id": clean(cluster.get("consensus_id")),
                "consensus_label": clean(cluster.get("consensus_label")) or format_consensus_label(cluster.get("consensus_id")),
                "contig": contig,
                "cluster_start": start,
                "cluster_end": end,
                "locus_tag": clean(gene.get("locus_tag")),
                "gene_start": int(gene["start"]),
                "gene_end": int(gene["end"]),
                "strand": clean(gene.get("strand")),
                "bakta_gene": clean(gene.get("gene")),
                "bakta_product": clean(gene.get("product")) or clean(gene.get("bakta_product")),
                "eggnog_description": clean(gene.get("eggnog_description")),
                "COG_category": clean(gene.get("COG_category")),
                "preferred_name": clean(gene.get("preferred_name")),
                "ec": clean(gene.get("ec")),
                "kegg_ko": clean(gene.get("KEGG_ko")) or clean(gene.get("kegg_ko")),
                "pfams": clean(gene.get("pfams")),
                "dbcan_hmm": clean(gene.get("dbcan_hmm")),
                "dbcan_subfamily": clean(gene.get("dbcan_subfamily")),
                "dbcan_diamond": clean(gene.get("dbcan_diamond")),
                "dbcan_recommendation": clean(gene.get("dbcan_recommendation")),
                "dbcan_substrate": clean(gene.get("dbcan_substrate")),
                "arts_evidence": arts_evidence,
                "predictor_core_evidence": core_by_locus.get(clean(gene.get("locus_tag")), ""),
            }
            row["gene_category"] = classify_gene(row)
            row["display_label"] = display_label(row)
            row["tooltip_text"] = tooltip_for(row)
            rows.append(row)
    if not rows:
        return pd.DataFrame(columns=GENE_COLUMNS)
    return pd.DataFrame(rows, columns=GENE_COLUMNS)


def gene_polygon(x1, x2, y, height, strand):
    arrow = min(12, max(5, abs(x2 - x1) * 0.35))
    y1 = y
    ymid = y + height / 2
    y2 = y + height
    if strand == "-":
        points = [
            (x2, y1), (x1 + arrow, y1), (x1, ymid),
            (x1 + arrow, y2), (x2, y2),
        ]
    else:
        points = [
            (x1, y1), (x2 - arrow, y1), (x2, ymid),
            (x2 - arrow, y2), (x1, y2),
        ]
    return " ".join("{0:.1f},{1:.1f}".format(px, py) for px, py in points)


def render_svg(cluster, genes, output_path):
    width = 1500
    margin_left = 44
    margin_right = 34
    map_width = width - margin_left - margin_right
    track_y = 96
    gene_height = 24
    label_y = track_y + gene_height + 18
    legend_height = 86
    height = 286

    cluster_start = int(cluster["start"])
    cluster_end = int(cluster["end"])
    span = max(1, cluster_end - cluster_start + 1)

    def scale(pos):
        return margin_left + ((int(pos) - cluster_start) / span) * map_width

    category_counts = genes["gene_category"].value_counts().to_dict() if not genes.empty else {}
    legend_items = []
    for category, label in CATEGORY_LABELS.items():
        if category in category_counts:
            legend_items.append(
                "<g><rect width='14' height='14' fill='{color}' rx='3'/><text x='20' y='12'>{label}</text></g>".format(
                    color=CATEGORY_COLORS[category],
                    label=escape(label),
                )
            )

    legend_groups = []
    x_cursor = margin_left
    y_cursor = height - legend_height + 30
    for index, item in enumerate(legend_items):
        if index and index % 4 == 0:
            x_cursor = margin_left
            y_cursor += 24
        legend_groups.append("<g transform='translate({0},{1})'>{2}</g>".format(x_cursor, y_cursor, item))
        x_cursor += 250

    gene_shapes = []
    for _, gene in genes.iterrows():
        x1 = max(margin_left, scale(gene["gene_start"]))
        x2 = min(width - margin_right, scale(gene["gene_end"]))
        if x2 - x1 < 4:
            x2 = x1 + 4
        color = CATEGORY_COLORS.get(clean(gene.get("gene_category")), CATEGORY_COLORS["other"])
        label = clean(gene.get("display_label"))
        title = escape(clean(gene.get("tooltip_text")))
        points = gene_polygon(x1, x2, track_y, gene_height, clean(gene.get("strand")))
        text = ""
        if (x2 - x1) >= 20:
            text = (
                "<text x='{x}' y='{y}' text-anchor='start' class='gene-label' "
                "transform='rotate(45 {x} {y})'>{label}</text>"
            ).format(
                x=(x1 + x2) / 2,
                y=label_y,
                label=escape(label[:30]),
            )
        gene_shapes.append(
            "<g class='gene'><title>{title}</title><polygon points='{points}' fill='{color}' />{text}</g>".format(
                title=title,
                points=points,
                color=color,
                text=text,
            )
        )

    ticks = []
    for frac in [0, 0.25, 0.5, 0.75, 1.0]:
        pos = cluster_start + int(span * frac)
        x = scale(pos)
        ticks.append(
            "<g><line x1='{x:.1f}' y1='82' x2='{x:.1f}' y2='90'/><text x='{x:.1f}' y='75' text-anchor='middle'>{label}</text></g>".format(
                x=x,
                label="{:,}".format(pos),
            )
        )

    cluster_label = clean(cluster.get("cluster_label")) or clean(cluster.get("consensus_id"))

    svg = """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">
<title id="title">{cluster_id} gene map</title>
<desc id="desc">Horizontal gene map colored by functional category. Hover genes for Bakta, eggNOG, dbCAN and ARTS details.</desc>
<style>
  .bg {{ fill: #fbfcff; stroke: rgba(94, 108, 152, 0.18); stroke-width: 1.25; }}
  .title {{ font: 700 24px -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; fill: #161f33; }}
  .meta {{ font: 14px -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; fill: #5e6c8f; }}
  .axis {{ stroke: #b9c3dc; stroke-width: 2; }}
  .gene polygon {{ stroke: rgba(39, 49, 73, 0.22); stroke-width: 1; }}
  .gene-label {{ font: 9.5px -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; fill: #44506b; pointer-events: none; }}
  text {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; fill: #5e6c8f; }}
</style>
<rect class="bg" x="0" y="0" width="{width}" height="{height}" rx="14"/>
<text class="title" x="{left}" y="34">{cluster_id}</text>
<text class="meta" x="{left}" y="58">{location} | {gene_count} genes | {length} bp</text>
<line class="axis" x1="{left}" y1="86" x2="{right}" y2="86"/>
{ticks}
{genes}
<text class="meta" x="{left}" y="{legend_title_y}">Gene category colors</text>
{legend}
</svg>
""".format(
        width=width,
        height=height,
        cluster_id=escape(cluster_label),
        left=margin_left,
        right=width - margin_right,
        location=escape("{0}:{1:,}-{2:,}".format(clean(cluster.get("contig")), cluster_start, cluster_end)),
        gene_count=len(genes),
        length="{:,}".format(span),
        ticks="\n".join(ticks),
        genes="\n".join(gene_shapes),
        legend_title_y=height - legend_height + 12,
        legend="\n".join(legend_groups),
    )
    ensure_parent(str(output_path))
    output_path.write_text(svg, encoding="utf-8")


def category_summary(genes):
    if genes.empty:
        return ""
    counts = genes["gene_category"].value_counts().to_dict()
    parts = []
    for category, count in sorted(counts.items(), key=lambda item: (-item[1], item[0])):
        parts.append("{0}: {1}".format(CATEGORY_LABELS.get(category, category), count))
    return "; ".join(parts)


def main(snakemake_obj):
    sample = snakemake_obj.wildcards.sample
    clusters = load_optional_table(snakemake_obj.input.evidence)
    dbcan_summary = load_optional_table(snakemake_obj.input.dbcan_summary)
    bakta = load_optional_table(snakemake_obj.input.bakta)
    eggnog = load_optional_table(snakemake_obj.input.eggnog)
    arts = load_optional_table(snakemake_obj.input.arts)
    dbcan_outdir = Path(snakemake_obj.input.dbcan_done).parent
    dbcan = load_dbcan_overview(dbcan_outdir)
    dbcan_sub = load_dbcan_sub(dbcan_outdir)
    map_dir = Path(snakemake_obj.output.map_assets)
    map_dir.mkdir(parents=True, exist_ok=True)

    consensus_regions = clusters.copy()
    if not consensus_regions.empty:
        consensus_regions["cluster_key"] = consensus_regions["consensus_id"].astype(str)
        if "consensus_label" not in consensus_regions.columns:
            consensus_regions["consensus_label"] = consensus_regions["consensus_id"].astype(str).map(format_consensus_label)
        consensus_regions["cluster_label"] = consensus_regions["consensus_label"].astype(str)
        consensus_regions["cluster_type"] = "consensus"

    dbcan_regions = pd.DataFrame()
    if not dbcan_summary.empty:
        dbcan_regions = dbcan_summary.copy()
        dbcan_regions["predicted_substrate"] = dbcan_regions.get("predicted_substrate", "").fillna("").astype(str).str.strip()
        dbcan_regions = dbcan_regions[
            (dbcan_regions["predicted_substrate"] != "")
            & (dbcan_regions["predicted_substrate"].str.lower() != "unresolved")
        ].copy()
        if not dbcan_regions.empty:
            dbcan_regions["start"] = pd.to_numeric(dbcan_regions["cluster_start"], errors="coerce")
            dbcan_regions["end"] = pd.to_numeric(dbcan_regions["cluster_end"], errors="coerce")
            dbcan_regions["cluster_key"] = dbcan_regions["cgc_id"].astype(str).map(
                lambda value: region_key(sample, "cgc", value)
            )
            dbcan_regions["cluster_label"] = dbcan_regions["cgc_id"].astype(str).map(
                lambda value: "{0}_{1}".format(sample, value)
            )
            dbcan_regions["cluster_type"] = "dbcan_cgc"
            dbcan_regions["consensus_id"] = ""

    region_frames = [frame for frame in [consensus_regions, dbcan_regions] if not frame.empty]
    all_regions = pd.concat(region_frames, ignore_index=True, sort=False) if region_frames else pd.DataFrame()

    cluster_genes = build_cluster_genes(sample, all_regions, bakta, eggnog, dbcan, dbcan_sub, arts)
    write_tsv(cluster_genes, snakemake_obj.output.cluster_genes)

    map_rows = []
    for _, cluster in all_regions.iterrows():
        cluster_key = clean(cluster.get("cluster_key")) or clean(cluster.get("consensus_id"))
        if not cluster_key:
            continue
        genes = cluster_genes[cluster_genes["cluster_key"] == cluster_key].copy()
        svg_path = map_dir / "{0}.svg".format(cluster_key)
        render_svg(cluster, genes, svg_path)
        map_rows.append({
            "sample": sample,
            "cluster_key": cluster_key,
            "cluster_label": clean(cluster.get("cluster_label")) or cluster_key,
            "cluster_type": clean(cluster.get("cluster_type")) or "consensus",
            "ai_cluster_id": clean(cluster.get("consensus_id")),
            "consensus_id": clean(cluster.get("consensus_id")),
            "consensus_label": clean(cluster.get("consensus_label")) or format_consensus_label(cluster.get("consensus_id")),
            "contig": clean(cluster.get("contig")),
            "start": int(float(cluster.get("start"))),
            "end": int(float(cluster.get("end"))),
            "length_bp": int(float(cluster.get("length_bp", 0))) if clean(cluster.get("length_bp")) else int(float(cluster.get("end"))) - int(float(cluster.get("start"))) + 1,
            "gene_count": len(genes),
            "category_summary": category_summary(genes),
            "map_path": "assets/cluster_maps/{0}.svg".format(cluster_key),
        })
    maps = pd.DataFrame(map_rows, columns=MAP_COLUMNS)
    write_tsv(maps, snakemake_obj.output.cluster_maps)


if "snakemake" in globals():
    main(snakemake)
