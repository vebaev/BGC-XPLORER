import json
import os
from datetime import datetime
from html import escape
from urllib.parse import quote

import pandas as pd

from common import html_page, load_table_if_exists, read_json
from report_branding import image_data_uri, report_home_link
from report_design import DONUT_COLORS, PRIMARY_GLANCE_LABELS, reproducibility_panel
from report_evidence_table import grouping_description, location_html, render_evidence_table


def _load_logo_data_uri():
    report_dir = os.path.dirname(os.path.abspath(str(snakemake.output[0])))
    return image_data_uri([
        os.path.join(report_dir, "assets", "logo.png"),
        os.path.join(report_dir, "assets", "logo@56.png"),
        "/app/logo.jpg",
        os.path.abspath("logo.jpg"),
    ])


LOGO_DATA_URI = _load_logo_data_uri()


def as_text(value, fallback="n/a"):
    if pd.isna(value):
        return fallback
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return fallback
    return text


def as_int(value, fallback=0):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return fallback


def render_html_table(df, html_columns=None):
    if df.empty:
        return "<p>No records found.</p>"
    html_columns = set(html_columns or [])
    header = "".join("<th>{0}</th>".format(escape(str(column))) for column in df.columns)
    body_rows = []
    for _, row in df.iterrows():
        cells = []
        for column in df.columns:
            value = row[column]
            if column in html_columns:
                cells.append("<td>{0}</td>".format(value))
            else:
                cells.append("<td>{0}</td>".format(escape("" if pd.isna(value) else str(value))))
        body_rows.append("<tr>{0}</tr>".format("".join(cells)))
    return "<div class='table-wrap'><table class='table'><thead><tr>{0}</tr></thead><tbody>{1}</tbody></table></div>".format(
        header,
        "".join(body_rows),
    )




def top_terms(series, limit=6):
    counts = {}
    for value in series.fillna(""):
        for part in str(value).split(","):
            text = part.strip()
            if not text or text.lower() == "nan":
                continue
            counts[text] = counts.get(text, 0) + 1
    ordered = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    return ordered[:limit]


STAT_TONES = {
    "antiSMASH": "violet",
    "GECCO": "green",
    "DeepBGC": "blue",
    "ARTS": "orange",
    "dbCAN CGC": "teal",
    "Total loci": "indigo",
    "Multi-caller loci": "lime",
    "MIBiG comparisons": "cyan",
    "ARTS known hits": "orange",
}

STAT_GLYPHS = {
    "antiSMASH": "○",
    "GECCO": "◇",
    "DeepBGC": "□",
    "ARTS": "⛨",
    "dbCAN CGC": "⌘",
    "Total loci": "◔",
    "Multi-caller loci": "◎",
    "MIBiG comparisons": "⬡",
    "ARTS known hits": "⛨",
}

def stat_card(label, value, note, tone=None, glyph=None):
    tone = tone or STAT_TONES.get(label, "violet")
    glyph = glyph or STAT_GLYPHS.get(label, "•")
    return (
        "<article class='metric-card metric-{tone}'>"
        "<div class='metric-icon' aria-hidden='true'>{glyph}</div>"
        "<div class='metric-copy'>"
        "<div class='metric-label'>{label}</div>"
        "<div class='metric-value'>{value}</div>"
        "<p>{note}</p>"
        "</div>"
        "</article>"
    ).format(
        label=escape(label),
        value=escape(str(value)),
        note=escape(note),
        tone=escape(tone),
        glyph=escape(glyph),
    )


def info_card(title, body, icon="i"):
    return (
        "<article class='info-card'>"
        "<div class='info-icon' aria-hidden='true'>{icon}</div>"
        "<div class='info-copy'>"
        "<h3>{title}</h3>"
        "<p>{body}</p>"
        "</div>"
        "</article>"
    ).format(icon=escape(icon), title=escape(title), body=escape(body))


def donut_style(pairs):
    total = sum(count for _, count in pairs)
    if total <= 0:
        return "background: #eef2ff;"
    cursor = 0.0
    stops = []
    for idx, (_, count) in enumerate(pairs):
        fraction = (float(count) / float(total)) * 100.0
        color = DONUT_COLORS[idx % len(DONUT_COLORS)]
        start = cursor
        end = cursor + fraction
        stops.append("{0} {1:.2f}% {2:.2f}%".format(color, start, end))
        cursor = end
    if cursor < 100:
        stops.append("#eef2ff {0:.2f}% 100%".format(cursor))
    return "background: conic-gradient({0});".format(", ".join(stops))


def render_donut_panel(title, pairs, total):
    if not pairs:
        return (
            "<section class='panel'>"
            "<div class='section-head'><h2>{title}</h2><span class='section-accent'></span></div>"
            "<p class='muted'>No dominant categories detected.</p>"
            "</section>"
        ).format(title=escape(title))
    legend = "".join(
        "<li><span class='legend-swatch' style='background:{color}'></span>"
        "<span class='legend-label'>{label}</span><span class='legend-value'>{count}</span></li>".format(
            color=DONUT_COLORS[idx % len(DONUT_COLORS)],
            label=escape(label),
            count=count,
        )
        for idx, (label, count) in enumerate(pairs)
    )
    shown_signals = sum(count for _, count in pairs)
    return (
        "<section class='panel chart-panel'>"
        "<div class='section-head'><h2>{title}</h2><span class='section-accent'></span></div>"
        "<div class='chart-grid'>"
        "<div class='donut-shell'>"
        "<div class='donut-chart' style='{style}'>"
        "<div class='donut-hole'><strong>{signals}</strong><span>Top signals</span></div>"
        "</div>"
        "</div>"
        "<ul class='legend-list'>{legend}</ul>"
        "</div>"
        "<p class='muted chart-note'>Six most frequent signals across {regions} grouped candidate loci. "
        "A region may have multiple signals or none.</p>"
        "</section>"
    ).format(
        title=escape(title), style=donut_style(pairs), signals=shown_signals,
        regions=escape(str(total)), legend=legend,
    )




def footer_strip(sample, generated_at):
    return (
        "<footer class='report-footer'>"
        "<span>BGC-XPLORER Report</span><span class='dot'></span>"
        "<span>Sample: {sample}</span><span class='dot'></span>"
        "<span>{generated_at}</span><span class='dot'></span>"
        "<span>Workflow 1.0</span>"
        "</footer>"
    ).format(sample=escape(sample), generated_at=escape(generated_at))








def gene_map_button(row, map_index):
    map_key = as_text(row.get("cluster_key"), "") or as_text(row.get("consensus_id"), "")
    item = map_index.get(str(map_key), {})
    map_path = as_text(item.get("map_path"), "")
    if not map_path:
        return "<span class='muted'>n/a</span>"
    interpretation = as_text(
        row.get("biological_interpretation"),
        "No biological interpretation available for this cluster.",
    )
    display_label = as_text(item.get("cluster_label"), map_key)
    ai_cluster_id = as_text(item.get("ai_cluster_id"), "")
    ai_enabled = "1" if ai_cluster_id else "0"
    report_dir = os.path.dirname(os.path.abspath(str(snakemake.output[0])))
    abs_path = os.path.abspath(os.path.join(report_dir, map_path))
    file_src = "file://{0}".format(quote(abs_path))
    return (
        "<button class='gene-map-btn' type='button' "
        "data-sample='{sample}' data-cluster='{cluster}' data-src='{src}' data-genes='{genes}' "
        "data-file-src='{file_src}' data-summary='{summary}' data-cluster-label='{cluster_label}' "
        "data-ai-cluster='{ai_cluster}' data-ai-enabled='{ai_enabled}' "
        "data-interpretation='{interpretation}'>Gene map</button>"
    ).format(
        sample=escape(str(sample), quote=True),
        cluster=escape(str(map_key), quote=True),
        src=escape(map_path, quote=True),
        file_src=escape(file_src, quote=True),
        genes=escape(str(item.get("gene_count", "")), quote=True),
        summary=escape(as_text(item.get("category_summary"), ""), quote=True),
        cluster_label=escape(display_label, quote=True),
        ai_cluster=escape(ai_cluster_id, quote=True),
        ai_enabled=escape(ai_enabled, quote=True),
        interpretation=escape(interpretation, quote=True),
    )


def load_gene_map_svg_data(cluster_maps):
    report_dir = os.path.dirname(os.path.abspath(str(snakemake.output[0])))
    svg_data = {}
    for _, row in cluster_maps.iterrows():
        cluster_key = as_text(row.get("cluster_key"), "") or as_text(row.get("consensus_id"), "")
        map_path = as_text(row.get("map_path"), "")
        if not cluster_key or not map_path:
            continue
        svg_path = os.path.abspath(os.path.join(report_dir, map_path))
        if not os.path.exists(svg_path):
            continue
        with open(svg_path, "r", encoding="utf-8") as handle:
            svg_data[cluster_key] = handle.read()
    return svg_data


def _cazyme_family(value):
    text = ("" if pd.isna(value) else str(value)).strip()
    if not text or text in {"-", "nan"}:
        return ""
    first = text.split("+")[0].split("|")[0].split(",")[0].split(";")[0].strip()
    family = first.split("(")[0].strip()
    return family if family and family != "-" else ""


def _best_dbcan_family(row):
    for column in ["dbcan_recommendation", "dbcan_hmm", "dbcan_subfamily", "dbcan_diamond"]:
        family = _cazyme_family(row.get(column))
        if family:
            return family
    return ""


def load_gene_table_data(cluster_genes_path, cluster_maps):
    """Build a {cluster_key: [gene_rows,...]} mapping for the popover gene table.

    Each gene row carries: locus, gene_name, product, eggnog, dbcan.
    The table is shared between consensus clusters (key = consensus_id) and
    dbCAN CGC clusters (key = {sample}_cgc_{cgc_id}); both keys appear in
    cluster_maps and in cluster_genes.tsv as `cluster_key`.
    """
    if not cluster_genes_path or not os.path.exists(cluster_genes_path):
        return {}
    df = pd.read_csv(cluster_genes_path, sep="\t")
    if df.empty or "cluster_key" not in df.columns:
        return {}
    keys = set()
    for _, row in cluster_maps.iterrows():
        key = as_text(row.get("cluster_key"), "") or as_text(row.get("consensus_id"), "")
        if key:
            keys.add(str(key))
    if not keys:
        return {}
    df = df[df["cluster_key"].astype(str).isin(keys)].copy()
    if df.empty:
        return {}
    table = {}
    for _, row in df.iterrows():
        key = str(as_text(row.get("cluster_key"), ""))
        if not key:
            continue
        gene_name = as_text(row.get("bakta_gene"), "")
        if not gene_name or gene_name.lower() == "nan":
            gene_name = as_text(row.get("preferred_name"), "")
        table.setdefault(key, []).append({
            "locus": as_text(row.get("locus_tag"), ""),
            "gene": gene_name,
            "product": as_text(row.get("bakta_product"), ""),
            "eggnog": as_text(row.get("eggnog_description"), ""),
            "dbcan": _best_dbcan_family(row),
        })
    return table


def dbcan_cluster_interpretation(row):
    substrate = as_text(row.get("predicted_substrate"), "unresolved")
    support = as_text(row.get("prediction_source"), "dbCAN")
    genes = as_int(row.get("genes"), 0)
    cazymes = as_int(row.get("cazyme_genes"), 0)
    signatures = as_int(row.get("signature_genes"), 0)
    return (
        "dbCAN carbohydrate gene cluster with predicted substrate {substrate}. "
        "Region contains {genes} genes, including {cazymes} CAZyme genes and {signatures} signature genes. "
        "Substrate support source: {support}."
    ).format(
        substrate=substrate,
        genes=genes,
        cazymes=cazymes,
        signatures=signatures,
        support=support,
    )






sample = snakemake.wildcards.sample
qc = read_json(snakemake.input.qc)
bakta = read_json(snakemake.input.bakta)
antismash = load_table_if_exists(snakemake.input.antismash, [
    "sample", "tool", "contig", "start", "end", "strand", "bgc_id", "bgc_type", "product", "score", "confidence", "source_file"
])
gecco = load_table_if_exists(snakemake.input.gecco, antismash.columns.tolist())
deepbgc = load_table_if_exists(snakemake.input.deepbgc, antismash.columns.tolist())
arts = load_table_if_exists(snakemake.input.arts, [
    "sample", "tool", "contig", "start", "end", "feature", "score", "evidence", "source_file"
])
dbcan = load_table_if_exists(snakemake.input.dbcan, [
    "sample", "cgcid", "cgc_id", "contig", "cluster_start", "cluster_end", "length_bp",
    "genes", "cazyme_genes", "tc_genes", "tf_genes", "stp_genes", "sulfatase_genes",
    "peptidase_genes", "signature_genes", "pul_id", "pul_substrate", "pul_bitscore",
    "signature_pairs", "dbcan_sub_substrate", "dbcan_sub_substrate_score",
    "predicted_substrate", "prediction_source", "source_file"
])
cluster_maps = load_table_if_exists(snakemake.input.cluster_maps, [
    "sample", "cluster_key", "cluster_label", "cluster_type", "ai_cluster_id",
    "consensus_id", "consensus_label", "contig", "start", "end", "length_bp",
    "gene_count", "category_summary", "map_path"
])
cluster_genes_path = str(snakemake.input.cluster_genes) if hasattr(snakemake.input, "cluster_genes") else ""
consensus = pd.read_csv(snakemake.input.consensus, sep="\t") if os.path.exists(snakemake.input.consensus) else pd.DataFrame()
evidence = pd.read_csv(snakemake.input.evidence, sep="\t") if os.path.exists(snakemake.input.evidence) else pd.DataFrame()
overlap = load_table_if_exists(snakemake.input.overlap, [
    "sample", "group_id", "tool_a", "bgc_id_a", "tool_b", "bgc_id_b",
    "contig", "overlap_bp", "overlap_relationship",
])
provenance = read_json(snakemake.input.provenance) if hasattr(snakemake.input, "provenance") else {}

map_index = {
    str(as_text(row.get("cluster_key"), "") or as_text(row.get("consensus_id"), "")): row.to_dict()
    for _, row in cluster_maps.iterrows()
    if as_text(row.get("cluster_key"), "") or as_text(row.get("consensus_id"), "")
}
gene_map_svg_data = load_gene_map_svg_data(cluster_maps)
gene_map_svg_json = json.dumps(gene_map_svg_data).replace("</", "<\\/")
gene_table_data = load_gene_table_data(cluster_genes_path, cluster_maps)
gene_table_json = json.dumps(gene_table_data).replace("</", "<\\/")
ai_config = snakemake.config.get("report", {}).get("ai", {})
ai_endpoint = ai_config.get("endpoint", "/analyze_cluster")

supported = consensus[consensus.get("support_count", pd.Series(dtype=int)).fillna(0).astype(int) > 1] if not consensus.empty else pd.DataFrame()

top_product_terms = top_terms(consensus["products"]) if "products" in consensus.columns else []
top_type_terms = top_terms(consensus["bgc_types"]) if "bgc_types" in consensus.columns else []
mibig_backed = evidence[evidence.get("best_mibig_id", pd.Series(dtype=str)).fillna("").astype(str).str.strip() != ""] if not evidence.empty else pd.DataFrame()
arts_known = evidence[pd.to_numeric(evidence.get("arts_known_hits", pd.Series(dtype=float)), errors="coerce").fillna(0) > 0] if not evidence.empty else pd.DataFrame()

hero_summary = (
    "Sample {sample} produced {consensus_count} grouped candidate loci."
).format(
    sample=sample,
    consensus_count=len(consensus),
)

generated_at = datetime.now().strftime("%b %d, %Y %H:%M")
glance_metrics = {
    "Total loci": (len(consensus), "candidate loci grouped across BGC predictors"),
    "Multi-caller loci": (len(supported), "loci containing predictions from at least two callers"),
    "MIBiG comparisons": (len(mibig_backed), "loci with a representative MIBiG comparison"),
    "ARTS known hits": (len(arts_known), "loci overlapping known-hit ARTS records"),
}
glance_cards = "".join(
    stat_card(label, glance_metrics[label][0], glance_metrics[label][1], tone="violet")
    for label in PRIMARY_GLANCE_LABELS
)

dbcan_table = pd.DataFrame()
if not dbcan.empty:
    dbcan_found = dbcan.copy()
    dbcan_found["predicted_substrate"] = dbcan_found["predicted_substrate"].fillna("").astype(str).str.strip()
    dbcan_found = dbcan_found[
        (dbcan_found["predicted_substrate"] != "")
        & (dbcan_found["predicted_substrate"].str.lower() != "unresolved")
    ].copy()
    dbcan_found["cluster_key"] = dbcan_found["cgc_id"].astype(str).map(
        lambda value: "{0}_cgc_{1}".format(sample, value)
    )
    dbcan_found["location"] = [
        location_html(row["contig"], row["cluster_start"], row["cluster_end"], row.get("length_bp"))
        for _, row in dbcan_found.iterrows()
    ]
    dbcan_table = dbcan_found[[
        "cluster_key",
        "cgc_id",
        "location",
        "genes",
        "cazyme_genes",
        "signature_genes",
        "predicted_substrate",
        "prediction_source",
        "dbcan_sub_substrate",
    ]].copy()
    dbcan_table = dbcan_table.rename(columns={
        "cgc_id": "CGC",
        "location": "Location",
        "genes": "Genes",
        "cazyme_genes": "CAZymes",
        "signature_genes": "Signatures",
        "predicted_substrate": "Predicted substrate",
        "prediction_source": "Support source",
        "dbcan_sub_substrate": "dbCAN-sub substrate",
    })
    if not dbcan_table.empty:
        dbcan_table["Gene map"] = [
            gene_map_button(
                {
                    "cluster_key": row["cluster_key"],
                    "biological_interpretation": dbcan_cluster_interpretation(row),
                },
                map_index,
            )
            for _, row in dbcan_found.iterrows()
        ]
        ordered = ["Gene map"] + [column for column in dbcan_table.columns if column != "Gene map" and column != "cluster_key"]
        dbcan_table = dbcan_table[ordered]

evidence_panel = render_evidence_table(
    evidence,
    map_index,
    gene_map_button,
    min_containment=float(snakemake.config.get("consensus", {}).get("min_containment", 0.80)),
    standalone=False,
    grouping_note=grouping_description(
        snakemake.config.get("consensus", {}).get("mode", "voting"),
        float(snakemake.config.get("consensus", {}).get("min_containment", 0.80)),
        int(snakemake.config.get("consensus", {}).get("min_callers_per_gene", 2)),
    ),
)
# The report shell already ships a tab strip (.table-panel / .tabs-nav / .tab-btn /
# .tab-pane) together with its click handler, so these panels reuse it instead of
# defining a second, competing tab style. Only the spacing inside a pane and the
# count pill are added here.
DOWNLOADS_CSS = """
<style>
.download-list { list-style: none; margin: 0; padding: 0; display: flex;
  flex-direction: column; gap: 10px; }
.download-item { display: flex; align-items: center; gap: 14px; padding: 15px 17px;
  border: 1px solid var(--line-strong); border-radius: 14px; background: var(--panel);
  color: inherit; text-decoration: none; transition: border-color var(--transition),
  box-shadow var(--transition); }
.download-item:hover, .download-item:focus-visible { border-color: var(--accent);
  box-shadow: var(--shadow-soft); }
.download-icon { flex: 0 0 auto; width: 38px; height: 38px; border-radius: 11px;
  display: grid; place-items: center; font-size: 19px; color: var(--accent);
  background: var(--accent-soft); }
.download-copy { display: flex; flex-direction: column; gap: 3px; min-width: 0; }
.download-name { font-weight: 700; color: var(--text-strong); font-size: 15px; }
.download-note { color: var(--muted); font-size: 13px; line-height: 1.45; }
.download-meta { margin-left: auto; text-align: right; color: var(--muted);
  font-size: 11.5px; font-family: var(--mono); white-space: nowrap; }
.download-missing { opacity: .55; pointer-events: none; }
</style>
"""


def human_size(num_bytes):
    for unit, step in (("MB", 1024 * 1024), ("KB", 1024)):
        if num_bytes >= step:
            return "{:.1f} {}".format(num_bytes / float(step), unit)
    return "{} B".format(num_bytes)


def downloads_panel(report_dir, sample):
    """Relative links to the tables behind the report.

    Bakta and eggNOG annotations are several MB together, so the files are
    linked rather than embedded; the paths resolve both from the served
    results tree and from the report opened directly off disk.
    """
    entries = [
        ("../summary/region_evidence.tsv", "BGC table",
         "One row per grouped candidate locus, as shown in the BGC tab"),
        ("../summary/dbcan.cgc.tsv", "CGC table",
         "All dbCAN carbohydrate gene clusters, including those without a resolved substrate"),
        ("../bakta/{0}.features.tsv".format(sample), "Bakta annotation",
         "Every annotated feature: locus tags, coordinates, products"),
        ("../summary/eggnog.annotations.tsv", "eggNOG annotation",
         "Orthology assignments with EC, KEGG KO and COG categories"),
    ]
    items = []
    available = 0
    for relative, title, note in entries:
        path = os.path.normpath(os.path.join(report_dir, relative))
        exists = os.path.isfile(path)
        if exists:
            available += 1
            meta = "{name} \u00b7 {size}".format(
                name=os.path.basename(path), size=human_size(os.path.getsize(path))
            )
        else:
            meta = "not produced"
        items.append(
            "<li><a class='download-item{missing}' href='{href}' download>"
            "<span class='download-icon' aria-hidden='true'>\u2913</span>"
            "<span class='download-copy'><span class='download-name'>{title}</span>"
            "<span class='download-note'>{note}</span></span>"
            "<span class='download-meta'>{meta}</span></a></li>".format(
                missing="" if exists else " download-missing",
                href=escape(relative, quote=True), title=escape(title),
                note=escape(note), meta=escape(meta),
            )
        )
    body = (
        DOWNLOADS_CSS
        + "<div class='downloads-panel'>"
        + "<ul class='download-list'>" + "".join(items) + "</ul>"
        + "</div>"
    )
    return body, available


CLUSTER_TABS_CSS = """
<style>
.cluster-tabs .tabs-nav { padding: 0 22px; gap: 4px; }
.cluster-tabs .tab-btn { font-size: 18px; font-weight: 700; letter-spacing: -0.01em;
  padding: 21px 28px 18px; gap: 12px; border-bottom-width: 4px; }
.cluster-tabs .tab-icon { font-size: 20px; }
.cluster-tabs .tab-pane > .evidence-panel,
.cluster-tabs .tab-pane > .cgc-panel,
.cluster-tabs .tab-pane > .downloads-panel { padding: 24px 24px 8px; }
.cluster-tabs .tab-count { font-size: 13px; font-weight: 700; line-height: 1;
  padding: 5px 11px; border-radius: 999px; color: var(--muted);
  background: rgba(94, 108, 152, 0.12); }
.cluster-tabs .tab-btn.is-active .tab-count { color: white; background: var(--accent); }
.cluster-tabs .tab-btn:hover { color: var(--accent-deep); }
@media (max-width: 720px) {
  .cluster-tabs .tabs-nav { padding: 0 12px; }
  .cluster-tabs .tab-btn { font-size: 16px; padding: 17px 18px 14px; gap: 9px; }
  .cluster-tabs .tab-pane > .evidence-panel,
  .cluster-tabs .tab-pane > .cgc-panel,
  .cluster-tabs .tab-pane > .downloads-panel { padding: 18px 16px 6px; }
}
</style>
"""


def cluster_tabs(panels):
    """Render the BGC and CGC tables as panes of the report's own tab strip."""
    buttons = []
    panes = []
    for index, (key, label, icon, count, body) in enumerate(panels):
        active = index == 0
        buttons.append(
            "<button class='tab-btn{active}' type='button' role='tab' "
            "aria-selected='{selected}' aria-controls='{key}' data-tab-target='{key}'>"
            "<span class='tab-icon' aria-hidden='true'>{icon}</span>"
            "<span>{label}</span><span class='tab-count'>{count}</span></button>".format(
                active=" is-active" if active else "", selected="true" if active else "false",
                key=escape(key, quote=True), icon=escape(icon), label=escape(label),
                count=escape(str(count)),
            )
        )
        panes.append(
            "<div class='tab-pane{active}' id='{key}' role='tabpanel'{hidden}>{body}</div>".format(
                active=" is-active" if active else "", key=escape(key, quote=True),
                hidden="" if active else " hidden", body=body,
            )
        )
    return (
        CLUSTER_TABS_CSS
        + "<section class='panel table-panel cluster-tabs'>"
        + "<div class='tabs-nav' role='tablist' aria-label='Cluster tables'>"
        + "".join(buttons) + "</div>" + "".join(panes) + "</section>"
    )


cgc_panel = ""
if not dbcan_table.empty:
    cgc_panel = (
        "<div class='cgc-panel'>"
        "<dl class='column-legend'>"
        "<div><dt>Gene map</dt><dd>opens the annotated gene diagram</dd></div>"
        "<div><dt>CGC</dt><dd>dbCAN carbohydrate gene cluster id</dd></div>"
        "<div><dt>Location</dt><dd>coordinates, contig and length</dd></div>"
        "<div><dt>Genes</dt><dd>annotated genes inside the cluster</dd></div>"
        "<div><dt>CAZymes</dt><dd>carbohydrate-active enzyme genes</dd></div>"
        "<div><dt>Signatures</dt><dd>signature genes: CAZyme, transporter, "
        "transcription factor, sulfatase or peptidase</dd></div>"
        "<div><dt>Predicted substrate</dt><dd>the carbohydrate the cluster likely acts on</dd></div>"
        "<div><dt>Support source</dt><dd>dbCAN-PUL homology, dbCAN-sub majority vote, or both</dd></div>"
        "<div><dt>dbCAN-sub substrate</dt><dd>the dbCAN-sub call on its own</dd></div>"
        "</dl>"
        "{table}</div>"
    ).format(table=render_html_table(dbcan_table, html_columns={"Gene map", "Location"}))

download_panel, download_count = downloads_panel(
    os.path.dirname(os.path.abspath(str(snakemake.output[0]))), sample
)

sections = [
    report_home_link(),
    (
        "<section class='hero-card'>"
        "<div class='hero-brand'>"
        + ("<img class='brand-logo' src='{src}' alt='BGC-XPLORER' />".format(src=LOGO_DATA_URI) if LOGO_DATA_URI else "<div class='brand-mark' aria-hidden='true'>⌬</div>")
        + "</div>"
        "<div class='hero-panel'>"
        "<h1>{title}</h1>"
        "<p>{summary}</p>"
        "</div>"
        "<div class='meta-badge'>"
        "<div class='meta-badge-label'>Report generated</div>"
        "<div class='meta-badge-value'>{generated}</div>"
        "<div class='meta-badge-foot'>University of Plovdiv</div>"
        "</div>"
        "</section>"
    ).format(
        title=escape("sample {0}".format(sample)),
        summary=escape(hero_summary),
        generated=escape(generated_at),
    ),
    "<section class='glance-shell'><div class='section-head'><h2>At a Glance</h2><span class='section-accent'></span></div><div class='metrics-grid'>{cards}</div></section>".format(
        cards=glance_cards
    ),
    (
        "<section class='info-grid summary-only'>"
        "{summary_card}"
        "</section>"
    ).format(
        summary_card=info_card(
            "Summary",
            "The table includes all {0} grouped candidate loci. {1} have predictions "
            "from at least two tools, {2} have a representative MIBiG comparison, "
            "and {3} overlap a known-hit ARTS record. Use the filters and column "
            "headers to explore each signal separately.".format(
                len(evidence), len(supported), len(mibig_backed), len(arts_known)
            ),
            icon="▣",
        ),
    ),
    "<section class='chart-panels'>{left}{right}</section>".format(
        left=render_donut_panel("Top Product Signals", top_product_terms, len(consensus)),
        right=render_donut_panel("Top BGC Class Signals", top_type_terms, len(consensus)),
    ),
    (
        "<section class='panel gene-map-viewer' id='gene-map-viewer' "
        "data-ai-endpoint='{endpoint}' hidden>"
        "<div class='gene-map-head'>"
        "<div><h2 id='gene-map-title'>Cluster Gene Map</h2>"
        "<p id='gene-map-meta' class='muted'></p></div>"
        "<button class='gene-map-close' type='button' aria-label='Close gene map'>Close</button>"
        "</div>"
        "<div class='gene-map-stage'>"
        "<div id='gene-map-inline' class='gene-map-inline'></div>"
        "<img id='gene-map-image' alt='' hidden>"
        "<div id='gene-map-interpretation' class='gene-map-interpretation'>"
        "<h3>Biological interpretation</h3>"
        "<p></p>"
        "</div>"
        "<div id='gene-map-gene-table' class='gene-map-gene-table' hidden>"
        "<h3>Genes in cluster</h3>"
        "<div class='gene-table-wrap'><table class='gene-table'>"
        "<thead><tr><th>Locus</th><th>Gene</th><th>Product</th>"
        "<th>eggNOG description</th><th>dbCAN family</th></tr></thead>"
        "<tbody></tbody></table></div>"
        "</div>"
        "<div class='gene-map-ai'>"
        "<div class='gene-map-ai-head'>"
        "<button id='gene-map-ai-button' class='ai-analysis-btn' type='button'>AI analysis</button>"
        "<span id='gene-map-ai-spinner' class='ai-spinner' hidden></span>"
        "<span id='gene-map-ai-status' class='muted'></span>"
        "</div>"
        "<p class='ai-analysis-note'>AI model outputs are computational interpretations "
        "that require expert review and experimental validation.</p>"
        "<div id='gene-map-ai-result' class='ai-analysis-result'></div>"
        "</div>"
        "</div>"
        "<p id='gene-map-error' class='gene-map-error' hidden></p>"
        "<div class='callout'><p id='gene-map-summary'></p></div>"
        "</section>"
    ).format(endpoint=escape(ai_endpoint, quote=True)),
    "<script type='application/json' id='gene-map-svg-data'>{0}</script>".format(gene_map_svg_json),
    "<script type='application/json' id='gene-map-table-data'>{0}</script>".format(gene_table_json),
    cluster_tabs(
        [("cluster-tab-bgc", "BGC", "\u25ce", len(evidence), evidence_panel)]
        + ([("cluster-tab-cgc", "CGC", "\u2318", len(dbcan_table), cgc_panel)] if cgc_panel else [])
        + [("cluster-tab-downloads", "DOWNLOADS", "\u2913", download_count, download_panel)]
    ),
]
if provenance:
    sections.append(reproducibility_panel(provenance))
sections.append(footer_strip(sample, generated_at))

html = html_page("BGC-XPLORER : sample {0}".format(sample), sections)
with open(snakemake.output[0], "w", encoding="utf-8") as handle:
    handle.write(html)
