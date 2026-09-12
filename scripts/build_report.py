import base64
import json
import os
from datetime import datetime
from html import escape
from urllib.parse import quote

import pandas as pd

from common import df_to_html_table, format_consensus_label, html_page, load_table_if_exists, read_json


def _load_logo_base64():
    """Load the project logo as a base64 PNG data URI for embedding in the report.

    Prefer the higher-resolution logo.png in the report assets directory
    (next to cluster_maps) and fall back to logo@56.png when only the smaller
    raster is available. Returns an empty string when the logo is not available
    so the report keeps rendering without it.
    """
    report_dir = os.path.dirname(os.path.abspath(str(snakemake.output[0])))
    candidates = [
        os.path.join(report_dir, "assets", "logo.png"),
        os.path.join(report_dir, "assets", "logo@56.png"),
    ]
    for path in candidates:
        if os.path.exists(path):
            with open(path, "rb") as handle:
                data = base64.b64encode(handle.read()).decode("ascii")
            return "data:image/png;base64,{0}".format(data)
    return ""


LOGO_DATA_URI = _load_logo_base64()


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


def derep_source_label(value):
    text = as_text(value, "").strip().lower()
    mapping = {
        "comparippson_html": "CompariPPson",
        "knownclusterblast": "KnownClusterBlast",
        "clustercompare": "ClusterCompare",
    }
    return mapping.get(text, text or "n/a")


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
    "Consensus": "indigo",
    "Multi-tool": "lime",
    "High-confidence": "sky",
    "High-interest": "amber",
    "MIBiG hits": "cyan",
}

STAT_GLYPHS = {
    "antiSMASH": "○",
    "GECCO": "◇",
    "DeepBGC": "□",
    "ARTS": "⛨",
    "dbCAN CGC": "⌘",
    "Consensus": "◔",
    "Multi-tool": "◎",
    "High-confidence": "✦",
    "High-interest": "✧",
    "MIBiG hits": "◌",
}

DONUT_COLORS = ["#635bff", "#73c95b", "#4f8dfd", "#f59b38", "#41bcc7", "#f176c5", "#7d6bff", "#86d992"]


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
    return (
        "<section class='panel chart-panel'>"
        "<div class='section-head'><h2>{title}</h2><span class='section-accent'></span></div>"
        "<div class='chart-grid'>"
        "<div class='donut-shell'>"
        "<div class='donut-chart' style='{style}'>"
        "<div class='donut-hole'><strong>{total}</strong><span>Total</span></div>"
        "</div>"
        "</div>"
        "<ul class='legend-list'>{legend}</ul>"
        "</div>"
        "</section>"
    ).format(title=escape(title), style=donut_style(pairs), total=escape(str(total)), legend=legend)


def table_tabs(items):
    nav = []
    panes = []
    for index, item in enumerate(items):
        slug = "tab-{0}".format(index)
        nav.append(
            "<button class='tab-btn{active}' type='button' data-tab-target='{slug}' aria-selected='{selected}'>"
            "<span class='tab-icon'>{icon}</span><span>{label}</span></button>".format(
                active=" is-active" if index == 0 else "",
                slug=slug,
                selected="true" if index == 0 else "false",
                icon=escape(item.get("icon", "•")),
                label=escape(item["label"]),
            )
        )
        panes.append(
            "<div class='tab-pane{active}' id='{slug}'{hidden}>{content}</div>".format(
                active=" is-active" if index == 0 else "",
                slug=slug,
                hidden="" if index == 0 else " hidden",
                content=item["content"],
            )
        )
    return (
        "<section class='panel table-panel'>"
        "<div class='tabs-nav' role='tablist'>{nav}</div>"
        "<div class='tabs-body'>{panes}</div>"
        "</section>"
    ).format(nav="".join(nav), panes="".join(panes))


def footer_strip(sample, generated_at):
    return (
        "<footer class='report-footer'>"
        "<span>BGC-XPLORER Report</span><span class='dot'></span>"
        "<span>Sample: {sample}</span><span class='dot'></span>"
        "<span>{generated_at}</span><span class='dot'></span>"
        "<span>Workflow 1.0</span>"
        "</footer>"
    ).format(sample=escape(sample), generated_at=escape(generated_at))


def chip(text, kind="accent"):
    return "<span class='chip chip-{kind}'>{text}</span>".format(
        kind=escape(kind), text=escape(str(text))
    )


def build_why_prioritized(row):
    existing = str(row.get("why_prioritized", "")).strip()
    if existing and existing.lower() != "nan":
        return existing
    support_count = as_int(row.get("support_count"))
    arts_hits = as_int(row.get("arts_hits"))
    support_tools = [part.strip() for part in as_text(row.get("support_tools"), "").split(",") if part.strip()]
    notes = str(row.get("notes", "")).strip()
    reasons = []
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

    if support_tools:
        reasons.append("callers: {0}".format(", ".join(support_tools)))

    if support_count <= 1:
        reasons.append("lower confidence")
    elif support_count >= 3 and arts_hits > 0:
        reasons.append("strong follow-up candidate")
    elif support_count >= 3:
        reasons.append("high cross-tool confidence")
    else:
        reasons.append("moderate cross-tool confidence")

    if notes and notes.lower() != "nan":
        compact_notes = notes.replace("ARTS overlap: ", "")
        reasons.append(compact_notes)

    return " + ".join(reasons)


def cluster_card(row):
    location = "{contig}:{start:,}-{end:,}".format(
        contig=as_text(row.get("contig")),
        start=as_int(row.get("start")),
        end=as_int(row.get("end")),
    )
    support_tools = [part.strip() for part in as_text(row.get("support_tools"), "").split(",") if part.strip()]
    support_label = ", ".join(support_tools) if support_tools else "single-tool"
    products = as_text(row.get("products"), "not annotated")
    bgc_types = as_text(row.get("bgc_types"), "not classified")
    interpretation = as_text(row.get("biological_interpretation"), "No interpretation text available.")
    notes = str(row.get("notes", "")).strip()
    note_block = ""
    if notes and notes.lower() != "nan":
        note_block = "<div class='callout'><p>{notes}</p></div>".format(notes=escape(notes))
    return (
        "<article class='cluster-card'>"
        "<h3>{cluster_id}</h3>"
        "<p>{location}</p>"
        "<div class='chip-row'>"
        "{support_chip}{tools_chip}{arts_chip}{score_chip}"
        "</div>"
        "<div class='kv'>"
        "<div class='k'>Products</div><div class='v'>{products}</div>"
        "<div class='k'>BGC types</div><div class='v'>{bgc_types}</div>"
        "<div class='k'>Meaning</div><div class='v'>{interpretation}</div>"
        "</div>"
        "{note_block}"
        "</article>"
    ).format(
        cluster_id=escape(as_text(row.get("consensus_id"))),
        location=escape(location),
        support_chip=chip("{0} tools".format(as_int(row.get("support_count"))), "accent"),
        tools_chip=chip(support_label, "blue"),
        arts_chip=chip("ARTS hits: {0}".format(as_int(row.get("arts_hits"))), "green"),
        score_chip=chip("Priority: {0}".format(as_text(row.get("priority_score"), "0")), "accent"),
        products=escape(products),
        bgc_types=escape(bgc_types),
        interpretation=escape(interpretation),
        note_block=note_block,
    )


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


def df_to_cluster_table(df, hidden, map_index, limit=None):
    if df.empty:
        return "<p>No records found.</p>"
    source = df.head(limit).copy() if limit else df.copy()
    if "consensus_id" in source.columns:
        if "consensus_label" in source.columns:
            source["CONSENSUS"] = source["consensus_label"].astype(str)
        else:
            source["CONSENSUS"] = source["consensus_id"].astype(str).map(format_consensus_label)
        source["GENES"] = source["consensus_id"].astype(str).map(
            lambda value: as_int(map_index.get(value, {}).get("gene_count"), "")
        )
        ordered = []
        for column in source.columns:
            if column in {"CONSENSUS", "GENES", "consensus_label"}:
                continue
            ordered.append(column)
            if column == "contig":
                ordered.append("CONSENSUS")
            if column == "length_bp":
                ordered.append("GENES")
        source = source[[column for column in ordered if column in source.columns]]
    buttons = [
        gene_map_button(row, map_index)
        for _, row in source.iterrows()
    ]
    display = source.drop(columns=[col for col in hidden if col in source.columns], errors="ignore")
    display.insert(0, "Gene map", buttons)
    header = "".join(
        "<th>{label}</th>".format(label=escape(str(column)))
        for column in display.columns
    )
    body_rows = []
    for _, row in display.iterrows():
        cells = []
        for column in display.columns:
            value = row[column]
            if column == "Gene map":
                cells.append("<td>{0}</td>".format(value))
            else:
                cells.append("<td>{0}</td>".format(escape(as_text(value, ""))))
        body_rows.append("<tr>{0}</tr>".format("".join(cells)))
    return (
        "<div class='table-wrap'><table class='table'>"
        "<thead><tr>{header}</tr></thead><tbody>{rows}</tbody>"
        "</table></div>"
    ).format(header=header, rows="".join(body_rows))


def df_to_mibig_table(df, hidden, map_index, limit=None):
    """MIBiG-specific cluster table: inserts MIBiG product + Dereplication
    status columns right after GENES, and hides the other MIBiG detail fields.
    """
    if df.empty:
        return "<p>No records found.</p>"
    source = df.head(limit).copy() if limit else df.copy()
    if "consensus_id" in source.columns:
        if "consensus_label" in source.columns:
            source["CONSENSUS"] = source["consensus_label"].astype(str)
        else:
            source["CONSENSUS"] = source["consensus_id"].astype(str).map(format_consensus_label)
        source["GENES"] = source["consensus_id"].astype(str).map(
            lambda value: as_int(map_index.get(value, {}).get("gene_count"), "")
        )
        ordered = []
        for column in source.columns:
            if column in {"CONSENSUS", "GENES", "consensus_label"}:
                continue
            ordered.append(column)
            if column == "contig":
                ordered.append("CONSENSUS")
            if column == "length_bp":
                ordered.append("GENES")
        source = source[[column for column in ordered if column in source.columns]]
    buttons = [
        gene_map_button(row, map_index)
        for _, row in source.iterrows()
    ]
    display = source.drop(columns=[col for col in hidden if col in source.columns], errors="ignore")
    display.insert(0, "Gene map", buttons)

    # Reorder: move MIBiG product and Dereplication status right after GENES
    movable = ["MIBiG product", "Dereplication status"]
    cols = [c for c in display.columns if c not in movable]
    insert_at = cols.index("GENES") + 1 if "GENES" in cols else len(cols)
    for offset, column in enumerate(movable):
        if column in display.columns:
            cols.insert(insert_at + offset, column)
    display = display[cols]

    header = "".join(
        "<th>{label}</th>".format(label=escape(str(column)))
        for column in display.columns
    )
    body_rows = []
    for _, row in display.iterrows():
        cells = []
        for column in display.columns:
            value = row[column]
            if column == "Gene map":
                cells.append("<td>{0}</td>".format(value))
            else:
                cells.append("<td>{0}</td>".format(escape(as_text(value, ""))))
        body_rows.append("<tr>{0}</tr>".format("".join(cells)))
    return (
        "<div class='table-wrap'><table class='table'>"
        "<thead><tr>{header}</tr></thead><tbody>{rows}</tbody>"
        "</table></div>"
    ).format(header=header, rows="".join(body_rows))


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
prioritized = pd.read_csv(snakemake.input.prioritized, sep="\t") if os.path.exists(snakemake.input.prioritized) else pd.DataFrame()
overlap = pd.read_csv(snakemake.input.overlap, sep="\t") if os.path.exists(snakemake.input.overlap) else pd.DataFrame()

prioritized_table = prioritized.copy()
if not prioritized_table.empty:
    prioritized_table["why_prioritized"] = prioritized_table.apply(build_why_prioritized, axis=1)
    if "evidence_source" in prioritized_table.columns:
        prioritized_table["dereplication_source"] = prioritized_table["evidence_source"].apply(derep_source_label)
    else:
        prioritized_table["dereplication_source"] = ""

hidden_columns = [
    "consensus_id",
    "contig_id",
    "supporting_tools",
    "support_count",
    "num_supporting_tools",
    "candidate_ids",
    "bgc_types",
    "product_annotations",
    "biological_interpretation",
    "overlap_relationship",
    "core_gene_support",
    "contig_length",
    "distance_to_left_edge",
    "distance_to_right_edge",
    "edge_truncated",
    "possible_partial_BGC",
    "boundary_confidence",
    "priority_score",
    "priority_class",
    "confidence_category",
    "interest_category",
    "arts_hits",
    "notes",
    "why_prioritized",
    "why_not_prioritized",
    "recommended_followup",
    "best_mibig_id",
    "best_mibig_product",
    "best_mibig_class",
    "mibig_similarity",
    "dereplication_status",
    "novelty_score",
    "evidence_source",
    "dereplication_source",
]
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

tool_counts = [
    ("antiSMASH", len(antismash), "rule-based cluster calls and product classes"),
    ("GECCO", len(gecco), "machine learning BGC candidates"),
    ("DeepBGC", len(deepbgc), "domain-driven BGC calls and activity hints"),
    ("ARTS", len(arts), "resistance-linked genomic evidence"),
    ("dbCAN CGC", len(dbcan), "carbohydrate gene clusters with substrate calls"),
    ("Consensus", len(consensus), "merged loci across predictors"),
]

supported = consensus[consensus.get("support_count", pd.Series(dtype=int)).fillna(0).astype(int) > 1] if not consensus.empty else pd.DataFrame()
triple_supported = consensus[consensus.get("support_count", pd.Series(dtype=int)).fillna(0).astype(int) >= 3] if not consensus.empty else pd.DataFrame()

top_product_terms = top_terms(consensus["products"]) if "products" in consensus.columns else []
top_type_terms = top_terms(consensus["bgc_types"]) if "bgc_types" in consensus.columns else []
high_confidence = prioritized[prioritized.get("confidence_category", pd.Series(dtype=str)).astype(str) == "high-confidence BGC"] if not prioritized.empty else pd.DataFrame()
high_interest = prioritized[prioritized.get("interest_category", pd.Series(dtype=str)).astype(str) == "high-interest / potentially novel"] if not prioritized.empty else pd.DataFrame()
partial_regions = prioritized[prioritized.get("possible_partial_BGC", pd.Series(dtype=bool)).fillna(False).astype(bool)] if not prioritized.empty else pd.DataFrame()
mibig_backed = prioritized[prioritized.get("best_mibig_id", pd.Series(dtype=str)).fillna("").astype(str).str.strip() != ""] if not prioritized.empty else pd.DataFrame()
arts_supported = prioritized[pd.to_numeric(prioritized.get("arts_hits", pd.Series(dtype=float)), errors="coerce").fillna(0) > 0] if not prioritized.empty else pd.DataFrame()
featured_ids = set()
if not high_confidence.empty:
    featured_ids.update(high_confidence["consensus_id"].astype(str).tolist())
if not high_interest.empty:
    featured_ids.update(high_interest["consensus_id"].astype(str).tolist())
if not mibig_backed.empty:
    featured_ids.update(mibig_backed["consensus_id"].astype(str).tolist())
remaining_prioritized = prioritized_table[
    ~prioritized_table["consensus_id"].astype(str).isin(featured_ids)
] if not prioritized_table.empty else prioritized_table

top_cards = ""
if not prioritized.empty:
    top_cards = "".join(cluster_card(row) for _, row in prioritized.head(8).iterrows())

hero_summary = (
    "Sample {sample} produced {consensus_count} consensus BGC regions."
).format(
    sample=sample,
    consensus_count=len(consensus),
)

generated_at = datetime.now().strftime("%b %d, %Y %H:%M")
glance_cards = "".join(
    stat_card(label, value, note) for label, value, note in tool_counts
) + "".join([
    stat_card("Multi-tool", len(supported), "regions supported by at least two callers"),
    stat_card("High-confidence", len(high_confidence), "strongly supported consensus BGCs"),
    stat_card("High-interest", len(high_interest), "potentially novel or ARTS-rich candidates"),
    stat_card("MIBiG hits", len(mibig_backed), "regions with explicit dereplication evidence"),
])

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
    dbcan_table = dbcan_found[[
        "cluster_key",
        "contig",
        "cgc_id",
        "cluster_start",
        "cluster_end",
        "length_bp",
        "genes",
        "cazyme_genes",
        "signature_genes",
        "predicted_substrate",
        "prediction_source",
        "dbcan_sub_substrate",
    ]].copy()
    dbcan_table = dbcan_table.rename(columns={
        "contig": "Contig",
        "cgc_id": "CGC",
        "cluster_start": "Start",
        "cluster_end": "End",
        "length_bp": "Length (bp)",
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

table_items = []
if not high_confidence.empty:
    table_items.append({
        "label": "High-confidence BGCs",
        "icon": "✦",
        "content": df_to_cluster_table(high_confidence, hidden_columns, map_index, limit=15),
    })
if not high_interest.empty:
    table_items.append({
        "label": "High-interest / Potentially Novel BGCs",
        "icon": "✧",
        "content": df_to_cluster_table(high_interest, hidden_columns, map_index, limit=15),
    })
if not mibig_backed.empty:
    mibig_hidden = [c for c in hidden_columns if c not in {
        "best_mibig_product", "dereplication_status",
    }]
    mibig_display = mibig_backed.copy()
    rename_map = {
        "best_mibig_id": "MIBiG ID",
        "best_mibig_product": "MIBiG product",
        "best_mibig_class": "MIBiG class",
        "mibig_similarity": "Similarity (%)",
        "dereplication_status": "Dereplication status",
        "novelty_score": "Novelty score",
        "evidence_source": "Evidence source",
    }
    for old, new in rename_map.items():
        if old in mibig_display.columns:
            mibig_display = mibig_display.rename(columns={old: new})
    mibig_hidden = [rename_map.get(c, c) for c in mibig_hidden]
    table_items.append({
        "label": "MIBiG hits",
        "icon": "⬡",
        "content": (
            "<p class='muted table-note'>Consensus regions with explicit MIBiG dereplication evidence.</p>"
            + df_to_mibig_table(mibig_display, mibig_hidden, map_index)
        ),
    })
table_items.append({
    "label": "Other BGCs",
    "icon": "≣",
    "content": (
        "<p class='muted table-note'>This table shows the remaining ranked loci after removing entries already shown above.</p>"
        + df_to_cluster_table(remaining_prioritized, hidden_columns, map_index)
    ),
})
if not dbcan_table.empty:
    table_items.append({
        "label": "CGC substrate clusters",
        "icon": "⌘",
        "content": (
            "<p class='muted table-note'>dbCAN CGC substrate calls with resolved substrate predictions.</p>"
            + render_html_table(dbcan_table.head(50), html_columns={"Gene map"})
        ),
    })

sections = [
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
    "<section class='glance-shell'><div class='section-head'><h2>At a Glance<span class='think-dots' aria-hidden='true'><span class='dot'></span><span class='dot'></span><span class='dot'></span></span></h2><span class='section-accent'></span></div><div class='metrics-grid'>{cards}</div></section>".format(
        cards=glance_cards
    ),
    (
        "<section class='info-grid'>"
        "{summary_card}{howto_card}"
        "</section>"
    ).format(
        summary_card=info_card(
            "Executive Summary",
            "We currently separate regions into strong consensus BGCs and high-interest candidates. "
            "{0} regions are tagged as high-confidence, {1} as high-interest / potentially novel, "
            "and {2} have explicit antiSMASH-to-MIBiG dereplication evidence.".format(
                len(high_confidence), len(high_interest), len(mibig_backed)
            ),
            icon="▣",
        ),
        howto_card=info_card(
            "How to Read This",
            "Consensus regions are merged genomic intervals from antiSMASH, GECCO, and DeepBGC. They represent loci, not individual tool rows.",
            icon="?",
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
        "<div id='gene-map-ai-result' class='ai-analysis-result'></div>"
        "</div>"
        "</div>"
        "<p id='gene-map-error' class='gene-map-error' hidden></p>"
        "<div class='callout'><p id='gene-map-summary'></p></div>"
        "</section>"
    ).format(endpoint=escape(ai_endpoint, quote=True)),
    "<script type='application/json' id='gene-map-svg-data'>{0}</script>".format(gene_map_svg_json),
    "<script type='application/json' id='gene-map-table-data'>{0}</script>".format(gene_table_json),
    table_tabs(table_items),
]
sections.append(footer_strip(sample, generated_at))

html = html_page("BGC-XPLORER : sample {0}".format(sample), sections)
with open(snakemake.output[0], "w", encoding="utf-8") as handle:
    handle.write(html)
