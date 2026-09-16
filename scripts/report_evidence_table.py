"""One searchable table of observable BGC evidence for the HTML report."""

from html import escape

import pandas as pd


def clean(value):
    if pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"nan", "none"} else text


def number(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def mibig_method(value):
    return {
        "knownclusterblast": "KnownClusterBlast",
        "clustercompare_mibig": "ClusterCompare MIBiG",
        "comparippson_html": "CompariPPson",
    }.get(clean(value), clean(value))


def mibig_cell(row):
    accession = clean(row.get("best_mibig_id"))
    if not accession:
        return "<span class='muted'>—</span>"
    details = []
    for label, key in (
        ("Product", "best_mibig_product"), ("Method", "evidence_source"),
        ("Score metric", "score_metric"), ("Match score", "match_score"),
        ("Peptide similarity (%)", "mibig_similarity"),
        ("Matched genes", "matched_genes"), ("Core-gene hits", "core_gene_hits"),
    ):
        value = mibig_method(row.get(key)) if key == "evidence_source" else clean(row.get(key))
        if value:
            details.append("<div><strong>{}</strong>: {}</div>".format(escape(label), escape(value)))
    return (
        "<details class='mibig-detail'><summary>{accession}</summary>"
        "<div class='mibig-detail-body'>{details}</div></details>"
    ).format(accession=escape(accession), details="".join(details))


def caller_markers(row):
    callers = [part.strip().lower() for part in clean(row.get("support_tools")).split(",")]
    names = (("antismash", "antiSMASH"), ("gecco", "GECCO"), ("deepbgc", "DeepBGC"))
    return "".join(
        "<span class='caller-marker'>{}</span>".format(escape(label))
        for key, label in names if key in callers
    ) or "<span class='muted'>—</span>"


EVIDENCE_CSS = """
<style>
.evidence-panel { overflow: hidden; }
.evidence-description { margin-bottom: 16px; }
.evidence-toolbar { display: flex; flex-wrap: wrap; gap: 12px; align-items: end; margin-bottom: 12px; }
.evidence-toolbar label { display: flex; flex-direction: column; gap: 5px; color: var(--muted); font-size: 12px; font-weight: 600; }
.evidence-toolbar label.evidence-check { flex-direction: row; align-items: center; gap: 7px; padding: 8px 2px; font-size: 13px; }
.evidence-toolbar input[type='search'], .evidence-toolbar select { font: inherit; font-size: 14px; color: var(--text-strong); background: white; border: 1px solid var(--line-strong); border-radius: 9px; padding: 9px 10px; min-height: 39px; }
.evidence-toolbar input[type='search'] { width: min(290px, 80vw); }
.evidence-toolbar input[type='checkbox'] { accent-color: var(--accent); }
.evidence-count { margin: 0 0 12px; color: var(--muted); font-size: 13px; }
.evidence-panel .table-wrap { border: 1px solid var(--line); border-radius: 12px; }
.evidence-table { min-width: 1420px; }
.evidence-table th, .evidence-table td { padding: 11px 12px; vertical-align: middle; }
.evidence-table tbody tr:hover td { background: var(--accent-soft); }
.evidence-table th { white-space: nowrap; }
.sort-header { appearance: none; background: none; border: 0; color: inherit; font: inherit; font-weight: 700; cursor: pointer; padding: 0; text-align: left; }
.sort-header:hover, .sort-header:focus-visible { color: var(--accent); }
.sort-header::after { content: ' ↕'; color: var(--muted-2); font-weight: 400; }
th[aria-sort='ascending'] .sort-header::after { content: ' ↑'; color: var(--accent); }
th[aria-sort='descending'] .sort-header::after { content: ' ↓'; color: var(--accent); }
.caller-markers { display: flex; gap: 4px; flex-wrap: wrap; min-width: 112px; }
.caller-marker { display: inline-block; border: 1px solid rgba(95,87,255,.16); background: var(--accent-soft); color: var(--accent-deep); border-radius: 7px; padding: 3px 6px; font-size: 11px; font-weight: 600; white-space: nowrap; }
.evidence-table .region-name { font-weight: 700; white-space: nowrap; }
.evidence-table .signal-cell { max-width: 230px; min-width: 170px; font-size: 12px; line-height: 1.45; }
.mibig-detail summary { cursor: pointer; color: var(--accent-deep); font-weight: 700; white-space: nowrap; }
.mibig-detail-body { min-width: 230px; font-size: 12px; line-height: 1.55; padding: 8px 0 0; }
.evidence-table tr[hidden] { display: none; }
</style>
"""


EVIDENCE_JS = """
<script>
(function () {
  const table = document.getElementById('cluster-evidence-table');
  if (!table) return;
  const body = table.tBodies[0];
  const rows = Array.from(body.rows);
  const search = document.getElementById('evidence-search');
  const toolCount = document.getElementById('evidence-tool-count');
  const caller = document.getElementById('evidence-caller');
  const mibig = document.getElementById('evidence-mibig');
  const known = document.getElementById('evidence-known');
  const duf = document.getElementById('evidence-duf');
  const count = document.getElementById('evidence-count');
  function filterRows() {
    const query = search.value.trim().toLowerCase();
    let shown = 0;
    rows.forEach(function (row) {
      const visible = (!query || row.dataset.search.includes(query))
        && (!toolCount.value || Number(row.dataset.tools) === Number(toolCount.value))
        && (!caller.value || row.dataset.callers.split(',').includes(caller.value))
        && (!mibig.checked || row.dataset.mibig === '1')
        && (!known.checked || Number(row.dataset.known) > 0)
        && (!duf.checked || Number(row.dataset.duf) > 0);
      row.hidden = !visible;
      if (visible) shown += 1;
    });
    count.textContent = shown + ' of ' + rows.length + ' regions shown';
  }
  [search, toolCount, caller, mibig, known, duf].forEach(function (control) {
    control.addEventListener(control === search ? 'input' : 'change', filterRows);
  });
  table.querySelectorAll('.sort-header').forEach(function (button) {
    button.addEventListener('click', function () {
      const column = Number(button.dataset.column);
      const numeric = button.dataset.type === 'number';
      const descending = button.parentElement.getAttribute('aria-sort') === 'ascending';
      table.querySelectorAll('th[aria-sort]').forEach(function (cell) { cell.removeAttribute('aria-sort'); });
      button.parentElement.setAttribute('aria-sort', descending ? 'descending' : 'ascending');
      const sorted = rows.slice().sort(function (a, b) {
        const left = a.cells[column].dataset.sort || '';
        const right = b.cells[column].dataset.sort || '';
        const comparison = numeric ? Number(left) - Number(right)
          : left.localeCompare(right, undefined, {numeric: true, sensitivity: 'base'});
        return descending ? -comparison : comparison;
      });
      sorted.forEach(function (row) { body.appendChild(row); });
    });
  });
  filterRows();
})();
</script>
"""


def render_evidence_table(frame, map_index, gene_map_button, overlap_fraction=0.30):
    columns = (
        ("Gene map", False), ("Region", False), ("Contig", False),
        ("Start", True), ("End", True), ("Length (bp)", True),
        ("Genes", True), ("Callers", True), ("ARTS known", True),
        ("ARTS DUF", True), ("MIBiG comparison", False),
        ("Core-gene evidence", False),
        ("BGC class signals", False), ("Product signals", False),
        ("Edge distance (bp)", True),
    )
    header = "".join(
        "<th>{label}</th>".format(label=escape(label)) if index == 0 else
        "<th><button class='sort-header' type='button' data-column='{index}' "
        "data-type='{type}'>{label}</button></th>".format(
            index=index, type="number" if numeric else "text", label=escape(label)
        )
        for index, (label, numeric) in enumerate(columns)
    )
    rendered_rows = []
    for _, row in frame.iterrows():
        accession = clean(row.get("best_mibig_id"))
        tools = number(row.get("support_count"))
        known = number(row.get("arts_known_hits"))
        duf = number(row.get("arts_duf_hits"))
        callers = clean(row.get("support_tools")).lower()
        label = clean(row.get("consensus_label")) or clean(row.get("consensus_id"))
        fields = [
            (gene_map_button(row, map_index), ""),
            ("<span class='region-name'>{}</span>".format(escape(label)), label),
            (escape(clean(row.get("contig"))), clean(row.get("contig"))),
            (escape(clean(row.get("start"))), number(row.get("start"))),
            (escape(clean(row.get("end"))), number(row.get("end"))),
            (escape(clean(row.get("length_bp"))), number(row.get("length_bp"))),
            (str(number(map_index.get(clean(row.get("consensus_id")), {}).get("gene_count"))),
             number(map_index.get(clean(row.get("consensus_id")), {}).get("gene_count"))),
            ("<div class='caller-markers'>{}</div>".format(caller_markers(row)), tools),
            (str(known), known), (str(duf), duf),
            (mibig_cell(row), accession),
            ("<span class='signal-cell'>{}</span>".format(escape(clean(row.get("core_gene_evidence")))), clean(row.get("core_gene_evidence"))),
            ("<span class='signal-cell'>{}</span>".format(escape(clean(row.get("bgc_types")))), clean(row.get("bgc_types"))),
            ("<span class='signal-cell'>{}</span>".format(escape(clean(row.get("products")))), clean(row.get("products"))),
            (escape(clean(row.get("nearest_contig_edge_bp"))), number(row.get("nearest_contig_edge_bp"))),
        ]
        cells = "".join(
            "<td data-sort='{sort}'>{body}</td>".format(sort=escape(str(value), quote=True), body=html)
            for html, value in fields
        )
        searchable = " ".join(clean(row.get(key)) for key in (
            "consensus_label", "consensus_id", "contig", "bgc_types", "products",
            "best_mibig_id", "best_mibig_product", "core_gene_evidence",
        )).lower()
        rendered_rows.append(
            "<tr data-search='{search}' data-tools='{tools}' data-callers='{callers}' "
            "data-mibig='{mibig}' data-known='{known}' data-duf='{duf}'>{cells}</tr>".format(
                search=escape(searchable, quote=True), tools=tools,
                callers=escape(callers, quote=True), mibig="1" if accession else "0",
                known=known, duf=duf, cells=cells,
            )
        )
    panel = (
        "<section class='panel evidence-panel'><div class='section-head'>"
        "<h2>Cluster evidence</h2><span class='section-accent'></span></div>"
        "<p class='muted evidence-description'>One row per grouped candidate locus. "
        "Cross-caller groups use reciprocal overlap of at least {overlap}, nesting, "
        "or a shared core gene. Core genes come first from explicit antiSMASH "
        "biosynthetic roles or scaffold-forming Pfam evidence from GECCO/DeepBGC; "
        "a strict Bakta annotation fallback is used "
        "only when a predictor has no resolvable core-gene role. "
        "Merged boundaries are approximate. ARTS overlaps and MIBiG comparisons "
        "are separate observations.</p>"
        "<div class='evidence-toolbar'>"
        "<label>Search<input id='evidence-search' type='search' placeholder='Region, contig, core gene, class, product, MIBiG'></label>"
        "<label>Tool count<select id='evidence-tool-count'><option value=''>All</option>"
        "<option value='1'>1 tool</option><option value='2'>2 tools</option><option value='3'>3 tools</option></select></label>"
        "<label>Caller<select id='evidence-caller'><option value=''>Any</option>"
        "<option value='antismash'>antiSMASH</option><option value='gecco'>GECCO</option>"
        "<option value='deepbgc'>DeepBGC</option></select></label>"
        "<label class='evidence-check'><input id='evidence-mibig' type='checkbox'>MIBiG comparison</label>"
        "<label class='evidence-check'><input id='evidence-known' type='checkbox'>ARTS known hit</label>"
        "<label class='evidence-check'><input id='evidence-duf' type='checkbox'>ARTS DUF hit</label>"
        "</div><p class='evidence-count' id='evidence-count'></p>"
        "<div class='table-wrap'><table class='table evidence-table' id='cluster-evidence-table'>"
        "<thead><tr>{header}</tr></thead><tbody>{rows}</tbody></table></div></section>"
    ).format(
        header=header,
        rows="".join(rendered_rows),
        overlap=escape("{:.2f}".format(overlap_fraction)),
    )
    return EVIDENCE_CSS + panel + EVIDENCE_JS
