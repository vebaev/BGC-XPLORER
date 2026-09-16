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
        ("Method", "evidence_source"), ("Product", "best_mibig_product"),
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
.column-legend { display: grid; grid-template-columns: repeat(auto-fit, minmax(255px, 1fr));
  gap: 5px 22px; margin: 0 0 20px; padding: 14px 16px; border-radius: 12px;
  background: var(--accent-soft); font-size: 12.5px; line-height: 1.45; }
.column-legend > div { display: flex; gap: 7px; align-items: baseline; }
.column-legend dt { font-weight: 700; color: var(--text-strong); white-space: nowrap; }
.column-legend dd { margin: 0; color: var(--muted); }
.evidence-toolbar { display: flex; flex-wrap: wrap; gap: 12px; align-items: end; margin-bottom: 12px; }
.evidence-toolbar label { display: flex; flex-direction: column; gap: 5px; color: var(--muted); font-size: 12px; font-weight: 600; }
.evidence-toolbar label.evidence-check { flex-direction: row; align-items: center; gap: 7px; padding: 8px 2px; font-size: 13px; }
.evidence-toolbar input[type='search'], .evidence-toolbar select { font: inherit; font-size: 14px; color: var(--text-strong); background: white; border: 1px solid var(--line-strong); border-radius: 9px; padding: 9px 10px; min-height: 39px; }
.evidence-toolbar input[type='search'] { width: min(290px, 80vw); }
.evidence-toolbar input[type='checkbox'] { accent-color: var(--accent); }
.evidence-count { margin: 0 0 12px; color: var(--muted); font-size: 13px; }
.evidence-panel .table-wrap { border: 1px solid var(--line); border-radius: 12px; }
.evidence-table { min-width: 980px; }
.evidence-table td { vertical-align: top; }
.loc-cell { display: flex; flex-direction: column; gap: 2px; min-width: 150px; }
.loc-range { font-variant-numeric: tabular-nums; white-space: nowrap; font-weight: 600; }
.loc-meta { font-size: 12px; color: var(--muted); white-space: nowrap; }
.arts-cell { display: flex; flex-direction: column; gap: 4px; align-items: flex-start; }
.arts-badge { font-size: 11px; font-weight: 700; line-height: 1; white-space: nowrap;
  padding: 4px 7px; border-radius: 6px; }
.arts-known { color: #8a4a06; background: rgba(255, 123, 29, 0.14); }
.arts-duf { color: #1d6b74; background: rgba(65, 188, 199, 0.16); }
.signal-badge { display: inline-block; font-size: 11px; font-weight: 700; line-height: 1;
  padding: 4px 7px; border-radius: 6px; margin: 3px 4px 0 0; white-space: nowrap; }
.signal-activity { color: #8f2743; background: rgba(235, 91, 120, 0.14); }
.signal-other { color: var(--muted); background: rgba(94, 108, 152, 0.12); }
.signal-classes { display: block; }
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
.evidence-table .signal-cell { max-width: 330px; min-width: 210px; font-size: 12.5px; line-height: 1.5; }
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


ACTIVITY_TERMS = {
    "antibacterial", "antifungal", "cytotoxic", "inhibitor",
    "antibacterial-cytotoxic", "antibacterial-inhibitor",
}


def split_terms(value):
    return [term.strip() for term in clean(value).split(",") if term.strip()]


def format_span(value):
    return "{:,}".format(value).replace(",", "\u202f")


def format_length(value):
    if value >= 1000:
        return "{:.1f} kb".format(value / 1000.0)
    return "{} bp".format(format_span(value))


def location_html(contig, start, end, length=0):
    """Contig, coordinates and length in one cell instead of four columns."""
    start = number(start)
    end = number(end)
    length = number(length) or max(end - start + 1, 0)
    contig = clean(contig)
    return (
        "<div class='loc-cell'>"
        "<span class='loc-range'>{start}\u2013{end}</span>"
        "<span class='loc-meta'>{contig}{sep}{length}</span>"
        "</div>"
    ).format(
        start=escape(format_span(start)), end=escape(format_span(end)),
        contig=escape(contig), sep=" \u00b7 " if contig else "",
        length=escape(format_length(length)),
    )


def location_cell(row):
    return location_html(
        row.get("contig"), row.get("start"), row.get("end"), row.get("length_bp")
    )


def arts_cell(known, duf):
    """ARTS known and DUF counts as two badges in a single column."""
    badges = []
    if known:
        badges.append("<span class='arts-badge arts-known'>{} known</span>".format(known))
    if duf:
        badges.append("<span class='arts-badge arts-duf'>{} DUF</span>".format(duf))
    if not badges:
        return "<span class='muted'>\u2014</span>"
    return "<div class='arts-cell'>{}</div>".format("".join(badges))


def signal_cell(row):
    """BGC classes plus the activity terms that only DeepBGC contributes.

    The old Product signals column repeated the class list verbatim in 32 of 53
    regions; the only thing it added was the predicted activity, which now reads
    as a badge instead of a second wide column.
    """
    seen = set()
    classes = []
    for term in split_terms(row.get("bgc_types")):
        if term.lower() not in seen:
            seen.add(term.lower())
            classes.append(term)
    extras = []
    for term in split_terms(row.get("products")):
        if term.lower() not in seen:
            seen.add(term.lower())
            extras.append(term)
    parts = []
    if classes:
        parts.append("<span class='signal-classes'>{}</span>".format(escape(", ".join(classes))))
    for term in extras:
        tone = "activity" if term.lower() in ACTIVITY_TERMS else "other"
        parts.append("<span class='signal-badge signal-{tone}'>{term}</span>".format(
            tone=tone, term=escape(term)
        ))
    if not parts:
        return "<span class='muted'>\u2014</span>", ""
    return "<div class='signal-cell'>{}</div>".format("".join(parts)), ", ".join(classes + extras)


def render_evidence_table(frame, map_index, gene_map_button, min_containment=0.80, standalone=True):
    columns = (
        ("Gene map", False), ("Region", False), ("Location", True),
        ("Genes", True), ("Callers", True), ("ARTS", True),
        ("MIBiG", False), ("Signals", False),
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
        genes = number(map_index.get(clean(row.get("consensus_id")), {}).get("gene_count"))
        signals_html, signals_sort = signal_cell(row)
        fields = [
            (gene_map_button(row, map_index), ""),
            ("<span class='region-name'>{}</span>".format(escape(label)), label),
            (location_cell(row), number(row.get("start"))),
            (str(genes), genes),
            ("<div class='caller-markers'>{}</div>".format(caller_markers(row)), tools),
            # Sorting ranks known hits first and uses DUF only to break ties.
            (arts_cell(known, duf), known * 1000 + duf),
            (mibig_cell(row), accession),
            (signals_html, signals_sort),
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
    # Inside a tab the surrounding card supplies the chrome, so the panel drops
    # its own <section class='panel'> wrapper.
    open_tag = "<section class='panel evidence-panel'>" if standalone else "<div class='evidence-panel'>"
    close_tag = "</section>" if standalone else "</div>"
    panel = (
        "{open_tag}"
        "<dl class='column-legend'>"
        "<div><dt>Gene map</dt><dd>opens the annotated gene diagram</dd></div>"
        "<div><dt>Region</dt><dd>grouped locus; members cover \u2265{overlap} of the shorter region</dd></div>"
        "<div><dt>Location</dt><dd>coordinates, contig and length</dd></div>"
        "<div><dt>Genes</dt><dd>annotated genes inside the region</dd></div>"
        "<div><dt>Callers</dt><dd>tools that predicted it</dd></div>"
        "<div><dt>ARTS</dt><dd>known resistance-model hits and DUF hits</dd></div>"
        "<div><dt>MIBiG</dt><dd>closest characterised cluster; expand for method and product</dd></div>"
        "<div><dt>Signals</dt><dd>predicted BGC classes, activity as a badge</dd></div>"
        "</dl>"
        "<div class='evidence-toolbar'>"
        "<label>Search<input id='evidence-search' type='search' placeholder='Region, contig, class, product, MIBiG'></label>"
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
        "<thead><tr>{header}</tr></thead><tbody>{rows}</tbody></table></div>{close_tag}"
    ).format(
        open_tag=open_tag,
        close_tag=close_tag,
        header=header,
        rows="".join(rendered_rows),
        overlap=escape("{:.0f}%".format(min_containment * 100)),
    )
    return EVIDENCE_CSS + panel + EVIDENCE_JS
