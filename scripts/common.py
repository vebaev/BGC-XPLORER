import json
import os
from html import escape
from pathlib import Path

import pandas as pd


AI_FETCH_RETRY_JS = Path(__file__).with_name("ai_fetch_retry.js").read_text(encoding="utf-8")


BGC_COLUMNS = [
    "sample",
    "tool",
    "contig",
    "start",
    "end",
    "strand",
    "bgc_id",
    "bgc_type",
    "product",
    "score",
    "confidence",
    "core_gene_records",
    "source_file",
]


ARTS_COLUMNS = [
    "sample",
    "tool",
    "contig",
    "start",
    "end",
    "feature",
    "score",
    "evidence",
    "source_file",
]


def ensure_parent(path):
    os.makedirs(os.path.dirname(path), exist_ok=True)


def empty_bgc_frame():
    return pd.DataFrame(columns=BGC_COLUMNS)


def empty_arts_frame():
    return pd.DataFrame(columns=ARTS_COLUMNS)


def load_table_if_exists(path, columns):
    if os.path.exists(path):
        try:
            df = pd.read_csv(path, sep="\t")
        except pd.errors.EmptyDataError:
            return pd.DataFrame(columns=columns)
        missing = [col for col in columns if col not in df.columns]
        for col in missing:
            df[col] = ""
        return df[columns]
    return pd.DataFrame(columns=columns)


def first_existing_path(base_dir, patterns):
    for pattern in patterns:
        matches = sorted(base_dir.glob(pattern))
        if matches:
            return matches[0]
    return None


def first_existing_column(df, names, default=""):
    for name in names:
        if name in df.columns:
            return df[name]
    if len(df.index) == 0:
        return pd.Series(dtype=object)
    return pd.Series([default] * len(df.index), index=df.index)


def coerce_frame_columns(df, columns, default=""):
    for column in columns:
        if column not in df.columns:
            df[column] = default
    return df[columns]


def format_consensus_label(value):
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return ""
    return text.replace("_consensus_", "_")


def add_consensus_label(df, source_col="consensus_id", target_col="consensus_label"):
    frame = df.copy()
    if target_col not in frame.columns:
        frame[target_col] = ""
    if source_col in frame.columns:
        frame[target_col] = frame[source_col].map(format_consensus_label)
    return frame


def write_tsv(df, path):
    ensure_parent(path)
    df.to_csv(path, sep="\t", index=False)


def write_json(data, path):
    ensure_parent(path)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2)


def read_json(path):
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def df_to_html_table(df):
    if df.empty:
        return "<p>No records found.</p>"
    return '<div class="table-wrap">{0}</div>'.format(
        df.to_html(index=False, classes="table", border=0, escape=True)
    )


def html_page(title, sections):
    blocks = "\n".join(sections)
    return """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{title}</title>
  <style>
    :root {{
      --bg: #f5f7ff;
      --bg-soft: #fbfcff;
      --panel: rgba(255, 255, 255, 0.94);
      --panel-2: #ffffff;
      --panel-3: #eef2ff;
      --line: rgba(94, 108, 152, 0.12);
      --line-strong: rgba(94, 108, 152, 0.22);
      --text: #273149;
      --text-strong: #161f33;
      --muted: #5e6c8f;
      --muted-2: #7d88a7;
      --accent: #5f57ff;
      --accent-soft: rgba(95, 87, 255, 0.1);
      --accent-deep: #4a42eb;
      --violet: #6a5cff;
      --emerald: #73c95b;
      --amber: #f59b38;
      --rose: #eb5b78;
      --teal: #41bcc7;
      --blue: #4f8dfd;
      --orange: #ff7b1d;
      --sky: #4a7fff;
      --lime: #7cc85d;
      --cyan: #35bfd3;
      --indigo: #7c6bff;
      --radius: 20px;
      --radius-sm: 14px;
      --transition: 0.18s ease;
      --sans: "Sora", "Avenir Next", "Segoe UI", sans-serif;
      --mono: "SF Mono", "JetBrains Mono", "Fira Code", ui-monospace, Menlo, monospace;
      --shadow-card: 0 12px 35px rgba(89, 104, 146, 0.08);
      --shadow-soft: 0 8px 24px rgba(89, 104, 146, 0.06);
      --shadow-modal: 0 30px 80px rgba(55, 69, 112, 0.22);
    }}
    * {{
      box-sizing: border-box;
    }}
    *::-webkit-scrollbar {{
      width: 10px;
      height: 10px;
    }}
    *::-webkit-scrollbar-thumb {{
      background: rgba(255, 255, 255, 0.12);
      border-radius: 999px;
      border: 2px solid transparent;
      background-clip: padding-box;
    }}
    *::-webkit-scrollbar-thumb:hover {{
      background: rgba(255, 255, 255, 0.22);
      background-clip: padding-box;
    }}
    * {{
      scrollbar-width: thin;
      scrollbar-color: rgba(255, 255, 255, 0.15) transparent;
    }}
    ::selection {{
      background: rgba(34, 211, 238, 0.25);
      color: var(--text-strong);
    }}
    body {{
      font-family: var(--sans);
      margin: 0;
      color: var(--text);
      background:
        radial-gradient(900px 440px at 8% 0%, rgba(95, 87, 255, 0.09), transparent 62%),
        linear-gradient(180deg, #fcfdff 0%, #f3f6ff 100%),
        var(--bg);
      -webkit-font-smoothing: antialiased;
      -moz-osx-font-smoothing: grayscale;
      text-rendering: optimizeLegibility;
    }}
    .page {{
      max-width: 1220px;
      margin: 0 auto;
      padding: 32px 28px 48px;
    }}
    h1, h2, h3 {{
      color: var(--text-strong);
      margin: 0;
      line-height: 1.15;
      letter-spacing: 0;
    }}
    p {{
      color: var(--muted);
      line-height: 1.6;
      margin: 0;
    }}
    .hero-card {{
      display: flex;
      justify-content: space-between;
      gap: 22px;
      align-items: stretch;
      background: transparent;
      margin-bottom: 18px;
    }}
    .report-nav {{
      display: flex;
      justify-content: flex-start;
      margin-bottom: 18px;
    }}
    .report-home-link {{
      display: inline-flex;
      align-items: center;
      gap: 9px;
      padding: 10px 14px;
      border: 1px solid var(--line);
      border-radius: 12px;
      background: rgba(255, 255, 255, 0.82);
      color: var(--accent-deep);
      font-weight: 600;
      text-decoration: none;
      box-shadow: var(--shadow-soft);
    }}
    .report-home-link:hover {{
      border-color: rgba(95, 87, 255, 0.32);
      transform: translateY(-1px);
    }}
    .hero-brand {{
      display: flex;
      flex-direction: column;
      justify-content: center;
      align-items: center;
      min-width: 0;
      flex-shrink: 0;
    }}
    .brand-mark {{
      width: 56px;
      height: 56px;
      border-radius: 14px;
      border: 1px solid var(--line);
      display: grid;
      place-items: center;
      color: var(--accent);
      font-size: 24px;
      font-weight: 700;
      background: rgba(255, 255, 255, 0.75);
      box-shadow: var(--shadow-card);
    }}
    .brand-logo {{
      height: 112px;
      width: auto;
      border-radius: 14px;
      flex-shrink: 0;
      object-fit: contain;
    }}
    .eyebrow {{
      color: var(--accent);
      font-weight: 600;
    }}
    .hero-panel {{
      flex: 1;
      min-height: 112px;
      box-sizing: border-box;
      background: rgba(255, 255, 255, 0.8);
      border: 1px solid var(--line);
      border-radius: 16px;
      padding: 20px 24px;
      box-shadow: var(--shadow-card);
      display: flex;
      flex-direction: column;
      justify-content: center;
      min-width: 0;
    }}
    .hero-panel h1 {{
      font-size: 28px;
      font-weight: 400;
      line-height: 1.2;
      margin-bottom: 12px;
      letter-spacing: 0.01em;
      color: var(--text);
    }}
    .hero-panel p {{
      max-width: 620px;
      font-size: 16px;
      color: var(--text);
    }}
    .meta-badge {{
      min-width: 190px;
      height: 112px;
      box-sizing: border-box;
      background: rgba(255, 255, 255, 0.8);
      border: 1px solid var(--line);
      border-radius: 16px;
      padding: 18px 20px;
      box-shadow: var(--shadow-card);
      display: flex;
      flex-direction: column;
      justify-content: center;
    }}
    .meta-badge-label {{
      color: var(--muted);
      font-size: 13px;
      margin-bottom: 8px;
    }}
    .meta-badge-value {{
      color: var(--text-strong);
      font-size: 18px;
      font-weight: 600;
      line-height: 1.4;
      margin-bottom: 8px;
    }}
    .meta-badge-foot {{
      color: var(--muted);
      font-size: 12px;
      font-weight: 500;
      letter-spacing: 0.04em;
      text-transform: uppercase;
    }}
    .panel,
    .glance-shell,
    .info-card,
    .table-panel {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 18px;
      box-shadow: var(--shadow-card);
    }}
    .panel,
    .glance-shell,
    .table-panel {{
      padding: 22px 24px;
      margin-bottom: 18px;
    }}
    .section-head {{
      display: flex;
      flex-direction: row;
      align-items: center;
      justify-content: space-between;
      gap: 16px;
      margin-bottom: 18px;
    }}
    .section-head h2 {{
      font-size: 30px;
      font-weight: 600;
      letter-spacing: -0.02em;
      color: var(--text-strong);
    }}
    .section-accent {{
      width: 28px;
      height: 4px;
      border-radius: 999px;
      background: linear-gradient(90deg, var(--accent), rgba(95, 87, 255, 0.15));
      flex: 0 0 auto;
    }}
    .metrics-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
      gap: 14px;
    }}
    .metric-card {{
      min-height: 136px;
      padding: 18px 18px 16px;
      border-radius: 16px;
      border: 1px solid var(--line);
      display: flex;
      gap: 14px;
      align-items: flex-start;
      background: linear-gradient(180deg, rgba(255,255,255,0.98), rgba(244,247,255,0.96));
    }}
    .metric-icon {{
      width: 48px;
      height: 48px;
      border-radius: 14px;
      display: grid;
      place-items: center;
      font-size: 24px;
      font-weight: 600;
      flex: 0 0 auto;
    }}
    .metric-label {{
      color: var(--text-strong);
      font-size: 14px;
      font-weight: 700;
      text-transform: uppercase;
      margin-bottom: 10px;
    }}
    .metric-value {{
      font-size: 28px;
      font-weight: 700;
      color: var(--text-strong);
      line-height: 1;
      font-variant-numeric: tabular-nums;
      margin-bottom: 8px;
    }}
    .metric-copy p {{
      font-size: 13px;
      color: var(--muted);
      line-height: 1.45;
    }}
    .metric-violet .metric-icon {{ color: var(--violet); background: rgba(106, 92, 255, 0.1); }}
    .metric-green .metric-icon {{ color: var(--emerald); background: rgba(115, 201, 91, 0.12); }}
    .metric-blue .metric-icon {{ color: var(--blue); background: rgba(79, 141, 253, 0.12); }}
    .metric-orange .metric-icon {{ color: var(--orange); background: rgba(255, 123, 29, 0.12); }}
    .metric-teal .metric-icon {{ color: var(--teal); background: rgba(65, 188, 199, 0.12); }}
    .metric-indigo .metric-icon {{ color: var(--indigo); background: rgba(124, 107, 255, 0.12); }}
    .metric-lime .metric-icon {{ color: var(--lime); background: rgba(124, 200, 93, 0.12); }}
    .metric-sky .metric-icon {{ color: var(--sky); background: rgba(74, 127, 255, 0.12); }}
    .metric-amber .metric-icon {{ color: var(--amber); background: rgba(245, 155, 56, 0.12); }}
    .metric-cyan .metric-icon {{ color: var(--cyan); background: rgba(53, 191, 211, 0.12); }}
    .glance-shell .metric-icon {{
      color: var(--accent) !important;
      background: var(--accent-soft) !important;
      border: 1px solid rgba(95, 87, 255, 0.16);
    }}
    .model-status {{
      margin-top: 10px;
    }}
    .info-grid,
    .chart-panels {{
      display: grid;
      grid-template-columns: repeat(2, minmax(0, 1fr));
      gap: 18px;
      margin-bottom: 18px;
    }}
    /* Cards inside a grid are spaced by the grid gap; their own bottom margin
       would double it and make those rows sit further apart than the rest. */
    .info-grid > *,
    .chart-panels > * {{
      margin-bottom: 0;
    }}
    .info-grid.summary-only {{
      grid-template-columns: 1fr;
    }}
    .info-card {{
      display: flex;
      gap: 18px;
      align-items: flex-start;
      padding: 24px;
      min-height: 150px;
    }}
    .info-icon {{
      width: 56px;
      height: 56px;
      border-radius: 16px;
      display: grid;
      place-items: center;
      color: var(--accent);
      background: rgba(95, 87, 255, 0.12);
      flex: 0 0 auto;
      font-size: 25px;
      font-weight: 700;
    }}
    .summary-only .info-card {{
      padding-left: 42px;
    }}
    .summary-only .info-icon {{
      width: 48px;
      height: 48px;
      border-radius: 14px;
      font-size: 24px;
      font-weight: 600;
      background: var(--accent-soft);
      border: 1px solid rgba(95, 87, 255, 0.16);
    }}
    .info-copy h3 {{
      font-size: 26px;
      letter-spacing: -0.03em;
      margin-bottom: 10px;
    }}
    .info-copy p {{
      color: var(--text);
      font-size: 15px;
    }}
    .chart-panel {{
      margin: 0;
    }}
    .chart-grid {{
      display: grid;
      grid-template-columns: minmax(220px, 0.95fr) minmax(220px, 1fr);
      gap: 24px;
      align-items: center;
    }}
    .donut-shell {{
      display: flex;
      justify-content: center;
      align-items: center;
    }}
    .donut-chart {{
      width: 220px;
      height: 220px;
      border-radius: 50%;
      display: grid;
      place-items: center;
    }}
    .donut-hole {{
      width: 128px;
      height: 128px;
      border-radius: 50%;
      background: #ffffff;
      display: flex;
      flex-direction: column;
      align-items: center;
      justify-content: center;
      box-shadow: inset 0 0 0 1px rgba(94, 108, 152, 0.08);
    }}
    .donut-hole strong {{
      font-size: 42px;
      color: var(--text-strong);
      line-height: 1;
      font-variant-numeric: tabular-nums;
    }}
    .donut-hole span {{
      color: var(--muted);
      font-size: 14px;
    }}
    .legend-list {{
      list-style: none;
      margin: 0;
      padding: 0;
      display: grid;
      gap: 14px;
    }}
    .legend-list li {{
      display: grid;
      grid-template-columns: 14px 1fr auto;
      gap: 12px;
      align-items: center;
      color: var(--text);
      font-size: 15px;
    }}
    .legend-swatch {{
      width: 10px;
      height: 10px;
      border-radius: 50%;
      display: inline-block;
    }}
    .legend-label {{
      color: var(--text-strong);
    }}
    .legend-value {{
      color: var(--text);
      font-variant-numeric: tabular-nums;
    }}
    .table-panel {{
      padding: 0;
      overflow: hidden;
    }}
    .tabs-nav {{
      display: flex;
      gap: 0;
      border-bottom: 1px solid var(--line);
      padding: 0 18px;
      overflow-x: auto;
    }}
    .tab-btn {{
      appearance: none;
      border: 0;
      background: transparent;
      color: var(--muted);
      font-family: var(--sans);
      font-size: 15px;
      font-weight: 600;
      padding: 16px 18px;
      cursor: pointer;
      display: inline-flex;
      align-items: center;
      gap: 10px;
      white-space: nowrap;
      border-bottom: 3px solid transparent;
    }}
    .tab-btn.is-active {{
      color: var(--accent);
      border-bottom-color: var(--accent);
    }}
    .tab-pane {{
      padding: 0 0 8px;
    }}
    .tab-pane[hidden] {{
      display: none;
    }}
    .tab-icon {{
      font-size: 16px;
      line-height: 1;
    }}
    .table-note {{
      padding: 18px 24px 8px;
    }}
    .table-wrap {{
      overflow-x: auto;
      border-radius: 0 0 18px 18px;
      background: transparent;
    }}
    table {{
      border-collapse: collapse;
      width: 100%;
      font-size: 14px;
      min-width: 900px;
    }}
    th, td {{
      border-bottom: 1px solid var(--line);
      padding: 15px 18px;
      text-align: left;
      vertical-align: top;
      color: var(--text);
      background: rgba(255, 255, 255, 0.72);
    }}
    th {{
      position: sticky;
      top: 0;
      background: #f9faff;
      color: var(--muted);
      font-size: 11px;
      font-weight: 600;
      text-transform: uppercase;
      letter-spacing: 0.08em;
      border-bottom: 1px solid var(--line-strong);
    }}
    td {{
      font-variant-numeric: tabular-nums;
    }}
    tr:nth-child(even) td {{
      background: rgba(246, 248, 255, 0.88);
    }}
    tr:hover td {{
      background: rgba(95, 87, 255, 0.05);
    }}
    pre {{
      white-space: pre-wrap;
      word-break: break-word;
      margin: 0;
      font-size: 13px;
      font-family: var(--mono);
      color: var(--text);
      background: var(--panel-2);
      border: 1px solid var(--line);
      border-radius: var(--radius-sm);
      padding: 14px;
      line-height: 1.5;
    }}
    code {{
      font-family: var(--mono);
      font-size: 0.88em;
      background: rgba(95, 87, 255, 0.08);
      color: var(--accent);
      padding: 1px 6px;
      border-radius: 5px;
      border: 1px solid var(--line);
    }}
    .muted {{
      color: var(--muted);
    }}
    .report-footer {{
      display: flex;
      justify-content: center;
      gap: 18px;
      align-items: center;
      color: var(--muted);
      font-size: 15px;
      padding: 16px 0 0;
    }}
    .report-footer .dot {{
      width: 4px;
      height: 4px;
      border-radius: 50%;
      background: var(--muted-2);
      display: inline-block;
    }}
    .gene-map-btn,
    .gene-map-close,
    .ai-analysis-btn {{
      appearance: none;
      border: 1px solid var(--line);
      background: #ffffff;
      color: var(--accent);
      border-radius: 12px;
      padding: 10px 14px;
      font-family: var(--sans);
      font-size: 13px;
      font-weight: 600;
      cursor: pointer;
      white-space: nowrap;
      transition: background var(--transition), border-color var(--transition), color var(--transition), transform var(--transition);
    }}
    .gene-map-btn:hover {{
      background: rgba(95, 87, 255, 0.08);
      border-color: rgba(95, 87, 255, 0.35);
      color: var(--accent);
      transform: translateY(-1px);
    }}
    .gene-map-close:hover {{
      background: rgba(235, 91, 120, 0.08);
      border-color: rgba(235, 91, 120, 0.35);
      color: var(--rose);
    }}
    .ai-analysis-btn {{
      border-color: rgba(95, 87, 255, 0.24);
      background: rgba(95, 87, 255, 0.1);
      color: var(--accent);
    }}
    .ai-analysis-btn:hover {{
      background: rgba(95, 87, 255, 0.16);
      border-color: rgba(95, 87, 255, 0.42);
    }}
    .ai-analysis-btn:disabled {{
      opacity: 0.5;
      cursor: not-allowed;
    }}
    .ai-analysis-note {{
      margin: 8px 0 0;
      color: var(--muted);
      font-size: 11px;
      line-height: 1.45;
    }}
    .gene-map-viewer {{
      position: fixed;
      z-index: 1000;
      top: clamp(10px, 3vh, 28px);
      left: clamp(10px, 3vw, 34px);
      right: clamp(10px, 3vw, 34px);
      max-height: calc(100vh - 56px);
      overflow: auto;
      margin: 0;
      background: rgba(255, 255, 255, 0.98);
      border: 1px solid var(--line-strong);
      border-radius: 22px;
      box-shadow: var(--shadow-modal);
      padding: 22px 24px;
    }}
    body.gene-map-open {{
      overflow: hidden;
    }}
    body.gene-map-open::before {{
      content: "";
      position: fixed;
      inset: 0;
      background: rgba(82, 94, 133, 0.18);
      backdrop-filter: blur(3px);
      -webkit-backdrop-filter: blur(3px);
      z-index: 999;
    }}
    .gene-map-head {{
      display: flex;
      justify-content: space-between;
      gap: 16px;
      align-items: flex-start;
      margin-bottom: 16px;
    }}
    .gene-map-head h2 {{
      font-size: 20px;
    }}
    .gene-map-stage {{
      overflow: auto;
      border: 1px solid var(--line);
      border-radius: 12px;
      background: #f8faff;
      margin-bottom: 14px;
      max-height: calc(100vh - 210px);
    }}
    .gene-map-stage img {{
      display: block;
      width: 100%;
      min-width: 980px;
      height: auto;
    }}
    .gene-map-stage img[hidden] {{
      display: none;
    }}
    .gene-map-inline {{
      min-width: 980px;
    }}
    .gene-map-inline svg {{
      display: block;
      width: 100%;
      height: auto;
    }}
    .gene-map-interpretation {{
      border-top: 1px solid var(--line);
      padding: 16px 18px 18px;
      background: #ffffff;
      min-width: 980px;
    }}
    .gene-map-interpretation h3 {{
      font-size: 13px;
      margin: 0 0 8px;
      color: var(--text-strong);
      text-transform: uppercase;
      letter-spacing: 0.08em;
      font-weight: 600;
    }}
    .gene-map-interpretation p {{
      max-width: 110ch;
      color: var(--text);
      font-size: 14px;
      line-height: 1.6;
    }}
    .gene-map-gene-table {{
      border-top: 1px solid var(--line);
      padding: 16px 18px 18px;
      background: #ffffff;
      min-width: 980px;
    }}
    .gene-map-gene-table h3 {{
      font-size: 13px;
      margin: 0 0 8px;
      color: var(--text-strong);
      text-transform: uppercase;
      letter-spacing: 0.08em;
      font-weight: 600;
    }}
    .gene-map-gene-table[hidden] {{
      display: none;
    }}
    .gene-table-wrap {{
      max-height: 320px;
      overflow: auto;
      border: 1px solid var(--line);
      border-radius: var(--radius-sm);
    }}
    .gene-table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 13px;
    }}
    .gene-table thead th {{
      position: sticky;
      top: 0;
      background: #f5f7fb;
      color: var(--text-strong);
      text-align: left;
      padding: 8px 10px;
      border-bottom: 1px solid var(--line);
      font-weight: 600;
      font-size: 12px;
      text-transform: uppercase;
      letter-spacing: 0.04em;
      white-space: nowrap;
      z-index: 1;
    }}
    .gene-table tbody td {{
      padding: 7px 10px;
      border-bottom: 1px solid rgba(94, 108, 152, 0.10);
      color: var(--text);
      vertical-align: top;
    }}
    .gene-table tbody tr:hover td {{
      background: rgba(99, 91, 255, 0.05);
    }}
    .gene-table tbody td.gene-locus {{
      font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      font-size: 12px;
      white-space: nowrap;
    }}
    .gene-table tbody td.gene-name {{
      font-weight: 600;
      color: var(--text-strong);
      white-space: nowrap;
    }}
    .gene-table tbody td.gene-product {{
      max-width: 240px;
    }}
    .gene-table tbody td.gene-eggnog {{
      max-width: 280px;
    }}
    .gene-table tbody td.gene-dbcan {{
      white-space: nowrap;
    }}
    .gene-map-ai {{
      border-top: 1px solid var(--line);
      padding: 16px 18px 18px;
      background: #ffffff;
      min-width: 980px;
    }}
    .gene-map-ai-head {{
      display: flex;
      gap: 10px;
      align-items: center;
      flex-wrap: wrap;
    }}
    .ai-spinner {{
      display: inline-block;
      width: 16px;
      height: 16px;
      border: 2px solid rgba(95, 87, 255, 0.18);
      border-top-color: var(--accent);
      border-radius: 50%;
      animation: ai-spin 0.7s linear infinite;
      vertical-align: middle;
    }}
    .ai-spinner[hidden] {{
      display: none;
    }}
    @keyframes ai-spin {{
      to {{
        transform: rotate(360deg);
      }}
    }}
    .gene-map-error {{
      color: #b4234d;
      font-size: 13px;
      font-weight: 600;
      margin: 0 0 12px;
      padding: 10px 12px;
      background: rgba(235, 91, 120, 0.08);
      border: 1px solid rgba(235, 91, 120, 0.24);
      border-radius: 8px;
      line-height: 1.5;
    }}
    .ai-analysis-result {{
      display: grid;
      gap: 12px;
      margin-top: 14px;
    }}
    .ai-analysis-result:empty {{
      display: none;
    }}
    .ai-analysis-result h3 {{
      font-size: 13px;
      margin: 0 0 6px;
      color: var(--text-strong);
      text-transform: uppercase;
      letter-spacing: 0.08em;
      font-weight: 600;
    }}
    .ai-analysis-result p {{
      color: var(--text);
      font-size: 14px;
      line-height: 1.6;
    }}
    .ai-analysis-result ul {{
      margin: 0;
      padding-left: 18px;
      color: var(--text);
      font-size: 14px;
    }}
    .ai-analysis-result li {{
      margin: 5px 0;
      line-height: 1.5;
    }}
    .ai-analysis-result li::marker {{
      color: var(--accent);
    }}
    .ai-status {{
      border: 1px solid var(--line);
      border-left: 3px solid var(--accent);
      border-radius: var(--radius-sm);
      background: #ffffff;
      padding: 14px 16px;
      color: var(--text);
    }}
    .ai-status-error {{
      border-left-color: var(--rose);
      background: rgba(235, 91, 120, 0.06);
      color: #b4234d;
    }}
    .callout {{
      background: #ffffff;
      border: 1px solid var(--line);
      border-left: 3px solid var(--accent);
      border-radius: var(--radius-sm);
      padding: 14px 16px;
    }}
    .callout p {{
      color: var(--text);
      font-size: 14px;
      line-height: 1.6;
      margin: 0;
    }}
    .list {{
      display: grid;
      gap: 10px;
    }}
    @media (max-width: 980px) {{
      .page {{
        padding: 20px 14px 40px;
      }}
      .hero-card,
      .info-grid,
      .chart-panels,
      .chart-grid {{
        grid-template-columns: 1fr;
        display: grid;
      }}
      .hero-card {{
        flex-wrap: wrap;
        gap: 14px;
      }}
      .hero-brand,
      .hero-panel,
      .meta-badge {{
        width: 100%;
        min-height: auto;
      }}
      .hero-panel h1 {{
        font-size: 20px;
      }}
      .glance-shell,
      .panel,
      .table-panel {{
        padding: 18px 16px;
      }}
      .summary-only .info-card {{
        padding-left: 34px;
      }}
      .tabs-nav {{
        padding: 0 10px;
      }}
      .table-note {{
        padding: 16px 16px 6px;
      }}
      .gene-map-viewer {{
        left: 8px;
        right: 8px;
        top: 8px;
        padding: 16px;
      }}
      .report-footer {{
        flex-wrap: wrap;
        gap: 10px;
      }}
    }}
  </style>
</head>
<body>
  <div class="page">
    {blocks}
  </div>
  <script>
    {ai_fetch_retry_js}
    (function () {{
      var viewer = document.getElementById('gene-map-viewer');
      var image = document.getElementById('gene-map-image');
      var inline = document.getElementById('gene-map-inline');
      var interpretation = document.getElementById('gene-map-interpretation');
      var title = document.getElementById('gene-map-title');
      var meta = document.getElementById('gene-map-meta');
      var summary = document.getElementById('gene-map-summary');
      var error = document.getElementById('gene-map-error');
      var aiButton = document.getElementById('gene-map-ai-button');
      var aiSpinner = document.getElementById('gene-map-ai-spinner');
      var aiStatus = document.getElementById('gene-map-ai-status');
      var aiResult = document.getElementById('gene-map-ai-result');
      var activeAiModel = document.getElementById('active-ai-model');
      var aiModelLabel = document.getElementById('ai-model-label');
      var svgDataNode = document.getElementById('gene-map-svg-data');
      var svgData = {{}};
      var activeButton = null;
      if (svgDataNode) {{
        try {{
          svgData = JSON.parse(svgDataNode.textContent || '{{}}');
        }} catch (parseError) {{
          svgData = {{}};
        }}
      }}
      var geneTableDataNode = document.getElementById('gene-map-table-data');
      var geneTableData = {{}};
      if (geneTableDataNode) {{
        try {{
          geneTableData = JSON.parse(geneTableDataNode.textContent || '{{}}');
        }} catch (parseError) {{
          geneTableData = {{}};
        }}
      }}
      var geneTableContainer = document.getElementById('gene-map-gene-table');
      var geneTableBody = geneTableContainer ? geneTableContainer.querySelector('tbody') : null;
      function renderGeneTable(clusterKey) {{
        if (!geneTableContainer || !geneTableBody) {{
          return;
        }}
        geneTableBody.innerHTML = '';
        var rows = (geneTableData && geneTableData[clusterKey]) || [];
        if (!rows.length) {{
          geneTableContainer.hidden = true;
          return;
        }}
        rows.forEach(function (row) {{
          var tr = document.createElement('tr');
          tr.innerHTML = '<td class="gene-locus">' + escapeHtml(row.locus) + '</td>'
            + '<td class="gene-name">' + escapeHtml(row.gene) + '</td>'
            + '<td class="gene-product">' + escapeHtml(row.product) + '</td>'
            + '<td class="gene-eggnog">' + escapeHtml(row.eggnog) + '</td>'
            + '<td class="gene-dbcan">' + escapeHtml(row.dbcan) + '</td>';
          geneTableBody.appendChild(tr);
        }});
        geneTableContainer.hidden = false;
      }}
      if (!viewer || !image || !inline || !interpretation || !title || !meta || !summary || !error || !aiButton || !aiSpinner || !aiStatus || !aiResult) {{
        var tabsOnly = true;
      }} else {{
        var tabsOnly = false;
      }}
      function activateTab(button) {{
        var target = button.getAttribute('data-tab-target');
        if (!target) {{
          return;
        }}
        document.querySelectorAll('.tab-btn[data-tab-target]').forEach(function (item) {{
          var active = item === button;
          item.classList.toggle('is-active', active);
          item.setAttribute('aria-selected', active ? 'true' : 'false');
        }});
        document.querySelectorAll('.tab-pane').forEach(function (pane) {{
          var active = pane.id === target;
          pane.classList.toggle('is-active', active);
          pane.hidden = !active;
        }});
      }}
      document.querySelectorAll('.tab-btn[data-tab-target]').forEach(function (button) {{
        button.addEventListener('click', function () {{
          activateTab(button);
        }});
      }});
      if (tabsOnly) {{
        return;
      }}
      function escapeHtml(value) {{
        return String(value || '').replace(/[&<>"']/g, function (char) {{
          return ({{'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;'}})[char];
        }});
      }}
      function listItems(items) {{
        if (!Array.isArray(items) || !items.length) {{
          return '<p class="muted">No entries returned.</p>';
        }}
        return '<ul>' + items.map(function (item) {{
          return '<li>' + escapeHtml(item) + '</li>';
        }}).join('') + '</ul>';
      }}
      function renderAnalysis(data) {{
        if (data.model) {{
          showActiveAiModel(data.model);
        }}
        var analysis = data.analysis || data;
        var sections = [
          ['Summary', analysis.summary],
          ['Likely product / function', analysis.likely_product_or_function],
          ['Biosynthetic logic', analysis.biosynthetic_logic],
          ['Resistance, transport and regulation', analysis.resistance_transport_regulation]
        ];
        var html = sections.map(function (section) {{
          if (!section[1]) {{
            return '';
          }}
          return '<div class="ai-status"><h3>' + escapeHtml(section[0]) + '</h3><p>' + escapeHtml(section[1]) + '</p></div>';
        }}).join('');
        html += '<div class="ai-status"><h3>Key genes</h3>' + listItems(analysis.key_genes) + '</div>';
        html += '<div class="ai-status"><h3>Recommended follow-up</h3>' + listItems(analysis.recommended_followup) + '</div>';
        aiResult.innerHTML = html;
      }}
      function statusEndpoint(endpoint, jobId) {{
        return endpoint.replace(/\/analyze_cluster\/?$/, '/analyze_cluster_status/' + encodeURIComponent(jobId));
      }}
      function diagnosticMessage(data) {{
        if (!data || typeof data !== 'object') {{
          return 'AI analysis request failed.';
        }}
        var parts = [];
        if (data.error) {{
          parts.push(String(data.error));
        }}
        ['error_type', 'stage', 'elapsed_seconds', 'timeout_seconds', 'payload_chars', 'gene_count', 'representative_gene_count', 'guided_json'].forEach(function (key) {{
          if (data[key] !== undefined && data[key] !== null && data[key] !== '') {{
            parts.push(key + '=' + data[key]);
          }}
        }});
        return parts.join(' | ') || data.message || 'AI analysis request failed.';
      }}
      function waitForAiJob(endpoint, jobId, attempt) {{
        aiStatus.textContent = 'Analyzing...';
        return fetchWithNetworkRetry(statusEndpoint(endpoint, jobId), {{
          method: 'GET',
          headers: {{'Accept': 'application/json'}}
        }}, 3, 750)
          .then(function (response) {{
            return response.json().then(function (data) {{
              if (response.status === 202 || data.status === 'queued' || data.status === 'running') {{
                if (attempt >= 120) {{
                  throw new Error('AI analysis is still running after 4 minutes. Try again in a moment or check the AI service diagnostics.');
                }}
                return new Promise(function (resolve) {{
                  window.setTimeout(resolve, 2000);
                }}).then(function () {{
                  return waitForAiJob(endpoint, jobId, attempt + 1);
                }});
              }}
              if (!response.ok) {{
                throw new Error(diagnosticMessage(data));
              }}
              return data;
            }});
          }});
      }}
      function aiRequestErrorMessage(endpoint, error) {{
        var detail = error && error.message ? error.message : 'request failed';
        if (/NVIDIA API error 503\\b/i.test(detail) && /Service temporarily overloaded/i.test(detail)) {{
          return 'NVIDIA API: 503 Service temporarily overloaded. The selected AI model is temporarily unavailable; please try this cluster again later.';
        }}
        if (detail === 'Load failed' || detail === 'Failed to fetch') {{
          if (window.location.protocol === 'file:') {{
            return 'Open this report from the BGC-XPLORER Docker application so its AI controls can reach the built-in service.';
          }}
          return (
            'The connection to the BGC-XPLORER Docker service was interrupted after three attempts. ' +
            'Keep the report open, wait a moment, and try the AI analysis again. Browser detail: ' + detail
          );
        }}
        if (detail.toLowerCase().indexOf('timed out') !== -1) {{
          return (
            'The AI model took too long to answer. The request reached the service, but the upstream analysis timed out. ' +
            'Try the same cluster again; cached or shorter follow-up requests usually return faster. Browser detail: ' + detail
          );
        }}
        return detail;
      }}
      function resolveAiEndpoint(endpoint) {{
        try {{
          var url = new URL(endpoint, window.location.href);
          if ((url.hostname === '127.0.0.1' || url.hostname === 'localhost') && window.location.hostname) {{
            url.hostname = window.location.hostname;
            url.protocol = window.location.protocol === 'https:' ? 'https:' : 'http:';
          }}
          return url.toString();
        }} catch (error) {{
          return endpoint;
        }}
      }}
      function showActiveAiModel(model) {{
        if (!activeAiModel || !model) {{
          return;
        }}
        activeAiModel.textContent = model;
        if (aiModelLabel) {{
          aiModelLabel.textContent = 'Active AI model:';
        }}
      }}
      function loadActiveAiModel() {{
        if (!activeAiModel) {{
          return;
        }}
        var generationModel = activeAiModel.getAttribute('data-generation-model') || 'not-configured';
        var healthUrl = new URL('/ai_health', window.location.href).toString();
        fetch(healthUrl, {{headers: {{'Accept': 'application/json'}}}})
          .then(function (response) {{
            if (!response.ok) {{
              throw new Error('AI health unavailable');
            }}
            return response.json();
          }})
          .then(function (data) {{
            if (!data.model) {{
              throw new Error('AI model missing');
            }}
            showActiveAiModel(data.model);
          }})
          .catch(function () {{
            activeAiModel.textContent = generationModel;
            if (aiModelLabel) {{
              aiModelLabel.textContent = 'AI model at report generation:';
            }}
          }});
      }}
      loadActiveAiModel();
      function resetAiAnalysis() {{
        aiButton.disabled = false;
        aiButton.removeAttribute('data-sample');
        aiButton.removeAttribute('data-cluster');
        aiButton.setAttribute('aria-expanded', 'false');
        aiSpinner.hidden = true;
        aiStatus.textContent = '';
        aiResult.innerHTML = '';
      }}
      function setAiLoading(isLoading) {{
        aiButton.disabled = isLoading;
        aiSpinner.hidden = !isLoading;
        aiStatus.textContent = isLoading ? 'Analyzing...' : '';
      }}
      function closeViewer() {{
        viewer.hidden = true;
        document.body.classList.remove('gene-map-open');
        inline.innerHTML = '';
        interpretation.querySelector('p').textContent = '';
        if (geneTableBody) {{
          geneTableBody.innerHTML = '';
        }}
        if (geneTableContainer) {{
          geneTableContainer.hidden = true;
        }}
        resetAiAnalysis();
        image.hidden = true;
        image.removeAttribute('src');
        image.removeAttribute('data-file-src');
        image.removeAttribute('data-tried-fallback');
        error.hidden = true;
        error.textContent = '';
        if (activeButton) {{
          activeButton.setAttribute('aria-expanded', 'false');
          activeButton = null;
        }}
      }}
      image.addEventListener('error', function () {{
        var fallback = image.getAttribute('data-file-src') || '';
        var triedFallback = image.getAttribute('data-tried-fallback') === '1';
        if (fallback && !triedFallback) {{
          image.setAttribute('data-tried-fallback', '1');
          image.setAttribute('src', fallback);
          return;
        }}
        error.textContent = 'Gene map image could not be loaded from the report assets. Open the SVG directly: ' + (fallback || image.getAttribute('src') || '');
        error.hidden = false;
      }});
      image.addEventListener('load', function () {{
        error.hidden = true;
        error.textContent = '';
      }});
      document.addEventListener('click', function (event) {{
        var button = event.target.closest('.gene-map-btn');
        if (button) {{
          if (activeButton && activeButton !== button) {{
            activeButton.setAttribute('aria-expanded', 'false');
          }}
          activeButton = button;
          activeButton.setAttribute('aria-expanded', 'true');
          var cluster = button.getAttribute('data-cluster') || 'cluster';
          var clusterLabel = button.getAttribute('data-cluster-label') || cluster;
          var sample = button.getAttribute('data-sample') || '';
          var src = button.getAttribute('data-src') || '';
          var fileSrc = button.getAttribute('data-file-src') || '';
          var genes = button.getAttribute('data-genes') || '0';
          var categorySummary = button.getAttribute('data-summary') || 'No gene category summary available.';
          var biologicalInterpretation = button.getAttribute('data-interpretation') || 'No biological interpretation available for this cluster.';
          var aiCluster = button.getAttribute('data-ai-cluster') || '';
          var aiEnabled = (button.getAttribute('data-ai-enabled') || '') === '1';
          title.textContent = clusterLabel;
          meta.textContent = genes + ' genes. Hover over genes in the SVG for Bakta, eggNOG, dbCAN and ARTS details.';
          summary.textContent = categorySummary;
          interpretation.querySelector('p').textContent = biologicalInterpretation;
          renderGeneTable(cluster);
          resetAiAnalysis();
          aiButton.hidden = !aiEnabled;
          aiSpinner.hidden = true;
          aiButton.disabled = !aiEnabled;
          aiButton.setAttribute('data-sample', sample);
          aiButton.setAttribute('data-cluster', aiCluster);
          aiStatus.textContent = aiEnabled ? '' : 'AI analysis is available only for consensus BGC maps.';
          error.hidden = true;
          error.textContent = '';
          inline.innerHTML = '';
          image.hidden = true;
          image.removeAttribute('src');
          image.removeAttribute('data-file-src');
          image.removeAttribute('data-tried-fallback');
          if (svgData[cluster]) {{
            inline.innerHTML = svgData[cluster];
          }} else {{
            image.hidden = false;
            image.setAttribute('data-tried-fallback', '0');
            image.setAttribute('data-file-src', fileSrc);
            image.setAttribute('src', src);
            image.setAttribute('alt', cluster + ' gene map');
          }}
          viewer.hidden = false;
          document.body.classList.add('gene-map-open');
          return;
        }}
        var aiClick = event.target.closest('#gene-map-ai-button');
        if (aiClick) {{
          var endpoint = resolveAiEndpoint(viewer.getAttribute('data-ai-endpoint') || 'http://127.0.0.1:8787/analyze_cluster');
          var aiSample = aiClick.getAttribute('data-sample') || '';
          var aiCluster = aiClick.getAttribute('data-cluster') || '';
          if (!aiSample || !aiCluster) {{
            aiResult.innerHTML = '<div class="ai-status ai-status-error">Open a cluster gene map before requesting AI analysis.</div>';
            return;
          }}
          aiClick.setAttribute('aria-expanded', 'true');
          setAiLoading(true);
          aiResult.innerHTML = '<div class="ai-status">Sending gene data for ' + escapeHtml(aiCluster) + ' to the remote NVIDIA AI model. Clusters with more than 30 genes may require more time.</div>';
          fetchWithNetworkRetry(endpoint, {{
            method: 'POST',
            headers: {{'Content-Type': 'application/json', 'X-AI-Async': '1'}},
            body: JSON.stringify({{sample: aiSample, consensus_id: aiCluster}})
          }}, 3, 750)
            .then(function (response) {{
              return response.json().then(function (data) {{
                if ((response.status === 202 || data.status === 'queued' || data.status === 'running') && data.job_id) {{
                  return waitForAiJob(endpoint, data.job_id, 0);
                }}
                if (!response.ok) {{
                  throw new Error(diagnosticMessage(data));
                }}
                return data;
              }});
            }})
            .then(renderAnalysis)
            .catch(function (error) {{
              aiResult.innerHTML = '<div class="ai-status ai-status-error">' + escapeHtml(aiRequestErrorMessage(endpoint, error)) + '</div>';
            }})
            .finally(function () {{
              setAiLoading(false);
            }});
          return;
        }}
        if (event.target.closest('.gene-map-close')) {{
          closeViewer();
          return;
        }}
        if (!viewer.hidden && !event.target.closest('.gene-map-viewer')) {{
          closeViewer();
        }}
      }});
      document.addEventListener('keydown', function (event) {{
        if (event.key === 'Escape' && !viewer.hidden) {{
          closeViewer();
        }}
      }});
    }})();
  </script>
</body>
</html>
""".format(title=escape(title), blocks=blocks, ai_fetch_retry_js=AI_FETCH_RETRY_JS)
