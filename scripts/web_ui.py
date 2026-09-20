#!/usr/bin/env python3
"""BGC-XPLORER NiceGUI web interface.

Upload a bacterial genome FASTA, run the integrated BGC workflow, and review the report.
Served on port 8778 alongside the AI cluster server on 8787.
"""

import argparse
import asyncio
import csv
import hashlib
import logging
import os
import shutil
import subprocess
import threading
import urllib.error
import urllib.request
import uuid
from pathlib import Path

import yaml
from nicegui import ui, app
from fastapi import Request, Response

from fasta_input import FASTA_EXTENSIONS, MAX_FASTA_BYTES, fasta_suffix, normalize_sample_name
from homepage_content import HERO_BODY, HERO_TITLE, RESULT_FEATURES, WORKFLOW_STEPS, section_header
from report_branding import image_data_uri
from workflow_progress import progress_value, unread_lines
from thread_config import configured_threads

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s %(message)s")
app_logger = logging.getLogger("bgc_xplorer")

WORK_DIR = Path("/work")
APP_DIR = Path("/app")
APP_LOGO_DATA_URI = image_data_uri((APP_DIR / "logo.jpg", Path.cwd() / "logo.jpg"))

workflow_state: dict = {}
AI_SERVICE_URL = os.environ.get("AI_SERVICE_URL", "http://127.0.0.1:8484")
AI_JOBS: dict = {}
AI_JOBS_LOCK = threading.Lock()

THEME_CSS = """
<style>
  :root {
    --bg: #f5f7ff;
    --bg-soft: #fbfcff;
    --panel: rgba(255, 255, 255, 0.94);
    --panel-2: #ffffff;
    --panel-3: #eef2ff;
    --line: rgba(94, 108, 152, 0.12);
    --line-strong: rgba(94, 108, 152, 0.2);
    --text: #273149;
    --text-strong: #161f33;
    --muted: #5e6c8f;
    --accent: #6c63f6;
    --accent-soft: rgba(108, 99, 246, 0.12);
    --accent-deep: #5b52eb;
    --violet: #6c63f6;
    --emerald: #73c95b;
    --amber: #f59b38;
    --teal: #41bcc7;
    --blue: #4f8dfd;
    --orange: #ff7b1d;
    --sky: #4a7fff;
    --lime: #7cc85d;
    --cyan: #35bfd3;
    --indigo: #7c6bff;
    --radius: 20px;
    --radius-sm: 14px;
    --shadow-card: 0 12px 35px rgba(89, 104, 146, 0.08);
    --shadow-soft: 0 8px 24px rgba(89, 104, 146, 0.06);
    --sans: "Sora", "Avenir Next", "Segoe UI", sans-serif;
  }
  body, .nicegui-content {
    font-family: var(--sans);
    color: var(--text);
    background:
      radial-gradient(900px 440px at 8% 0%, rgba(95, 87, 255, 0.09), transparent 62%),
      radial-gradient(920px 480px at 100% 16%, rgba(53, 191, 211, 0.09), transparent 58%),
      linear-gradient(180deg, #fcfdff 0%, #f3f6ff 100%),
      var(--bg);
  }
  .nicegui-content,
  .q-page,
  .q-page-container,
  .q-layout__section--main {
    width: 100%;
  }
  .app-shell {
    width: min(1680px, calc(100vw - 56px));
    max-width: none;
    margin: 0 auto;
    padding: 32px 28px 48px;
    box-sizing: border-box;
    align-self: center;
  }
  .hero-card {
    display: flex;
    justify-content: space-between;
    gap: 22px;
    align-items: flex-start;
    width: 100%;
    margin-bottom: 28px;
  }
  .hero-brand {
    display: flex;
    gap: 18px;
    align-items: flex-start;
    min-width: 0;
    flex: 1;
  }
  .brand-mark {
    width: 56px;
    height: 56px;
    border-radius: 16px;
    border: 1px solid rgba(95, 87, 255, 0.28);
    background: linear-gradient(180deg, rgba(255,255,255,0.98), rgba(241,244,255,0.96));
    display: inline-flex;
    align-items: center;
    justify-content: center;
    color: var(--accent);
    font-size: 26px;
    box-shadow: var(--shadow-soft);
    flex: 0 0 auto;
  }
  .app-brand-logo {
    width: clamp(180px, 18vw, 250px);
    height: auto;
    border-radius: 14px;
    object-fit: contain;
    flex: 0 0 auto;
  }
  .hero-copy {
    min-width: 0;
  }
  .upload-hero { padding: 12px 0 8px; }
  .upload-hero .hero-brand {
    display: grid;
    grid-template-columns: clamp(240px, 26vw, 360px) minmax(0, 1fr);
    align-items: center;
    gap: 28px;
    width: 100%;
  }
  .upload-hero .app-brand-logo {
    display: block;
    width: 100%;
    max-width: 100%;
    height: auto;
    object-fit: contain;
  }
  .upload-hero .hero-copy { width: 100%; }
  .upload-hero .hero-title { text-wrap: balance; }
  .upload-hero .hero-text { margin: 16px 0 0; }
  @media (max-width: 720px) {
    .upload-hero .hero-brand { grid-template-columns: 1fr; gap: 20px; }
    .upload-hero .app-brand-logo { max-width: 280px; }
  }
  .eyebrow {
    text-transform: uppercase;
    letter-spacing: 0;
    font-size: 15px;
    font-weight: 700;
    color: var(--accent);
    margin-bottom: 10px;
  }
  .hero-title {
    font-size: clamp(34px, 4vw, 52px);
    line-height: 1.02;
    margin: 0;
    color: var(--text-strong);
  }
  .hero-text {
    margin-top: 16px;
    max-width: 920px;
    color: var(--muted);
    font-size: 16px;
    line-height: 1.7;
  }
  .meta-badge {
    min-width: 210px;
    padding: 18px 20px;
    border-radius: 18px;
    border: 1px solid var(--line);
    background: rgba(255, 255, 255, 0.78);
    box-shadow: var(--shadow-soft);
    color: var(--muted);
  }
  .meta-badge-label {
    font-size: 13px;
    color: var(--muted);
    margin-bottom: 8px;
  }
  .meta-badge-value {
    font-weight: 700;
    color: var(--text-strong);
    line-height: 1.4;
  }
  .panel {
    width: 100%;
    background: var(--panel);
    border: 1px solid var(--line);
    border-radius: 22px;
    box-shadow: var(--shadow-card);
    padding: 24px;
    backdrop-filter: blur(12px);
  }
  .section-head {
    display: flex;
    align-items: center;
    justify-content: space-between;
    gap: 16px;
    margin-bottom: 18px;
  }
  .section-head h2 {
    margin: 0;
    font-size: 30px;
    color: var(--text-strong);
  }
  .section-accent {
    width: 28px;
    height: 4px;
    border-radius: 999px;
    background: linear-gradient(90deg, var(--accent), rgba(95, 87, 255, 0.15));
  }
  .section-description {
    margin: -8px 0 18px;
  }
  .metrics-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
    gap: 16px;
    width: 100%;
  }
  .metrics-grid-4 {
    grid-template-columns: repeat(4, minmax(0, 1fr));
  }
  .metric-card {
    display: flex;
    align-items: flex-start;
    gap: 14px;
    min-height: 150px;
    padding: 18px 18px 16px;
    border-radius: 18px;
    border: 1px solid var(--line);
    background: linear-gradient(180deg, rgba(255,255,255,0.98), rgba(244,247,255,0.96));
  }
  .metric-icon {
    width: 48px;
    height: 48px;
    border-radius: 14px;
    display: inline-flex;
    align-items: center;
    justify-content: center;
    font-size: 21px;
    flex: 0 0 auto;
    border: 1px solid transparent;
  }
  .metric-copy {
    min-width: 0;
  }
  .metric-label {
    font-size: 12px;
    font-weight: 700;
    letter-spacing: 0;
    text-transform: uppercase;
    color: var(--text-strong);
    margin-bottom: 8px;
  }
  .metric-value {
    font-size: 22px;
    font-weight: 700;
    color: var(--text-strong);
    line-height: 1;
    margin-bottom: 10px;
  }
  .metric-copy p {
    margin: 0;
    font-size: 13px;
    color: var(--muted);
    line-height: 1.55;
  }
  .metric-violet .metric-icon { color: var(--violet); background: rgba(106,92,255,0.1); border-color: rgba(106,92,255,0.18); }
  .metric-green .metric-icon { color: var(--emerald); background: rgba(115,201,91,0.12); border-color: rgba(115,201,91,0.18); }
  .metric-blue .metric-icon { color: var(--blue); background: rgba(79,141,253,0.1); border-color: rgba(79,141,253,0.18); }
  .metric-orange .metric-icon { color: var(--orange); background: rgba(255,123,29,0.1); border-color: rgba(255,123,29,0.18); }
  .metric-teal .metric-icon { color: var(--teal); background: rgba(65,188,199,0.1); border-color: rgba(65,188,199,0.18); }
  .metric-indigo .metric-icon { color: var(--indigo); background: rgba(124,107,255,0.1); border-color: rgba(124,107,255,0.18); }
  .metric-lime .metric-icon { color: var(--lime); background: rgba(124,200,93,0.12); border-color: rgba(124,200,93,0.18); }
  .metric-sky .metric-icon { color: var(--sky); background: rgba(74,127,255,0.1); border-color: rgba(74,127,255,0.18); }
  .metric-amber .metric-icon { color: var(--amber); background: rgba(245,155,56,0.12); border-color: rgba(245,155,56,0.18); }
  .metric-cyan .metric-icon { color: var(--cyan); background: rgba(53,191,211,0.1); border-color: rgba(53,191,211,0.18); }
  .info-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(360px, 1fr));
    gap: 18px;
    width: 100%;
  }
  .info-card {
    display: flex;
    gap: 16px;
    align-items: flex-start;
  }
  .info-icon {
    width: 54px;
    height: 54px;
    border-radius: 16px;
    display: inline-flex;
    align-items: center;
    justify-content: center;
    background: linear-gradient(180deg, rgba(95,87,255,0.16), rgba(95,87,255,0.06));
    color: var(--accent);
    font-size: 24px;
    flex: 0 0 auto;
  }
  .info-copy h3 {
    margin: 0 0 8px 0;
    font-size: 24px;
    color: var(--text-strong);
  }
  .info-copy p, .muted {
    color: var(--muted);
    line-height: 1.65;
  }
  .upload-drop {
    border-radius: 18px;
    border: 1.5px dashed rgba(108, 99, 246, 0.12);
    background: linear-gradient(180deg, rgba(108,99,246,0.02), rgba(255,255,255,0.84));
  }
  .analysis-launch-panel {
    margin-bottom: 28px;
    border-color: rgba(108, 99, 246, 0.2);
  }
  .workflow-grid {
    display: grid;
    grid-template-columns: repeat(5, minmax(0, 1fr));
    gap: 14px;
  }
  .workflow-step,
  .feature-card {
    border: 1px solid var(--line);
    border-radius: 18px;
    background: linear-gradient(180deg, rgba(255,255,255,0.98), rgba(244,247,255,0.96));
    padding: 18px;
  }
  .workflow-step-number {
    width: 34px;
    height: 34px;
    display: grid;
    place-items: center;
    border-radius: 10px;
    background: var(--accent-soft);
    color: var(--accent-deep);
    font-weight: 700;
    margin-bottom: 15px;
  }
  .workflow-step h3,
  .feature-card h3 {
    margin: 0 0 9px;
    color: var(--text-strong);
    font-size: 17px;
  }
  .workflow-step p,
  .feature-card p {
    margin: 0;
    color: var(--muted);
    font-size: 13px;
    line-height: 1.6;
  }
  .feature-grid {
    display: grid;
    grid-template-columns: repeat(3, minmax(0, 1fr));
    gap: 14px;
  }
  .feature-icon {
    color: var(--accent);
    font-size: 23px;
    margin-bottom: 12px;
  }
  .status-chip {
    display: inline-flex;
    align-items: center;
    min-height: 34px;
    padding: 7px 12px;
    border-radius: 999px;
    border: 1px solid var(--line);
    background: rgba(255, 255, 255, 0.78);
    color: var(--muted);
    font-size: 13px;
    font-weight: 600;
  }
  .status-chip.ready {
    color: #227857;
    background: rgba(115,201,91,0.12);
    border-color: rgba(115,201,91,0.25);
  }
  .cta-row {
    display: flex;
    gap: 12px;
    flex-wrap: wrap;
    align-items: center;
  }
  .previous-results {
    display: grid;
    grid-template-columns: repeat(3, minmax(0, 1fr));
    gap: 14px;
    width: 100%;
  }
  @media (max-width: 1120px) {
    .previous-results { grid-template-columns: repeat(2, minmax(0, 1fr)); }
  }
  @media (max-width: 720px) {
    .previous-results { grid-template-columns: minmax(0, 1fr); }
  }
  .result-card {
    display: grid;
    grid-template-columns: 48px minmax(0, 1fr);
    align-content: start;
    column-gap: 14px;
    row-gap: 12px;
    width: 100%;
    min-height: 190px;
    padding: 18px 18px 16px;
    border-radius: 18px;
    border: 1px solid var(--line);
    background: linear-gradient(180deg, rgba(255,255,255,0.98), rgba(244,247,255,0.96));
  }
  .result-card-copy {
    min-width: 0;
  }
  .result-card .metric-icon {
    width: 48px;
    height: 48px;
    font-size: 21px;
  }
  .result-card-title {
    font-size: 12px;
    font-weight: 700;
    letter-spacing: 0;
    text-transform: uppercase;
    color: var(--text-strong);
    margin-bottom: 8px;
  }
  .result-card-value {
    font-size: 22px;
    font-weight: 700;
    color: var(--text-strong);
    line-height: 1;
    margin-bottom: 10px;
  }
  .result-card-meta {
    font-size: 13px;
    color: var(--muted);
    line-height: 1.55;
  }
  .result-link {
    display: inline-flex;
    align-items: center;
    gap: 10px;
    padding: 10px 14px;
    border-radius: 14px;
    border: 1px solid var(--line);
    background: rgba(255,255,255,0.82);
    text-decoration: none;
    color: var(--accent-deep);
    font-weight: 600;
  }
  .result-actions {
    display: flex;
    flex-wrap: wrap;
    gap: 10px;
    margin-top: 2px;
    grid-column: 2;
  }
  .result-link:hover {
    border-color: rgba(95,87,255,0.28);
    transform: translateY(-1px);
  }
  .progress-shell {
    display: grid;
    gap: 18px;
    width: 100%;
  }
  .log-shell {
    border-radius: 18px;
    overflow: hidden;
    border: 1px solid rgba(31, 41, 55, 0.14);
    background: #0f172a;
  }
  .iframe-shell {
    width: 100%;
    border-radius: 22px;
    overflow: hidden;
    border: 1px solid var(--line);
    box-shadow: var(--shadow-card);
    background: #fff;
  }
  .summary-links {
    display: flex;
    flex-wrap: wrap;
    gap: 10px;
  }
  .summary-links a {
    display: inline-flex;
    align-items: center;
    padding: 10px 14px;
    border-radius: 999px;
    border: 1px solid var(--line);
    color: var(--accent-deep);
    background: rgba(255,255,255,0.84);
    text-decoration: none;
    font-weight: 600;
  }
  .q-btn {
    border-radius: 14px;
    text-transform: none;
    letter-spacing: 0;
    box-shadow: none;
  }
  .q-btn.bg-primary {
    background: linear-gradient(135deg, var(--accent) 0%, var(--accent-deep) 100%) !important;
  }
  .q-uploader {
    background: transparent !important;
  }
  .q-uploader__header {
    background: linear-gradient(135deg, rgba(108,99,246,0.62) 0%, rgba(91,82,235,0.7) 100%) !important;
  }
  .q-uploader__header-content,
  .q-uploader__title,
  .q-uploader__subtitle {
    color: #ffffff !important;
  }
  .q-uploader__list {
    background: transparent !important;
  }
  .soft-primary-btn {
    background: linear-gradient(135deg, rgba(108,99,246,0.62) 0%, rgba(91,82,235,0.7) 100%) !important;
    color: #ffffff !important;
  }
  .q-linear-progress {
    border-radius: 999px;
    overflow: hidden;
    height: 12px;
    background: rgba(95,87,255,0.08);
  }
  @media (max-width: 900px) {
    .app-shell {
      width: 100%;
      padding: 22px 18px 32px;
    }
    .hero-card {
      flex-direction: column;
    }
    .meta-badge {
      width: 100%;
      min-width: 0;
    }
  }
  @media (max-width: 1280px) {
    .metrics-grid-4 {
      grid-template-columns: repeat(2, minmax(0, 1fr));
    }
    .workflow-grid {
      grid-template-columns: repeat(2, minmax(0, 1fr));
    }
  }
  @media (max-width: 720px) {
    .metrics-grid-4 {
      grid-template-columns: 1fr;
    }
    .workflow-grid,
    .feature-grid {
      grid-template-columns: 1fr;
    }
  }
</style>
"""

STAT_TONES = {
    "Upload set": "violet",
    "Required files": "green",
    "Ready samples": "blue",
    "Workflow": "orange",
    "Total loci": "indigo",
    "Multi-caller loci": "lime",
    "MIBiG comparisons": "cyan",
    "ARTS known hits": "orange",
    "ARTS DUF hits": "teal",
    "antiSMASH": "violet",
    "GECCO": "green",
    "DeepBGC": "blue",
    "ARTS": "orange",
    "dbCAN CGC": "teal",
    "MIBiG hits": "cyan",
}

STAT_GLYPHS = {
    "Upload set": "⌂",
    "Required files": "✓",
    "Ready samples": "◌",
    "Workflow": "↺",
    "Total loci": "◔",
    "Multi-caller loci": "◎",
    "MIBiG comparisons": "⬡",
    "ARTS known hits": "⛨",
    "ARTS DUF hits": "◇",
    "antiSMASH": "○",
    "GECCO": "◇",
    "DeepBGC": "□",
    "ARTS": "⛨",
    "dbCAN CGC": "⌘",
    "MIBiG hits": "◌",
}


def inject_theme():
    ui.add_head_html(THEME_CSS)


def stat_card(label: str, value: str, note: str) -> str:
    tone = STAT_TONES.get(label, "violet")
    glyph = STAT_GLYPHS.get(label, "•")
    return (
        "<article class='metric-card metric-{tone}'>"
        "<div class='metric-icon' aria-hidden='true'>{glyph}</div>"
        "<div class='metric-copy'>"
        "<div class='metric-label'>{label}</div>"
        "<div class='metric-value'>{value}</div>"
        "<p>{note}</p>"
        "</div>"
        "</article>"
    ).format(tone=tone, glyph=glyph, label=label, value=value, note=note)


def hero_section(title: str, body: str, generated: str | None = None, upload: bool = False) -> str:
    meta = ""
    if generated:
        meta = (
            "<div class='meta-badge'>"
            "<div class='meta-badge-value'>{generated}</div>"
            "</div>"
        ).format(generated=generated)
    brand = (
        "<img class='app-brand-logo' src='{0}' alt='BGC-XPLORER' />".format(APP_LOGO_DATA_URI)
        if APP_LOGO_DATA_URI
        else "<div class='brand-mark' aria-hidden='true'>⌬</div>"
    )
    return (
        "<section class='hero-card{variant}'>"
        "<div class='hero-brand'>"
        "{brand}"
        "<div class='hero-copy'>"
        "<div class='eyebrow'>BGC-XPLORER</div>"
        "<h1 class='hero-title'>{title}</h1>"
        "<p class='hero-text'>{body}</p>"
        "</div>"
        "</div>"
        "{meta}"
        "</section>"
    ).format(title=title, body=body, meta=meta, brand=brand, variant=" upload-hero" if upload else "")


def info_card(title: str, body: str, icon: str = "i") -> str:
    return (
        "<article class='info-card'>"
        "<div class='info-icon' aria-hidden='true'>{icon}</div>"
        "<div class='info-copy'>"
        "<h3>{title}</h3>"
        "<p>{body}</p>"
        "</div>"
        "</article>"
    ).format(title=title, body=body, icon=icon)


def workflow_step_card(title: str, body: str, number: str) -> str:
    return (
        "<article class='workflow-step'>"
        "<div class='workflow-step-number'>{number}</div>"
        "<h3>{title}</h3><p>{body}</p>"
        "</article>"
    ).format(title=title, body=body, number=number)


def result_feature_card(title: str, body: str, icon: str) -> str:
    return (
        "<article class='feature-card'>"
        "<div class='feature-icon' aria-hidden='true'>{icon}</div>"
        "<h3>{title}</h3><p>{body}</p>"
        "</article>"
    ).format(title=title, body=body, icon=icon)


def ai_proxy_headers():
    return {
        "Access-Control-Allow-Origin": "*",
        "Access-Control-Allow-Headers": "Content-Type, X-AI-Async",
        "Access-Control-Allow-Methods": "GET, POST, OPTIONS",
    }


@app.options("/analyze_cluster")
async def proxy_analyze_cluster_options():
    return Response(status_code=204, headers=ai_proxy_headers())


@app.get("/analyze_cluster")
async def proxy_analyze_cluster_info():
    return Response(
        content='{"ok": true, "message": "AI analysis endpoint is ready. Use POST with sample and consensus_id."}',
        status_code=200,
        media_type="application/json",
        headers=ai_proxy_headers(),
    )


def forward_ai_request(body, content_type, timeout):
    upstream = urllib.request.Request(
        f"{AI_SERVICE_URL.rstrip('/')}/analyze_cluster",
        data=body,
        headers={"Content-Type": content_type},
        method="POST",
    )
    with urllib.request.urlopen(upstream, timeout=timeout) as response:
        return response.status, response.headers.get("Content-Type", "application/json"), response.read()


def ai_job_id(body):
    return hashlib.sha256(body).hexdigest()[:24]


def run_ai_job(job_id, body, content_type):
    with AI_JOBS_LOCK:
        AI_JOBS[job_id]["status"] = "running"
    try:
        status_code, media_type, content = forward_ai_request(body, content_type, 270)
        with AI_JOBS_LOCK:
            AI_JOBS[job_id].update(
                {
                    "status": "done",
                    "status_code": status_code,
                    "media_type": media_type,
                    "content": content,
                }
            )
    except urllib.error.HTTPError as error:
        with AI_JOBS_LOCK:
            AI_JOBS[job_id].update(
                {
                    "status": "error",
                    "status_code": error.code,
                    "media_type": error.headers.get("Content-Type", "application/json"),
                    "content": error.read(),
                }
            )
    except Exception as error:
        with AI_JOBS_LOCK:
            AI_JOBS[job_id].update(
                {
                    "status": "error",
                    "status_code": 503,
                    "media_type": "application/json",
                    "content": ('{"error": "AI service proxy failed: %s"}' % str(error).replace('"', "'")).encode(),
                }
            )


def start_ai_job(body, content_type):
    job_id = ai_job_id(body)
    with AI_JOBS_LOCK:
        existing = AI_JOBS.get(job_id)
        if existing and existing.get("status") in {"queued", "running", "done"}:
            return job_id, existing.get("status", "queued")
        AI_JOBS[job_id] = {
            "status": "queued",
            "status_code": 202,
            "media_type": "application/json",
            "content": b"",
        }
    thread = threading.Thread(target=run_ai_job, args=(job_id, body, content_type), daemon=True)
    thread.start()
    return job_id, "queued"


@app.post("/analyze_cluster")
async def proxy_analyze_cluster(request: Request):
    body = await request.body()
    content_type = request.headers.get("content-type", "application/json")
    if request.headers.get("x-ai-async") == "1":
        job_id, status = start_ai_job(body, content_type)
        return Response(
            content='{"status": "%s", "job_id": "%s"}' % (status, job_id),
            status_code=202,
            media_type="application/json",
            headers=ai_proxy_headers(),
        )
    try:
        status_code, media_type, content = await asyncio.to_thread(
            forward_ai_request,
            body,
            content_type,
            180,
        )
        return Response(
            content=content,
            status_code=status_code,
            media_type=media_type,
            headers=ai_proxy_headers(),
        )
    except urllib.error.HTTPError as error:
        return Response(
            content=error.read(),
            status_code=error.code,
            media_type=error.headers.get("Content-Type", "application/json"),
            headers=ai_proxy_headers(),
        )
    except Exception as error:
        return Response(
            content='{"error": "AI service proxy failed: %s"}' % str(error).replace('"', "'"),
            status_code=503,
            media_type="application/json",
            headers=ai_proxy_headers(),
        )


@app.get("/analyze_cluster_status/{job_id}")
async def proxy_analyze_cluster_status(job_id: str):
    with AI_JOBS_LOCK:
        job = AI_JOBS.get(job_id)
        if not job:
            return Response(
                content='{"error": "AI analysis job was not found."}',
                status_code=404,
                media_type="application/json",
                headers=ai_proxy_headers(),
            )
        status = job.get("status", "queued")
        if status in {"queued", "running"}:
            return Response(
                content='{"status": "%s", "job_id": "%s"}' % (status, job_id),
                status_code=202,
                media_type="application/json",
                headers=ai_proxy_headers(),
            )
        return Response(
            content=job.get("content", b""),
            status_code=job.get("status_code", 200),
            media_type=job.get("media_type", "application/json"),
            headers=ai_proxy_headers(),
        )


@app.get("/analyze_cluster_ping")
async def proxy_analyze_cluster_ping():
    return Response(
        content='{"ok": true, "proxy": "ready", "upstream": "%s"}' % AI_SERVICE_URL.rstrip("/"),
        status_code=200,
        media_type="application/json",
        headers=ai_proxy_headers(),
    )


@app.get("/ai_health")
async def proxy_ai_health():
    try:
        with urllib.request.urlopen(f"{AI_SERVICE_URL.rstrip('/')}/health", timeout=5) as response:
            return Response(
                content=response.read(),
                status_code=response.status,
                media_type=response.headers.get("Content-Type", "application/json"),
                headers=ai_proxy_headers(),
            )
    except Exception as error:
        return Response(
            content='{"error": "AI service health check failed: %s"}' % str(error).replace('"', "'"),
            status_code=503,
            media_type="application/json",
            headers=ai_proxy_headers(),
        )


def read_tsv_rows(path: Path):
    if not path.exists():
        return []
    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def summarize_sample(sample: str):
    summary_dir = WORK_DIR / "results" / sample / "summary"
    evidence = read_tsv_rows(summary_dir / "region_evidence.tsv")
    arts = read_tsv_rows(summary_dir / "arts.hits.tsv")
    dbcan = read_tsv_rows(summary_dir / "dbcan.cgc.tsv")
    mibig = read_tsv_rows(summary_dir / "mibig_dereplication.tsv")
    antismash = read_tsv_rows(summary_dir / "antismash.bgc.tsv")
    gecco = read_tsv_rows(summary_dir / "gecco.bgc.tsv")
    deepbgc = read_tsv_rows(summary_dir / "deepbgc.bgc.tsv")

    def nonempty(value):
        text = str(value or "").strip()
        return bool(text) and text.lower() not in {"nan", "n/a", "none"}

    return {
        "consensus": len(evidence),
        "multitool": sum(1 for row in evidence if int(float(row.get("support_count") or 0)) >= 2),
        "arts_known_loci": sum(1 for row in evidence if int(float(row.get("arts_known_hits") or 0)) > 0),
        "arts_duf_loci": sum(1 for row in evidence if int(float(row.get("arts_duf_hits") or 0)) > 0),
        "mibig_hits": sum(1 for row in mibig if nonempty(row.get("best_mibig_id"))),
        "antismash": len(antismash),
        "gecco": len(gecco),
        "deepbgc": len(deepbgc),
        "arts": len(arts),
        "dbcan": len(dbcan),
    }


# ─── Upload Page ──────────────────────────────────────────────────────────────

@ui.page("/")
async def upload_page():
    ui.page_title("BGC-XPLORER")
    inject_theme()

    state = {"uploaded_path": None}
    upload_token = str(uuid.uuid4())
    upload_dir = WORK_DIR / "data" / "uploads" / upload_token
    upload_dir.mkdir(parents=True, exist_ok=True)

    results_dir = WORK_DIR / "results"
    existing = []
    if results_dir.exists():
        existing = [
            d.name for d in sorted(results_dir.iterdir())
            if d.is_dir()
            and (d / "report" / f"{d.name}.html").exists()
            and d.name != "example_sample"
        ]

    with ui.column().classes("app-shell"):
        ui.html(
            hero_section(
                HERO_TITLE,
                HERO_BODY,
                upload=True,
            )
        )

        with ui.card().classes("panel analysis-launch-panel w-full"):
            with ui.column().classes("w-full gap-4"):
                ui.html(
                    section_header(
                        "Start a new analysis",
                        "Provide one assembled bacterial genome, MAG, or plasmid. "
                        "BGC-XPLORER will run the complete discovery, biological-context, "
                        "evidence aggregation, and reporting workflow.",
                    )
                )
                sample_input = ui.input(
                    "Sample name", placeholder="e.g. Soil_1"
                ).classes("w-full")
                sample_input.props("outlined standout")

                ui.separator()
                ui.label("Genome FASTA").classes("text-2xl font-semibold")
                ui.label(
                    "One uncompressed nucleotide file with extension .fa, .fasta, or .fna (maximum 30 MB)."
                ).classes("muted text-sm")

                status_label = ui.label("\u2399 Waiting for FASTA").classes("status-chip")

                async def handle_upload(e):
                    file = e.file
                    filename = file.name
                    ext = fasta_suffix(filename)
                    if not ext:
                        ui.notify(f"Skipped {filename}: use {', '.join(FASTA_EXTENSIONS)}.", type="warning")
                        return
                    content = await file.read()
                    if len(content) > MAX_FASTA_BYTES:
                        ui.notify(f"{filename} is larger than 30 MB. Choose a smaller FASTA file.", type="warning")
                        return
                    stored_path = upload_dir / filename
                    stored_path.write_bytes(content)
                    previous = state.get("uploaded_path")
                    if previous and Path(previous) != stored_path:
                        Path(previous).unlink(missing_ok=True)
                    state["uploaded_path"] = str(stored_path)
                    status_label.set_text(f"\u2705 {filename}")
                    status_label.classes(replace="status-chip ready")
                    ui.notify(f"Uploaded: {filename}", type="positive")

                ui.upload(
                    on_upload=handle_upload,
                    multiple=False,
                    auto_upload=True,
                    max_files=1,
                    max_file_size=MAX_FASTA_BYTES,
                    max_total_size=MAX_FASTA_BYTES,
                ).classes("w-full upload-drop").props(
                    'accept=".fa,.fasta,.fna"'
                )

                async def start_analysis():
                    try:
                        sample = normalize_sample_name(sample_input.value)
                    except ValueError:
                        ui.notify("Please enter a sample name", type="warning")
                        return
                    uploaded_path = state.get("uploaded_path")
                    if not uploaded_path or not Path(uploaded_path).is_file():
                        ui.notify("Please upload one FASTA file", type="warning")
                        return
                    if Path(uploaded_path).stat().st_size > MAX_FASTA_BYTES:
                        ui.notify("FASTA file is larger than 30 MB. Choose a smaller file.", type="warning")
                        return

                    fasta_dir = WORK_DIR / "data" / "fasta"
                    fasta_dir.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(uploaded_path, fasta_dir / f"{sample}.fasta")

                    update_samples_tsv(sample, default_taxon())
                    ui.notify(f"Starting workflow for '{sample}'...", type="positive")
                    state["uploaded_path"] = None
                    shutil.rmtree(str(upload_dir), ignore_errors=True)
                    ui.navigate.to(f"/progress/{sample}")

                with ui.row().classes("cta-row mt-4"):
                    ui.button("Start Analysis", on_click=start_analysis).props(
                        "size=lg unelevated"
                    ).classes("soft-primary-btn").style(
                        "background: linear-gradient(135deg, rgba(108,99,246,0.62) 0%, rgba(91,82,235,0.7) 100%) !important; "
                        "background-color: #7b73f6 !important; color: #ffffff !important; border: none !important;"
                    )

        ui.html(
            "<section class='panel'>"
            "{heading}"
            "<div class='workflow-grid'>{steps}</div>"
            "</section>".format(
                heading=section_header(
                    "From sequence to an integrated evidence report",
                    "Each stage contributes independent evidence to the final integrated result.",
                ),
                steps="".join(workflow_step_card(*step) for step in WORKFLOW_STEPS)
            )
        )

        ui.html(
            "<section class='panel'>"
            "{heading}"
            "<div class='feature-grid'>{features}</div>"
            "</section>".format(
                heading=section_header(
                    "What you get",
                    "A report designed for candidate review, comparison, and reproducible scientific use.",
                ),
                features="".join(result_feature_card(*feature) for feature in RESULT_FEATURES)
            )
        )

        if existing:
            with ui.card().classes("panel w-full"):
                with ui.column().classes("w-full gap-4"):
                    ui.html(
                        "<div class='section-head'><h2>Previous results</h2><span class='section-accent'></span></div>"
                    )
                    ui.html(
                        "<div class='previous-results'>" + "".join(
                                "<article class='result-card metric-indigo'>"
                                "<div class='metric-icon' aria-hidden='true'>◔</div>"
                                "<div class='result-card-copy'>"
                                "<div class='result-card-title'>Ready report</div>"
                                f"<div class='result-card-value'>{s}</div>"
                                "<div class='result-card-meta'>Generated report, tables and exported HTML assets are ready to open.</div>"
                                "</div>"
                                "<div class='result-actions'>"
                                f"<a class='result-link' href='/static_results/{s}/report/{s}.html' target='_blank'>Open report</a>"
                                "</div>"
                                "</article>"
                            for s in existing
                        ) + "</div>"
                    )


# ─── Progress Page ────────────────────────────────────────────────────────────

def start_workflow_process(sample: str):
    thread_count = configured_threads()
    cmd = [
        "/opt/conda/bin/snakemake",
        "-s", str(APP_DIR / "Snakefile"),
        "--configfile", str(WORK_DIR / "config" / "config.yaml"),
        "--cores", str(thread_count),
        "--directory", str(WORK_DIR),
        "--printshellcmds",
        "--rerun-incomplete",
        f"results/{sample}/report/{sample}.html",
    ]
    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        cwd=str(WORK_DIR),
        text=True,
    )
    state = {
        "process": process,
        "log": [],
        "status": "running",
    }
    workflow_state[sample] = state

    def reader():
        for line in process.stdout:
            workflow_state[sample]["log"].append(line.rstrip())
        process.wait()
        workflow_state[sample]["status"] = (
            "done" if process.returncode == 0 else "error"
        )

    t = threading.Thread(target=reader, daemon=True)
    t.start()
    state["thread"] = t


@ui.page("/progress/{sample}")
async def progress_page(sample: str):
    ui.page_title(f"BGC-XPLORER — {sample}")
    inject_theme()

    if sample not in workflow_state:
        start_workflow_process(sample)

    state = workflow_state[sample]

    with ui.column().classes("app-shell"):
        ui.html(
            hero_section(
                "Analyzing sample {0}".format(sample),
                "The workflow is running the detection, merging and report-generation steps. You can stay here and watch the live log update.",
                "Running now",
            )
        )

        with ui.row().classes("w-full info-grid mb-6"):
            ui.html(
                info_card(
                    "Live progress",
                    "This page follows the workflow state in real time and turns into a launch point for the finished report as soon as the job completes.",
                    "↺",
                )
            )
            ui.html(
                info_card(
                    "What to watch",
                    "If something fails, the log below is the first place to inspect. Otherwise you can simply wait for the result button to appear.",
                    "⌘",
                )
            )

        with ui.card().classes("panel w-full progress-shell"):
            with ui.row().classes("items-center justify-between w-full"):
                ui.html("<div class='section-head' style='margin:0'><h2>Workflow progress</h2><span class='section-accent'></span></div>")
                ui.button("New Upload", on_click=lambda: ui.navigate.to("/")).props(
                    "outline color=primary"
                )

            status_label = ui.label("Running...").classes("text-lg font-semibold")
            with ui.row().classes("items-center w-full gap-3"):
                progress_bar = ui.linear_progress(value=0).classes("grow")
                progress_label = ui.label("0%").classes("text-lg font-semibold min-w-12 text-right")
            report_url = f"/static_results/{sample}/report/{sample}.html"
            result_btn = ui.button("View Results", on_click=lambda: ui.navigate.to(
                report_url
            )).props("color=primary size=lg unelevated")
            result_btn.set_visibility(False)

        with ui.expansion("Workflow log", icon="terminal", value=False).classes("panel w-full"):
            ui.label("Live combined stdout and stderr from the Snakemake run.").classes("muted text-sm")
            log_widget = ui.log(max_lines=5000).classes(
                "w-full h-96 bg-gray-900 text-green-400 font-mono text-xs log-shell"
            )

        log_cursor = 0

        async def poll():
            nonlocal log_cursor
            new_lines, log_cursor = unread_lines(state["log"], log_cursor)
            for line in new_lines:
                log_widget.push(line)

            if state["status"] == "done":
                status_label.set_text("Analysis complete!")
                status_label.classes(replace="text-green-400 text-lg")
                progress_bar.set_value(1.0)
                progress_label.set_text("100%")
                result_btn.set_visibility(True)
            elif state["status"] == "error":
                status_label.set_text("Analysis failed. Check log below.")
                status_label.classes(replace="text-red-400 text-lg")
                progress_bar.set_value(1.0)

            if state["status"] == "running":
                value, label = progress_value(state["log"], total=15)
                progress_bar.set_value(value)
                progress_label.set_text(label)

        ui.timer(0.5, poll)


# ─── Results Page ─────────────────────────────────────────────────────────────

@ui.page("/results/{sample}")
async def results_page(sample: str):
    ui.page_title(f"BGC-XPLORER — {sample}")
    inject_theme()

    report_path = WORK_DIR / "results" / sample / "report" / f"{sample}.html"
    summary_dir = WORK_DIR / "results" / sample / "summary"
    metrics = summarize_sample(sample)

    with ui.column().classes("app-shell"):
        ui.html(
            hero_section(
                "Results for {0}".format(sample),
                "Open the full BGC-XPLORER report here, review the generated tables, and jump to the standalone HTML report when you need a shareable view.",
                "Report view",
                )
            )

        ui.html(
            "<section class='panel'>"
            "<div class='section-head'><h2>At a Glance</h2><span class='section-accent'></span></div>"
            "<div class='metrics-grid'>{cards}</div>"
            "</section>".format(
                cards="".join(
                    [
                        stat_card("Total loci", str(metrics["consensus"]), "Approximate candidate loci grouped across callers."),
                        stat_card("Multi-caller loci", str(metrics["multitool"]), "Loci containing predictions from at least two callers."),
                        stat_card("MIBiG comparisons", str(metrics["mibig_hits"]), "Loci with a representative computational MIBiG comparison."),
                        stat_card("ARTS known hits", str(metrics["arts_known_loci"]), "Loci overlapping known-hit ARTS records."),
                    ]
                )
            )
        )

        ui.html(
            "<section class='panel'>"
            "<div class='section-head'><h2>Tool Signals</h2><span class='section-accent'></span></div>"
            "<div class='metrics-grid'>{cards}</div>"
            "</section>".format(
                cards="".join(
                    [
                        stat_card("antiSMASH", str(metrics["antismash"]), "Rule-based BGC region calls and class assignments."),
                        stat_card("GECCO", str(metrics["gecco"]), "Machine-learning BGC candidates from the same Bakta sample."),
                        stat_card("DeepBGC", str(metrics["deepbgc"]), "Domain-driven BGC predictions and activity hints."),
                        stat_card("ARTS", str(metrics["arts"]), "Resistance-linked genomic evidence overlapping candidate loci."),
                        stat_card("dbCAN CGC", str(metrics["dbcan"]), "Carbohydrate gene cluster substrate prediction rows."),
                        stat_card("ARTS DUF hits", str(metrics["arts_duf_loci"]), "Loci overlapping ARTS domains of unknown function."),
                    ]
                )
            )
        )

        with ui.row().classes("w-full info-grid mb-6"):
            ui.html(
                info_card(
                    "Embedded report",
                    "The HTML report is shown directly below so you can browse it without leaving the app.",
                    "▣",
                )
            )
            ui.html(
                info_card(
                    "Standalone links",
                    "Use the extra links for direct HTML access or to download the summary tables that feed the report.",
                    "⇢",
                )
            )

        if report_path.exists():
            report_url = f"/static_results/{sample}/report/{sample}.html"
            with ui.card().classes("panel w-full"):
                with ui.row().classes("items-center justify-between w-full"):
                    ui.html("<div class='section-head' style='margin:0'><h2>Report viewer</h2><span class='section-accent'></span></div>")
                    with ui.row().classes("cta-row"):
                        ui.button("New Upload", on_click=lambda: ui.navigate.to("/")).props(
                            "outline color=primary"
                        )
                        ui.link("Open report in full page", report_url).props("target=_blank").classes("result-link")

                ui.html(
                    f'<div class="iframe-shell"><iframe src="{report_url}" '
                    f'style="width:100%;height:calc(100vh - 120px);border:none;background:#fff;"></iframe></div>'
                )
        else:
            with ui.card().classes("panel w-full"):
                ui.label("Report not found. Run analysis first.").classes(
                    "text-red-400 text-lg font-semibold"
                )

        if summary_dir.exists():
            with ui.card().classes("panel w-full"):
                ui.html(
                    "<div class='section-head'><h2>Summary files</h2><span class='section-accent'></span></div>"
                )
                with ui.row().classes("summary-links"):
                    for f in sorted(summary_dir.iterdir()):
                        if f.is_file():
                            url = f"/static_results/{sample}/summary/{f.name}"
                            ui.link(f.name, url).props("download")


# ─── Helpers ──────────────────────────────────────────────────────────────────

def update_samples_tsv(sample: str, taxon: str):
    tsv_path = WORK_DIR / "config" / "samples.tsv"
    fieldnames = ["sample", "taxon"]
    rows = []
    if tsv_path.exists():
        with open(tsv_path, "r") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                if row["sample"] != sample:
                    rows.append(row)
    rows.append({"sample": sample, "taxon": taxon})
    with open(tsv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def default_taxon() -> str:
    for config_path in (WORK_DIR / "config" / "config.yaml", APP_DIR / "config" / "config.yaml"):
        if not config_path.exists():
            continue
        try:
            with open(config_path, "r") as fh:
                config = yaml.safe_load(fh) or {}
            taxon = config.get("tools", {}).get("arts", {}).get("reference_set", "")
            if taxon:
                return str(taxon)
        except Exception:
            continue
    return "actinobacteria"


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="BGC-XPLORER NiceGUI web interface")
    parser.add_argument("--host", default="0.0.0.0")
    parser.add_argument("--port", type=int, default=8778)
    parser.add_argument("--work-dir", default="/work")
    parser.add_argument("--app-dir", default="/app")
    args = parser.parse_args()

    global WORK_DIR, APP_DIR
    WORK_DIR = Path(args.work_dir)
    APP_DIR = Path(args.app_dir)

    (WORK_DIR / "data" / "bakta").mkdir(parents=True, exist_ok=True)
    (WORK_DIR / "data" / "fasta").mkdir(parents=True, exist_ok=True)
    (WORK_DIR / "results").mkdir(parents=True, exist_ok=True)
    (WORK_DIR / "config").mkdir(parents=True, exist_ok=True)

    config_src = APP_DIR / "config"
    config_dst = WORK_DIR / "config"
    if not (config_dst / "config.yaml").exists() and (config_src / "config.yaml").exists():
        shutil.copytree(str(config_src), str(config_dst), dirs_exist_ok=True)

    app.add_static_files("/static_results", str(WORK_DIR / "results"))

    ui.run(
        host=args.host,
        port=args.port,
        dark=False,
        title="BGC-XPLORER",
        reload=False,
        show=False,
        storage_secret=os.environ.get("NICEGUI_STORAGE_SECRET", "bgc-xplorer-storage-secret"),
    )


if __name__ == "__main__":
    main()
