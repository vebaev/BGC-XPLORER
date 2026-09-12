#!/usr/bin/env python3
import argparse
import json
import os
import re
import socket
import time
import urllib.error
import urllib.request
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path
from socketserver import ThreadingMixIn
from urllib.parse import unquote

import pandas as pd


DEFAULT_BASE_URL = "https://integrate.api.nvidia.com/v1"
DEFAULT_MODEL = "deepseek-ai/deepseek-v4-pro"
MIN_REQUEST_INTERVAL_SECONDS = 60.0 / 40.0

ANALYSIS_JSON_SCHEMA = {
    "type": "object",
    "properties": {
        "summary": {"type": "string", "description": "Analytical BGC interpretation grounded in the representative genes, tool predictions, ARTS/dbCAN/MIBiG evidence, and relevant scientific literature."},
        "likely_product_or_function": {"type": "string"},
        "biosynthetic_logic": {"type": "string"},
        "key_genes": {"type": "array", "items": {"type": "string"}, "description": "List of strings: locus_tag (bakta_gene / eggnog_description)"},
        "resistance_transport_regulation": {"type": "string"},
        "novelty_assessment": {"type": "string"},
        "confidence": {"type": "string"},
        "caveats": {"type": "string"},
        "recommended_followup": {"type": "array", "items": {"type": "string"}},
    },
    "required": [
        "summary", "likely_product_or_function", "biosynthetic_logic",
        "key_genes", "resistance_transport_regulation", "novelty_assessment",
        "confidence", "caveats", "recommended_followup",
    ],
}

MIME_TYPES = {
    ".html": "text/html; charset=utf-8",
    ".htm": "text/html; charset=utf-8",
    ".svg": "image/svg+xml",
    ".js": "application/javascript; charset=utf-8",
    ".css": "text/css; charset=utf-8",
    ".json": "application/json; charset=utf-8",
    ".png": "image/png",
    ".jpg": "image/jpeg",
    ".jpeg": "image/jpeg",
    ".gif": "image/gif",
    ".ico": "image/x-icon",
    ".tsv": "text/tab-separated-values; charset=utf-8",
    ".txt": "text/plain; charset=utf-8",
}


def guess_mime(path):
    return MIME_TYPES.get(Path(path).suffix.lower(), "application/octet-stream")


SYSTEM_PROMPT = """You are a microbial natural products and biosynthetic gene cluster analyst.
Analyze the supplied BGC candidate using the provided evidence and broader scientific knowledge of natural product biosynthesis. Be useful but cautious.
Do not invent compounds, genes, or database hits. If evidence is weak, say so.
Return ONLY a JSON object with exactly these keys (no markdown, no prose, no code fences):
summary, likely_product_or_function, biosynthetic_logic, key_genes, resistance_transport_regulation,
novelty_assessment, confidence, caveats, recommended_followup.
- summary: a concise but analytically rich narrative. Interpret the representative genes, their annotations and categories, the supporting tools, ARTS/dbCAN/MIBiG evidence, and what kind of BGC this likely is. Ground the interpretation in relevant scientific literature on BGC types when possible.
- key_genes: array of strings in the exact form \"locus_tag (bakta_gene / eggnog_description)\". If bakta_gene is missing, use the eggnog_description only.
- recommended_followup: array of very concise strings.
- All other keys must be strings."""


def clean(value):
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if text.lower() in {"nan", "none", "-"}:
        return ""
    return text


def safe_name(value):
    text = clean(value)
    if not re.match(r"^[A-Za-z0-9_.-]+$", text):
        raise ValueError("Unsafe identifier: {0}".format(text))
    return text


def read_tsv(path):
    if path.exists():
        return pd.read_csv(path, sep="\t")
    return pd.DataFrame()


def row_by_cluster(df, consensus_id):
    if df.empty or "consensus_id" not in df.columns:
        return {}
    rows = df[df["consensus_id"].astype(str) == consensus_id]
    if rows.empty:
        return {}
    return rows.iloc[0].fillna("").to_dict()


def matching_regions(df, contig, start, end):
    if df.empty or not {"contig", "start", "end"}.issubset(df.columns):
        return []
    local = df[df["contig"].astype(str) == str(contig)].copy()
    if local.empty:
        return []
    local["start_num"] = pd.to_numeric(local["start"], errors="coerce")
    local["end_num"] = pd.to_numeric(local["end"], errors="coerce")
    local = local[
        (local["start_num"].fillna(-1).astype(int) <= int(end)) &
        (local["end_num"].fillna(-1).astype(int) >= int(start))
    ]
    return local.fillna("").to_dict(orient="records")


def compact_record(record, keys):
    return {key: clean(record.get(key)) for key in keys if clean(record.get(key))}


def shorten(value, limit=120):
    text = clean(value)
    if len(text) <= limit:
        return text
    return text[:limit].rsplit(" ", 1)[0].rstrip(" ,;:") + "..."


def category_counts(genes):
    counts = {}
    for gene in genes:
        category = clean(gene.get("gene_category")) or "unknown"
        counts[category] = counts.get(category, 0) + 1
    return dict(sorted(counts.items(), key=lambda item: (-item[1], item[0])))


def compact_gene_for_model(gene):
    keys = [
        "locus_tag", "display_label", "gene_category", "bakta_gene", "bakta_product",
        "preferred_name", "eggnog_description", "pfams", "dbcan_hmm",
        "dbcan_subfamily", "dbcan_recommendation", "arts_evidence",
    ]
    compact = {}
    for key in keys:
        value = clean(gene.get(key))
        if not value:
            continue
        compact[key] = shorten(value, 90 if key in {"eggnog_description", "bakta_product"} else 70)
    return compact


def select_model_genes(genes, limit=None):
    """Pick the most informative genes for the model.

    By default no fixed limit is applied; all genes are kept, ordered by informativeness.
    If a limit is supplied, only the top-N most informative genes are returned.

    Ranking priority:
      1. gene_category importance (biosynthetic core > resistance > tailoring > ...)
      2. presence of strong annotation signal (ARTS, dbCAN, Pfam, preferred_name, bakta_gene)
      3. original order (stable, deterministic)
    """
    priority = {
        "biosynthetic core": 0,
        "resistance evidence": 1,
        "tailoring enzyme": 2,
        "transporter": 3,
        "regulator": 4,
        "cazyme / carbohydrate": 5,
        "hypothetical": 8,
    }
    indexed = list(enumerate(genes))

    def sort_key(item):
        index, gene = item
        category = clean(gene.get("gene_category")).lower()
        has_signal = any(clean(gene.get(field)) for field in (
            "arts_evidence", "dbcan_hmm", "dbcan_subfamily", "dbcan_recommendation",
            "pfams", "preferred_name", "bakta_gene"
        ))
        return (priority.get(category, 6), 0 if has_signal else 1, index)

    sorted_genes = sorted(indexed, key=sort_key)
    if limit:
        sorted_genes = sorted_genes[:limit]
    selected = [gene for _, gene in sorted_genes]
    selected.sort(key=lambda gene: int(float(gene.get("gene_start") or 0)))
    return selected


def compact_payload_for_model(payload):
    genes = payload.get("genes", [])
    region = payload.get("region", {})
    compact_region_keys = [
        "contig", "start", "end", "length_bp", "support_tools", "support_count",
        "bgc_types", "products", "biological_interpretation",
        "overlap_relationship", "core_gene_support", "boundary_confidence",
        "confidence_category", "interest_category", "arts_hits", "why_prioritized",
        "why_not_prioritized", "recommended_followup", "best_mibig_id",
        "best_mibig_product", "mibig_similarity", "dereplication_status", "novelty_score",
    ]
    return {
        "sample": payload.get("sample"),
        "consensus_id": payload.get("consensus_id"),
        "region": compact_record(region, compact_region_keys),
        "gene_count": len(genes),
        "gene_category_counts": category_counts(genes),
        "representative_genes": [compact_gene_for_model(gene) for gene in select_model_genes(genes)],
        "tool_predictions": {
            tool: records[:4]
            for tool, records in payload.get("tool_predictions", {}).items()
        },
        "arts_hits": payload.get("arts_hits", [])[:8],
        "dbcan_cgc": payload.get("dbcan_cgc", [])[:4],
        "mibig_dereplication": compact_record(payload.get("mibig_dereplication", {}), [
            "best_mibig_id", "best_mibig_product", "best_mibig_class",
            "mibig_similarity", "dereplication_status", "novelty_score", "evidence_source",
        ]),
    }


def build_fast_preview(payload):
    """Build a deterministic, AI-free preview summary from the payload.

    Used as an instant response when the client passes preview=true (no AI call,
    no cache write). Lets the UI render something useful immediately while a
    deeper AI analysis runs in the background.
    """
    region = payload.get("region", {})
    genes = payload.get("genes", [])
    counts = category_counts(genes)
    top_categories = list(counts.items())[:4]
    tools = [tool for tool, recs in (payload.get("tool_predictions") or {}).items() if recs]
    mibig = payload.get("mibig_dereplication", {}) or {}
    arts = payload.get("arts_hits", []) or []

    summary_bits = []
    bgc_types = clean(region.get("bgc_types"))
    products = clean(region.get("products"))
    interest = clean(region.get("interest_category"))
    contig = clean(region.get("contig"))
    start = clean(region.get("start"))
    end = clean(region.get("end"))
    length = clean(region.get("length_bp"))

    if bgc_types:
        summary_bits.append("Predicted BGC class: {0}.".format(bgc_types))
    if products:
        summary_bits.append("Tool products: {0}.".format(products[:200]))
    if contig and start and end:
        loc = contig
        if length:
            loc = "{0} ({1}-{2}, {3} bp)".format(contig, start, end, length)
        else:
            loc = "{0} ({1}-{2})".format(contig, start, end)
        summary_bits.append("Location: {0}.".format(loc))
    if interest:
        summary_bits.append("Interest category: {0}.".format(interest))
    if top_categories:
        cat_str = ", ".join("{0}={1}".format(k, v) for k, v in top_categories)
        summary_bits.append("Gene categories ({0}): {1}.".format(len(genes), cat_str))
    if tools:
        summary_bits.append("Supporting tools: {0}.".format(", ".join(tools)))
    if arts:
        summary_bits.append("ARTS hits: {0}.".format(len(arts)))
    if mibig.get("best_mibig_id"):
        sim = clean(mibig.get("mibig_similarity"))
        sim_str = " (similarity {0})".format(sim) if sim else ""
        summary_bits.append(
            "Closest MIBiG: {0} - {1}{2}.".format(
                mibig["best_mibig_id"],
                clean(mibig.get("best_mibig_product")) or "unknown product",
                sim_str,
            )
        )

    summary_text = " ".join(summary_bits) or "No structured data available for this cluster."

    key_genes = []
    for gene in select_model_genes(genes, limit=5):
        label = clean(gene.get("display_label")) or clean(gene.get("preferred_name")) or clean(gene.get("bakta_gene"))
        category = clean(gene.get("gene_category"))
        if label and category:
            key_genes.append("{0} ({1})".format(label, category))
        elif label:
            key_genes.append(label)

    caveats_bits = []
    relevance = clean(region.get("why_prioritized"))
    if relevance:
        caveats_bits.append("Why prioritized: {0}".format(relevance[:160]))
    if not tools:
        caveats_bits.append("No tool predictions inside the region; treat with caution.")

    return {
        "summary": summary_text,
        "likely_product_or_function": products or "Unclassified",
        "biosynthetic_logic": bgc_types or "Not predicted",
        "key_genes": key_genes,
        "resistance_transport_regulation": "",
        "novelty_assessment": "Run a full AI analysis to assess novelty." if not mibig.get("best_mibig_id") else "Closest MIBiG: {0}".format(mibig["best_mibig_id"]),
        "confidence": "Preview (deterministic, no model inference).",
        "caveats": " ".join(caveats_bits) or "Preview summary derived directly from the structured tables; no AI call was made.",
        "recommended_followup": [
            "Trigger a full AI analysis (POST /analyze_cluster without preview) for narrative interpretation.",
            "Inspect the cluster map for biosynthetic logic details.",
        ],
        "_preview": True,
    }


def build_cluster_payload(results_dir, sample, consensus_id):
    summary_dir = results_dir / sample / "summary"
    prioritized = read_tsv(summary_dir / "prioritized_regions.tsv")
    consensus = read_tsv(summary_dir / "consensus_bgcs.tsv")
    cluster_genes = read_tsv(summary_dir / "cluster_genes.tsv")
    antismash = read_tsv(summary_dir / "antismash.bgc.tsv")
    gecco = read_tsv(summary_dir / "gecco.bgc.tsv")
    deepbgc = read_tsv(summary_dir / "deepbgc.bgc.tsv")
    arts = read_tsv(summary_dir / "arts.hits.tsv")
    dbcan_cgc = read_tsv(summary_dir / "dbcan.cgc.tsv")
    mibig = read_tsv(summary_dir / "mibig_dereplication.tsv")

    priority_row = row_by_cluster(prioritized, consensus_id)
    consensus_row = row_by_cluster(consensus, consensus_id)
    cluster_row = priority_row or consensus_row
    if not cluster_row:
        raise ValueError("Unknown consensus_id: {0}".format(consensus_id))

    contig = clean(cluster_row.get("contig")) or clean(cluster_row.get("contig_id"))
    start = int(float(cluster_row.get("start")))
    end = int(float(cluster_row.get("end")))

    genes = cluster_genes[cluster_genes["consensus_id"].astype(str) == consensus_id].copy()
    gene_keys = [
        "locus_tag", "gene_start", "gene_end", "strand", "bakta_gene", "bakta_product",
        "preferred_name", "eggnog_description", "COG_category", "ec", "kegg_ko", "pfams",
        "dbcan_hmm", "dbcan_subfamily", "dbcan_diamond", "dbcan_recommendation",
        "dbcan_substrate", "arts_evidence", "gene_category", "display_label",
    ]
    gene_records = [
        compact_record(row, gene_keys)
        for row in genes.fillna("").to_dict(orient="records")
    ]

    tool_keys = ["tool", "bgc_id", "bgc_type", "product", "score", "confidence", "start", "end"]
    arts_keys = ["feature", "score", "evidence", "start", "end"]
    dbcan_keys = [
        "cgc_id", "cluster_start", "cluster_end", "genes", "cazyme_genes",
        "signature_genes", "predicted_substrate", "prediction_source",
        "pul_id", "pul_substrate", "dbcan_sub_substrate",
    ]

    return {
        "sample": sample,
        "consensus_id": consensus_id,
        "region": compact_record(cluster_row, [
            "contig", "contig_id", "start", "end", "length_bp", "support_tools",
            "supporting_tools", "support_count", "bgc_types", "products",
            "biological_interpretation", "overlap_relationship", "core_gene_support",
            "boundary_confidence", "priority_score", "priority_class",
            "confidence_category", "interest_category", "arts_hits", "why_prioritized",
            "why_not_prioritized", "recommended_followup", "best_mibig_id",
            "best_mibig_product", "mibig_similarity", "dereplication_status", "novelty_score",
        ]),
        "tool_predictions": {
            "antismash": [compact_record(row, tool_keys) for row in matching_regions(antismash, contig, start, end)],
            "gecco": [compact_record(row, tool_keys) for row in matching_regions(gecco, contig, start, end)],
            "deepbgc": [compact_record(row, tool_keys) for row in matching_regions(deepbgc, contig, start, end)],
        },
        "arts_hits": [compact_record(row, arts_keys) for row in matching_regions(arts, contig, start, end)],
        "dbcan_cgc": [
            compact_record(row, dbcan_keys)
            for row in matching_regions(
                dbcan_cgc.rename(columns={"cluster_start": "start", "cluster_end": "end"}),
                contig,
                start,
                end,
            )
        ],
        "mibig_dereplication": row_by_cluster(mibig, consensus_id),
        "genes": gene_records,
    }


def _repair_json(text):
    """Attempt to repair truncated JSON by closing open strings and braces."""
    if not text or text.lstrip()[0] != '{':
        return None
    s = text.strip()
    depth = 0
    in_string = False
    escape = False
    for ch in s:
        if escape:
            escape = False
            continue
        if ch == '\\':
            escape = True
            continue
        if ch == '"':
            in_string = not in_string
            continue
        if in_string:
            continue
        if ch == '{':
            depth += 1
        elif ch == '}':
            depth -= 1
    repair = s
    if in_string:
        repair += '"'
    repair = repair.rstrip()
    repair = repair.rstrip(',')
    repair = repair.rstrip()
    if depth > 0:
        repair += '}' * depth
    return repair


def _extract_fields(text):
    """Best-effort field extraction from partial/truncated JSON."""
    fields = {}
    for key in ("summary", "likely_product_or_function", "biosynthetic_logic",
                "resistance_transport_regulation", "novelty_assessment",
                "confidence", "caveats"):
        m = re.search(r'"' + key + r'"\s*:\s*"((?:[^"\\]|\\.)*)', text, flags=re.S)
        if m:
            fields[key] = m.group(1).encode().decode("unicode_escape", errors="ignore")
    for key in ("key_genes", "recommended_followup"):
        m = re.search(r'"' + key + r'"\s*:\s*\[(.*?)\]', text, flags=re.S)
        if m:
            items = re.findall(r'"((?:[^"\\]|\\.)*)"', m.group(1))
            if items:
                fields[key] = items
    return fields


def extract_json(text):
    text = (text or "").strip()
    fenced = re.match(r"^```(?:json)?\s*(.*?)\s*```$", text, flags=re.S | re.I)
    if fenced:
        text = fenced.group(1).strip()
    # 1. direct parse
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        pass
    # 2. greedy {.*}
    match = re.search(r"\{.*\}", text, flags=re.S)
    if match:
        try:
            return json.loads(match.group(0))
        except json.JSONDecodeError:
            pass
    # 3. repair truncated JSON (close open strings/braces)
    repaired = _repair_json(text)
    if repaired and repaired != text:
        try:
            return json.loads(repaired)
        except json.JSONDecodeError:
            pass
    # 4. best-effort field extraction from partial JSON
    fields = _extract_fields(text)
    if fields:
        fields.setdefault("summary", "")
        fields.setdefault("key_genes", [])
        fields.setdefault("recommended_followup", [])
        return fields
    # 5. fallback
    return {
        "summary": limit_words(text or "The model returned an empty response.", 120),
        "likely_product_or_function": "",
        "biosynthetic_logic": "",
        "key_genes": [],
        "resistance_transport_regulation": "",
        "novelty_assessment": "",
        "confidence": "Not assessed from structured JSON because the model returned plain text.",
        "caveats": "Model response was not valid JSON; displayed as plain-text summary.",
        "recommended_followup": ["Repeat analysis or inspect the gene map manually."],
    }


def limit_words(text, limit):
    words = str(text).split()
    if len(words) <= limit:
        return str(text).strip()
    return " ".join(words[:limit]).strip() + "..."


def normalize_analysis(analysis):
    if not isinstance(analysis, dict):
        analysis = {"summary": str(analysis)}
    normalized = {
        "summary": clean(analysis.get("summary", "")),
        "likely_product_or_function": clean(analysis.get("likely_product_or_function", "")),
        "biosynthetic_logic": clean(analysis.get("biosynthetic_logic", "")),
        "key_genes": normalize_list(analysis.get("key_genes")),
        "resistance_transport_regulation": clean(analysis.get("resistance_transport_regulation", "")),
        "novelty_assessment": clean(analysis.get("novelty_assessment", "")),
        "confidence": clean(analysis.get("confidence", "")),
        "caveats": clean(analysis.get("caveats", "")),
        "recommended_followup": normalize_list(analysis.get("recommended_followup")),
    }
    if not normalized["summary"]:
        normalized["summary"] = "No concise summary returned by the model."
    return normalized


def normalize_list(value):
    if not value:
        return []
    if isinstance(value, str):
        value = [value]
    if not isinstance(value, list):
        value = [str(value)]
    return [clean(item) for item in value if clean(item)]


class RateLimiter:
    def __init__(self, min_interval):
        self.min_interval = min_interval
        self.last_request = 0.0

    def wait(self):
        elapsed = time.time() - self.last_request
        remaining = self.min_interval - elapsed
        if remaining > 0:
            time.sleep(remaining)
        self.last_request = time.time()


class ThreadingHTTPServer(ThreadingMixIn, HTTPServer):
    daemon_threads = True


class AIClusterServer(BaseHTTPRequestHandler):
    server_version = "AktinoAIClusterServer/0.1"

    def _headers(self, status=200):
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
        self.end_headers()

    def _json(self, data, status=200):
        self._headers(status)
        self.wfile.write(json.dumps(data, ensure_ascii=False, indent=2).encode("utf-8"))

    def do_OPTIONS(self):
        self._headers(204)

    def do_GET(self):
        if self.path == "/health":
            self._json({
                "ok": True,
                "model": self.server.model,
                "api_key_loaded": bool(self.server.api_key),
                "timeout_seconds": self.server.timeout,
                "max_tokens": self.server.max_tokens,
                "reasoning_effort": self.server.reasoning_effort,
                "guided_json": self.server.use_guided_json,
                "report_dir": str(getattr(self.server, "report_dir", "") or ""),
            })
            return
        self._serve_static(self.path)

    def _serve_static(self, raw_path):
        report_dir = getattr(self.server, "report_dir", None)
        if not report_dir:
            self._json({"error": "No report directory configured. Use --report-dir."}, status=404)
            return
        report_root = report_dir.resolve()
        path = unquote(raw_path.split("?", 1)[0].split("#", 1)[0].lstrip("/"))
        if not path or path.endswith("/"):
            htmls = sorted(report_root.glob("*.html"))
            if not htmls:
                self._json({"error": "No HTML report found in {0}.".format(report_root)}, status=404)
                return
            file_path = htmls[0]
        else:
            file_path = (report_root / path).resolve()
        try:
            file_path.relative_to(report_root)
        except ValueError:
            self._json({"error": "Forbidden."}, status=403)
            return
        if not file_path.is_file():
            self._json({"error": "Not found: {0}".format(path)}, status=404)
            return
        try:
            with open(file_path, "rb") as handle:
                data = handle.read()
        except OSError as exc:
            self._json({"error": str(exc)}, status=500)
            return
        self.send_response(200)
        self.send_header("Content-Type", guess_mime(file_path))
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Access-Control-Allow-Origin", "*")
        self.end_headers()
        self.wfile.write(data)

    def do_POST(self):
        if self.path != "/analyze_cluster":
            self._json({"error": "Not found."}, status=404)
            return
        try:
            length = int(self.headers.get("Content-Length", "0"))
            body = json.loads(self.rfile.read(length).decode("utf-8") or "{}")
            sample = safe_name(body.get("sample"))
            consensus_id = safe_name(body.get("consensus_id"))
            force = bool(body.get("force", False))
            preview = bool(body.get("preview", False))
            cache_path = self.server.results_dir / sample / "ai_cluster_reports" / "{0}.json".format(consensus_id)

            if cache_path.exists() and not force:
                with open(cache_path, "r", encoding="utf-8") as handle:
                    cached = json.load(handle)
                cached["cached"] = True
                print("[cache hit] {0}/{1}".format(sample, consensus_id), flush=True)
                self._json(cached)
                return

            payload = build_cluster_payload(self.server.results_dir, sample, consensus_id)

            if preview:
                preview_payload = build_fast_preview(payload)
                self._json({
                    "sample": sample,
                    "consensus_id": consensus_id,
                    "model": self.server.model,
                    "cached": False,
                    "preview": True,
                    "created_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                    "analysis": preview_payload,
                })
                return

            if not self.server.api_key:
                self._json({
                    "error": "NVIDIA_API_KEY is not set. Start the server with your NVIDIA API key in the environment.",
                    "sample": sample,
                    "consensus_id": consensus_id,
                }, status=503)
                return

            payload = build_cluster_payload(self.server.results_dir, sample, consensus_id)
            model_payload = compact_payload_for_model(payload)
            prompt = (
                "Analyze the biosynthetic gene cluster candidate described below. "
                "In the 'summary' field provide an analytical interpretation of the genes, "
                "their annotations and categories, the supporting tools, and any ARTS/dbCAN/MIBiG "
                "evidence. Ground the likely BGC type in relevant scientific literature on natural "
                "product biosynthesis. In the 'key_genes' field list the most important genes as: "
                "locus_tag (bakta_gene / eggnog_description). "
                "Reply with ONLY a JSON object containing exactly these 9 keys: summary, "
                "likely_product_or_function, biosynthetic_logic, key_genes, "
                "resistance_transport_regulation, novelty_assessment, confidence, caveats, "
                "recommended_followup. No prose or markdown outside the JSON.\n\n"
                + json.dumps(model_payload, ensure_ascii=False, separators=(",", ":"))
            )
            request_body = {
                "model": self.server.model,
                "messages": [
                    {"role": "system", "content": SYSTEM_PROMPT},
                    {"role": "user", "content": prompt},
                ],
                "temperature": self.server.temperature,
                "top_p": self.server.top_p,
                "max_tokens": self.server.max_tokens,
                "reasoning_effort": self.server.reasoning_effort,
                "seed": self.server.seed,
            }
            if self.server.use_guided_json:
                request_body["nvext"] = {"guided_json": ANALYSIS_JSON_SCHEMA}
            diagnostics = {
                "sample": sample,
                "consensus_id": consensus_id,
                "model": self.server.model,
                "timeout_seconds": self.server.timeout,
                "payload_chars": len(json.dumps(model_payload, ensure_ascii=False, separators=(",", ":"))),
                "gene_count": len(payload.get("genes", [])),
                "representative_gene_count": len(model_payload.get("representative_genes", [])),
                "guided_json": bool(self.server.use_guided_json),
            }
            self.server.rate_limiter.wait()
            started_at = time.time()
            try:
                raw_content = self.call_nvidia(request_body)
            except RuntimeError as exc:
                if self.server.use_guided_json and "guided_json" in str(exc).lower():
                    try:
                        request_body.pop("nvext", None)
                        diagnostics["guided_json_retry"] = True
                        raw_content = self.call_nvidia(request_body)
                    except Exception as retry_exc:
                        diagnostics.update({
                            "error": str(retry_exc),
                            "error_type": retry_exc.__class__.__name__,
                            "stage": "nvidia_retry_without_guided_json",
                            "elapsed_seconds": round(time.time() - started_at, 2),
                        })
                        self._json(diagnostics, status=504 if isinstance(retry_exc, (TimeoutError, socket.timeout)) else 502)
                        return
                else:
                    diagnostics.update({
                        "error": str(exc),
                        "error_type": exc.__class__.__name__,
                        "stage": "nvidia_api",
                        "elapsed_seconds": round(time.time() - started_at, 2),
                    })
                    self._json(diagnostics, status=502)
                    return
            except Exception as exc:
                diagnostics.update({
                    "error": str(exc),
                    "error_type": exc.__class__.__name__,
                    "stage": "nvidia_request",
                    "elapsed_seconds": round(time.time() - started_at, 2),
                })
                self._json(diagnostics, status=504 if isinstance(exc, (TimeoutError, socket.timeout)) else 502)
                return
            analysis = normalize_analysis(raw_content)
            response = {
                "sample": sample,
                "consensus_id": consensus_id,
                "model": self.server.model,
                "cached": False,
                "diagnostics": diagnostics,
                "created_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "analysis": analysis,
            }
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            with open(cache_path, "w", encoding="utf-8") as handle:
                json.dump(response, handle, ensure_ascii=False, indent=2)
            self._json(response)
        except Exception as exc:
            self._json({"error": str(exc)}, status=500)

    def call_nvidia(self, request_body):
        url = self.server.base_url.rstrip("/") + "/chat/completions"
        request = urllib.request.Request(
            url,
            data=json.dumps(request_body).encode("utf-8"),
            headers={
                "Authorization": "Bearer {0}".format(self.server.api_key),
                "Content-Type": "application/json",
            },
            method="POST",
        )
        try:
            with urllib.request.urlopen(request, timeout=self.server.timeout) as response:
                data = json.loads(response.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            detail = exc.read().decode("utf-8", errors="ignore")
            raise RuntimeError("NVIDIA API error {0}: {1}".format(exc.code, detail))
        content = data["choices"][0]["message"]["content"]
        return extract_json(content)


def main():
    parser = argparse.ArgumentParser(description="Local on-demand AI analysis service for Aktino BGC reports.")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8787)
    parser.add_argument("--results-dir", default="results")
    parser.add_argument("--report-dir", default=None,
                        help="Directory with the HTML report to serve (auto-detected from results-dir if omitted)")
    parser.add_argument("--model", default=os.environ.get("NVIDIA_MODEL", DEFAULT_MODEL))
    parser.add_argument("--base-url", default=os.environ.get("NVIDIA_API_BASE", DEFAULT_BASE_URL))
    parser.add_argument("--temperature", type=float, default=float(os.environ.get("AI_TEMPERATURE", "0.1")))
    parser.add_argument("--top-p", type=float, default=float(os.environ.get("AI_TOP_P", "0.95")))
    parser.add_argument("--max-tokens", type=int, default=int(os.environ.get("AI_MAX_TOKENS", "1800")))
    parser.add_argument("--reasoning-effort", default=os.environ.get("AI_REASONING_EFFORT", "none"),
                        choices=["none", "high", "max"],
                        help="DeepSeek reasoning mode: none (fast, strict JSON) or high/max (deeper analysis)")
    parser.add_argument("--seed", type=int, default=int(os.environ.get("AI_SEED", "42")))
    parser.add_argument("--no-guided-json", action="store_true",
                        help="Disable NVIDIA guided_json constrained decoding (fallback to prompt-based JSON)")
    parser.add_argument("--timeout", type=int, default=int(os.environ.get("AI_TIMEOUT", "240")))
    args = parser.parse_args()

    server = ThreadingHTTPServer((args.host, args.port), AIClusterServer)
    server.results_dir = Path(args.results_dir).resolve()
    report_dir_arg = args.report_dir
    if not report_dir_arg:
        candidates = sorted(server.results_dir.glob("*/report"))
        report_dir_arg = str(candidates[0]) if candidates else None
    server.report_dir = Path(report_dir_arg).resolve() if report_dir_arg else None
    server.api_key = os.environ.get("NVIDIA_API_KEY", "")
    server.model = args.model
    server.base_url = args.base_url
    server.temperature = args.temperature
    server.top_p = args.top_p
    server.max_tokens = args.max_tokens
    server.reasoning_effort = args.reasoning_effort
    server.seed = args.seed
    server.use_guided_json = not args.no_guided_json
    server.timeout = args.timeout
    server.rate_limiter = RateLimiter(MIN_REQUEST_INTERVAL_SECONDS)

    print("AI cluster server listening on http://{0}:{1}".format(args.host, args.port))
    print("Model: {0}".format(server.model))
    print("NVIDIA_API_KEY loaded: {0}".format("yes" if server.api_key else "no"))
    print("Parameters: temperature={t}, top_p={p}, max_tokens={m}, reasoning_effort={r}, seed={s}, guided_json={g}".format(
        t=server.temperature, p=server.top_p, m=server.max_tokens,
        r=server.reasoning_effort, s=server.seed, g=server.use_guided_json))
    if server.report_dir:
        print("Serving report from: {0}".format(server.report_dir))
        htmls = sorted(server.report_dir.glob("*.html"))
        if htmls:
            print("Open: http://{0}:{1}/{2}".format(args.host, args.port, htmls[0].name))
        else:
            print("Open: http://{0}:{1}/".format(args.host, args.port))
    else:
        print("No report directory found; only /health and /analyze_cluster are available.")
    server.serve_forever()


if __name__ == "__main__":
    main()
