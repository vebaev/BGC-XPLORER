import json
import re
from pathlib import Path

import pandas as pd

from common import format_consensus_label, load_table_if_exists, write_tsv


HEADER_RE = re.compile(r"<h4>MIBiG 4\.0 matches for ([^<]+)</h4>")
SIMILARITY_RE = re.compile(r'<span class="comparippson-similarity">([0-9.]+)%</span>')
LINK_RE = re.compile(
    r'<div class="comparippson-link"><a href="https://mibig\.secondarymetabolites\.org/go/([^"]+)">([^<]+)</a>:</div>'
)
DESC_RE = re.compile(r'<span class="comparippson-description">([^<]+)</span>')
ALLORF_RE = re.compile(r"allorf_(\d+)_(\d+)")
MIBIG_ID_RE = re.compile(r"(BGC\d{7}\.\d+)")
PERCENT_RE = re.compile(r"([0-9]+(?:\.[0-9]+)?)%")
REGION_RE = re.compile(r"region0*([0-9]+)", re.IGNORECASE)


def extract_regions_from_regions_js(path):
    if not path.exists():
        return []
    text = path.read_text(encoding="utf-8", errors="ignore")
    prefix = "var recordData = "
    if prefix not in text:
        return []
    payload_text = text[text.index(prefix) + len(prefix):]
    for marker in [";\nvar ", ";var ", "\nvar "]:
        if marker in payload_text:
            payload_text = payload_text.split(marker, 1)[0]
            break
    payload_text = payload_text.strip()
    if payload_text.endswith(";"):
        payload_text = payload_text[:-1]
    payload = json.loads(payload_text)
    regions = []
    for record in payload:
        contig_id = record.get("seq_id", "")
        for region in record.get("regions", []):
            region_id = str(region.get("idx", ""))
            locus_tags = {}
            for orf in region.get("orfs", []):
                locus = str(orf.get("locus_tag", "")).strip()
                if locus:
                    locus_tags[locus] = {
                        "start": int(orf.get("start", 0) or 0),
                        "end": int(orf.get("end", 0) or 0),
                    }
            regions.append({
                "contig_id": contig_id,
                "region_id": region_id,
                "start": int(region.get("start", 0) or 0),
                "end": int(region.get("end", 0) or 0),
                "products": [],
                "locus_tags": locus_tags,
            })
    return regions


def extract_regions_from_antismash_json(path):
    if not path.exists():
        return []
    payload = json.loads(path.read_text(encoding="utf-8", errors="ignore"))
    regions = []
    for record in payload.get("records", []):
        contig_id = record.get("id", "")
        features = record.get("features", [])
        locus_tags = {}
        for feature in features:
            locus = str(feature.get("locus_tag", "")).strip()
            if locus:
                locus_tags[locus] = {
                    "start": int(feature.get("start", 0) or 0),
                    "end": int(feature.get("end", 0) or 0),
                }
        for idx, area in enumerate(record.get("areas", []), start=1):
            regions.append({
                "contig_id": contig_id,
                "region_id": str(idx),
                "start": int(area.get("start", 0) or 0),
                "end": int(area.get("end", 0) or 0),
                "products": area.get("products", []) or [],
                "locus_tags": locus_tags,
            })
    return regions


def combine_region_sources(regions_js_regions, json_regions):
    by_key = {}
    for region in json_regions + regions_js_regions:
        key = (region["contig_id"], region["region_id"])
        if key not in by_key:
            by_key[key] = region
            continue
        existing = by_key[key]
        if not existing.get("locus_tags") and region.get("locus_tags"):
            existing["locus_tags"] = region["locus_tags"]
        if not existing.get("products") and region.get("products"):
            existing["products"] = region["products"]
        if not existing.get("start") and region.get("start"):
            existing["start"] = region["start"]
        if not existing.get("end") and region.get("end"):
            existing["end"] = region["end"]
    return list(by_key.values())


def parse_mibig_hits_from_html(path):
    if not path.exists():
        return []
    text = path.read_text(encoding="utf-8", errors="ignore")
    matches = []
    for header in HEADER_RE.finditer(text):
        query = header.group(1).strip()
        block = text[header.end():]
        next_header = HEADER_RE.search(block)
        if next_header:
            block = block[:next_header.start()]
        similarity = SIMILARITY_RE.search(block)
        link = LINK_RE.search(block)
        desc = DESC_RE.search(block)
        if not similarity or not link:
            continue
        matches.append({
            "query_label": query,
            "mibig_similarity": float(similarity.group(1)),
            "best_mibig_id": link.group(1).strip(),
            "best_mibig_label": link.group(2).strip(),
            "best_mibig_product": (desc.group(1).strip() if desc else ""),
            "evidence_source": "comparippson_html",
        })
    return matches


def infer_region_id_from_path(path):
    match = REGION_RE.search(str(path))
    if match:
        return str(int(match.group(1)))
    return ""


def parse_auxiliary_dereplication_files(base_dir):
    hits = []
    if not base_dir.exists():
        return hits
    patterns = [
        "*knownclusterblast*",
        "*KnownClusterBlast*",
        "*clustercompare*",
        "*ClusterCompare*",
    ]
    seen_files = set()
    for pattern in patterns:
        for path in base_dir.rglob(pattern):
            if path.is_dir():
                continue
            seen_files.add(path)
    for path in sorted(seen_files):
        text = path.read_text(encoding="utf-8", errors="ignore")
        region_id = infer_region_id_from_path(path)
        source_name = "knownclusterblast" if "knownclusterblast" in str(path).lower() else "clustercompare"
        mibig_ids = MIBIG_ID_RE.findall(text)
        if not mibig_ids:
            continue
        percents = [float(value) for value in PERCENT_RE.findall(text)]
        best_percent = max(percents) if percents else 0.0
        best_id = mibig_ids[0]
        hits.append({
            "query_label": "region_{0}".format(region_id) if region_id else path.stem,
            "region_id": region_id,
            "mibig_similarity": best_percent,
            "best_mibig_id": best_id,
            "best_mibig_label": best_id,
            "best_mibig_product": "",
            "evidence_source": source_name,
        })
    return hits


def map_query_to_region(query_label, regions):
    if str(query_label).startswith("region_"):
        region_id = str(query_label).split("_", 1)[1]
        for region in regions:
            if str(region["region_id"]) == str(int(region_id)):
                return region
    allorf = ALLORF_RE.fullmatch(query_label)
    if allorf:
        start = int(allorf.group(1))
        end = int(allorf.group(2))
        for region in regions:
            if region["start"] <= start <= region["end"] or region["start"] <= end <= region["end"]:
                return region
        return None
    for region in regions:
        if query_label in region["locus_tags"]:
            return region
    return None


def classify_dereplication(similarity):
    thresholds = snakemake.config["consensus"].get("dereplication_thresholds", {})
    known_like = float(thresholds.get("known_like", 80.0))
    related = float(thresholds.get("related", 50.0))
    divergent = float(thresholds.get("divergent", 20.0))
    if similarity >= known_like:
        return "known-like", 0.1
    if similarity >= related:
        return "related", 0.45
    if similarity >= divergent:
        return "divergent", 0.75
    return "novel_candidate", 1.0


def infer_mibig_class(hit, region, anti_row):
    product_text = str(hit.get("best_mibig_product", "")).strip()
    if ":" in product_text:
        return product_text.split(":", 1)[0].strip()
    anti_type = str(anti_row.get("bgc_type", "")).strip()
    if anti_type:
        return anti_type
    products = region.get("products", [])
    if products:
        return ",".join(str(item).strip() for item in products if str(item).strip())
    return ""


sample = snakemake.wildcards.sample
antismash = load_table_if_exists(snakemake.input.antismash, [
    "sample", "tool", "contig", "start", "end", "strand", "bgc_id", "bgc_type", "product", "score", "confidence", "source_file"
])
consensus = pd.read_csv(snakemake.input.consensus, sep="\t")
index_html = Path(snakemake.input.index_html)
regions_js = Path(snakemake.input.regions_js)
antismash_json = index_html.parent / "{sample}.json".format(sample=sample)

regions = combine_region_sources(
    extract_regions_from_regions_js(regions_js),
    extract_regions_from_antismash_json(antismash_json),
)
hits = parse_mibig_hits_from_html(index_html)
hits.extend(parse_auxiliary_dereplication_files(index_html.parent))

rows = []
for hit in hits:
    region = map_query_to_region(hit["query_label"], regions)
    if not region:
        continue
    region_id = region["region_id"]
    anti_rows = antismash[antismash["bgc_id"].astype(str) == region_id]
    if anti_rows.empty:
        anti_rows = antismash[antismash["bgc_id"].astype(str) == "{0}.0".format(region_id)]
    if anti_rows.empty:
        continue
    anti_row = anti_rows.iloc[0]
    consensus_rows = consensus[
        (consensus["contig"].astype(str) == str(anti_row["contig"]))
        & (pd.to_numeric(consensus["start"], errors="coerce") <= float(anti_row["end"]))
        & (pd.to_numeric(consensus["end"], errors="coerce") >= float(anti_row["start"]))
    ]
    if consensus_rows.empty:
        continue
    dereplication_status, novelty_score = classify_dereplication(hit["mibig_similarity"])
    mibig_class = infer_mibig_class(hit, region, anti_row)
    for _, consensus_row in consensus_rows.iterrows():
        rows.append({
            "sample": sample,
            "consensus_id": consensus_row["consensus_id"],
            "consensus_label": consensus_row.get("consensus_label", format_consensus_label(consensus_row["consensus_id"])),
            "best_mibig_id": hit["best_mibig_id"],
            "best_mibig_product": hit["best_mibig_product"],
            "best_mibig_class": mibig_class,
            "mibig_similarity": hit["mibig_similarity"],
            "dereplication_status": dereplication_status,
            "novelty_score": novelty_score,
            "evidence_source": hit.get("evidence_source", ""),
            "query_label": hit["query_label"],
        })

if rows:
    out = pd.DataFrame(rows).sort_values(
        ["consensus_id", "mibig_similarity", "evidence_source", "best_mibig_id"],
        ascending=[True, False, True, True],
    )
    out = out.drop_duplicates(subset=["consensus_id"], keep="first")
else:
    out = pd.DataFrame(columns=[
        "sample",
        "consensus_id",
        "consensus_label",
        "best_mibig_id",
        "best_mibig_product",
        "best_mibig_class",
        "mibig_similarity",
        "dereplication_status",
        "novelty_score",
        "evidence_source",
        "query_label",
    ])

write_tsv(out, snakemake.output[0])
