import json
import re
from pathlib import Path

import pandas as pd

from common import ARTS_COLUMNS, coerce_frame_columns, empty_arts_frame, write_tsv


LOC_RE = re.compile(r"loc\|(\d+)[ _](\d+)[ _](-?1)")


def normalize_contig_name(value):
    text = str(value).strip()
    if not text:
        return ""
    if text.startswith("scaffold_"):
        return text.replace("scaffold_", "contig_", 1)
    return text


def parse_loc_string(text, fallback_contig=""):
    value = str(text).strip()
    if not value:
        return None
    match = LOC_RE.search(value)
    if not match:
        return None
    start, end, strand = match.groups()
    contig = fallback_contig
    parts = value.split("|")
    for idx, part in enumerate(parts):
        if part == "source" and idx + 1 < len(parts):
            contig = parts[idx + 1]
            break
    contig = normalize_contig_name(contig)
    return {
        "contig": contig,
        "start": int(start),
        "end": int(end),
        "strand": strand,
    }


def add_row(rows, sample, contig, start, end, feature, score="", evidence="", source_file=""):
    if not contig or start is None or end is None:
        return
    rows.append({
        "sample": sample,
        "tool": "arts",
        "contig": contig,
        "start": int(start),
        "end": int(end),
        "feature": str(feature).strip(),
        "score": str(score).strip(),
        "evidence": str(evidence).strip(),
        "source_file": str(source_file),
    })


def parse_knownhits_tsv(path, sample, rows):
    if not path.exists():
        return
    df = pd.read_csv(path, sep="\t")
    for _, row in df.iterrows():
        info = parse_loc_string(row.get("Sequence description", ""))
        if not info:
            continue
        feature = row.get("Description") or row.get("#Model") or "known_hit"
        add_row(
            rows,
            sample,
            info["contig"],
            info["start"],
            info["end"],
            feature,
            score=row.get("bitscore", ""),
            evidence="known_hit",
            source_file=path,
        )


def parse_json_seqs(path, sample, rows, evidence_label):
    if not path.exists():
        return
    data = json.loads(path.read_text())
    seqs = data.get("seqs", {})
    for seq_id, values in seqs.items():
        if len(values) < 6:
            continue
        contig = normalize_contig_name(values[2])
        start = values[3]
        end = values[4]
        score = values[7] if len(values) > 7 else values[6] if len(values) > 6 else ""
        feature = values[1] if len(values) > 1 else seq_id
        add_row(
            rows,
            sample,
            contig,
            start,
            end,
            feature,
            score=score,
            evidence=evidence_label,
            source_file=path,
        )


def parse_hits_listed_table(path, sample, rows, feature_column, hits_column, evidence_label, score_column=None):
    if not path.exists():
        return
    df = pd.read_csv(path, sep="\t")
    for _, row in df.iterrows():
        feature = row.get(feature_column, "")
        score = row.get(score_column, "") if score_column else ""
        hits_text = str(row.get(hits_column, "")).strip()
        for item in hits_text.strip("[]").split(";"):
            info = parse_loc_string(item)
            if not info:
                continue
            add_row(
                rows,
                sample,
                info["contig"],
                info["start"],
                info["end"],
                feature,
                score=score,
                evidence=evidence_label,
                source_file=path,
            )


sample = snakemake.wildcards.sample
outdir = Path(snakemake.input.done).parent
rows = []

parse_knownhits_tsv(outdir / "tables" / "knownhits.tsv", sample, rows)
parse_json_seqs(outdir / "tables" / "dufhits.json", sample, rows, "duf_hit")
parse_hits_listed_table(outdir / "tables" / "coretable.tsv", sample, rows, "#Core_gene", "[Hits_listed]", "core_gene")
parse_hits_listed_table(outdir / "tables" / "duptable.tsv", sample, rows, "#Core_gene", "[Hits_listed]", "duplication_signal", "Count")

if rows:
    out = pd.DataFrame(rows).drop_duplicates(
        subset=["contig", "start", "end", "feature", "evidence", "source_file"]
    )
    out = coerce_frame_columns(out, ARTS_COLUMNS)
else:
    out = empty_arts_frame()

write_tsv(out, snakemake.output[0])
