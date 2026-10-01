#!/usr/bin/env python3
"""Summarise evidence-mode results over many loci and draw a sample for manual audit.

For every stored summary the evidence the model was sent is rebuilt from the
run (it is deterministic; the stored request fingerprint confirms it), the
statements are verified again with scripts/ai_evidence.py, and the totals per
locus and overall are written. A random sample of verified statements, and all
unverified ones, are written with their cited values for a person to judge:
the sample estimates how often "verified" is nonetheless wrong.

Usage:
  evaluate_evidence_mode.py --results RESULTS_DIR --sample Soil_1 \
      --summaries DIR [DIR ...] --out OUT_DIR [--audit 100]
"""
import argparse
import csv
import glob
import json
import os
import random
import sys
from pathlib import Path

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "scripts"))

import ai_cluster_server as server  # noqa: E402
import ai_evidence as evidence  # noqa: E402

SEED = 20261001


def cited_values(statement, payload):
    out = []
    for path in statement["evidence"]:
        try:
            _, value = evidence._resolve(payload, path)
            out.append("{0} = {1}".format(path, evidence._flat_original(value)[:200]))
        except KeyError as error:
            out.append("{0} -> {1}".format(path, error))
    return " | ".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--summaries", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--audit", type=int, default=100)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    per_locus, verified_rows, unverified_rows = [], [], []
    for folder in args.summaries:
        for path in sorted(glob.glob(os.path.join(folder, "*.json"))):
            if path.endswith(".payload.json"):
                continue
            record = json.load(open(path))
            if record.get("error") or record.get("mode") != "evidence":
                per_locus.append({"consensus_id": os.path.basename(path)[:-5], "error": record.get("error", "not evidence mode")})
                continue
            payload = server.compact_payload_for_model(
                server.build_cluster_payload(Path(args.results), args.sample, record["consensus_id"]))
            check = evidence.verify(record["analysis"], payload)
            cover = check["coverage"]
            per_locus.append({
                "consensus_id": record["consensus_id"], "set": os.path.basename(folder.rstrip("/")),
                "model": record["model"], "created_at": record["created_at"],
                "statements": check["statements"], "verified": check["verified"],
                **{"coverage_" + key: "{0}/{1}".format(value["covered"], value["expected"]) for key, value in cover.items()},
            })
            for result in check["results"]:
                statement = record["analysis"][result["section"]][result["index"]]
                row = {"consensus_id": record["consensus_id"], "section": result["section"],
                       "text": statement["text"], "cited": cited_values(statement, payload),
                       "reasons": "; ".join(result["reasons"]), "judgement": "", "note": ""}
                (verified_rows if result["verified"] else unverified_rows).append(row)

    fields = sorted({k for row in per_locus for k in row}, key=lambda k: (k != "consensus_id", k))
    with open(os.path.join(args.out, "per_locus.tsv"), "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(per_locus)
    rng = random.Random(SEED)
    sample = rng.sample(verified_rows, min(args.audit, len(verified_rows)))
    for name, rows in (("audit_verified_sample.tsv", sample), ("audit_unverified_all.tsv", unverified_rows)):
        with open(os.path.join(args.out, name), "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else ["consensus_id"], delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    ok = [r for r in per_locus if "statements" in r]
    total = sum(r["statements"] for r in ok)
    good = sum(r["verified"] for r in ok)
    print("loci: {0} ({1} with errors); statements: {2}; verified: {3} ({4:.1%})".format(
        len(ok), len(per_locus) - len(ok), total, good, good / total if total else 0))
    for key in ("callers", "core_genes", "mibig", "arts"):
        pairs = [r["coverage_" + key].split("/") for r in ok if "coverage_" + key in r]
        if pairs:
            print("  coverage {0}: {1}/{2}".format(key, sum(int(a) for a, _ in pairs), sum(int(b) for _, b in pairs)))


if __name__ == "__main__":
    main()
