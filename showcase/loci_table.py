#!/usr/bin/env python3
"""Supplementary table of every locus in a finished BGC-XPLORER run.

One row per locus, built only from the workflow's own outputs, so the table can
be regenerated after any re-run. Loci are sorted into the three groups a reader
has to treat differently: multi-caller loci that include antiSMASH, multi-caller
loci supported only by GECCO and DeepBGC (two Pfam-based callers, so weaker
evidence than their count suggests), and single-caller loci.

Usage:
  loci_table.py --results results --sample Soil_1 --out showcase/Soil_1_loci.xlsx
  (a .tsv path writes a tab-separated table instead)
"""
import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
from bgc_classes import activity_labels, harmonised_classes  # noqa: E402

CALLER_NAMES = {"antismash": "antiSMASH", "gecco": "GECCO", "deepbgc": "DeepBGC"}


def callers_of(value):
    return [tool.strip().lower() for tool in str(value).split(",") if tool.strip() and tool.strip().lower() != "nan"]


def category(tools):
    if len(tools) >= 2:
        return "Multi-caller, with antiSMASH" if "antismash" in tools else "Multi-caller, GECCO + DeepBGC only"
    return "Single caller"


def text(value):
    return "" if pd.isna(value) else str(value).strip()


def build_table(results, sample):
    summary = os.path.join(results, sample, "summary")
    loci = pd.read_csv(os.path.join(summary, "region_evidence.tsv"), sep="\t")
    genes = pd.read_csv(os.path.join(results, sample, "bakta", sample + ".features.tsv"), sep="\t")
    genes = genes[genes["type"].astype(str) == "cds"]
    support = pd.read_csv(os.path.join(summary, "gene_caller_support.tsv"), sep="\t")
    callers_per_gene = dict(zip(support["locus_tag"].astype(str), support["callers"]))

    rows = []
    for _, locus in loci.iterrows():
        tools = callers_of(locus.get("support_tools"))
        start, end = int(locus["start"]), int(locus["end"])
        inside = genes[(genes["contig"].astype(str) == str(locus["contig"]))
                       & (genes["start"] >= start) & (genes["end"] <= end)]
        tags = inside["locus_tag"].astype(str)
        multi = sum(1 for tag in tags if callers_per_gene.get(tag, 0) >= 2)
        score = text(locus.get("match_score"))
        rows.append({
            "Locus": text(locus.get("consensus_label")) or text(locus.get("consensus_id")),
            "Contig": text(locus.get("contig")),
            "Start": start,
            "End": end,
            "Agreed length (bp)": end - start + 1,
            "Gene-map window": "{0}-{1}".format(int(locus["union_start"]), int(locus["union_end"])),
            "Callers": ", ".join(CALLER_NAMES.get(tool, tool) for tool in tools),
            "Number of callers": len(tools),
            "Category": category(tools),
            "Genes in agreed span": len(tags),
            "Genes covered by >=2 callers": multi,
            "Biosynthetic core gene": "yes" if text(locus.get("core_gene_evidence")) else "no",
            "MIBiG representative": text(locus.get("best_mibig_id")),
            "MIBiG compound": text(locus.get("best_mibig_product")),
            "MIBiG score": "{0} ({1})".format(score, text(locus.get("score_metric"))) if score else "",
            "MIBiG core genes shared": text(locus.get("core_gene_hits")),
            "ARTS known hits": int(locus.get("arts_known_hits") or 0),
            "ARTS DUF hits": int(locus.get("arts_duf_hits") or 0),
            "Classes (MIBiG)": ", ".join(harmonised_classes(locus.get("bgc_types"))),
            "Predicted activity (DeepBGC)": ", ".join(activity_labels(locus.get("products"))),
        })
    order = {"Multi-caller, with antiSMASH": 0, "Multi-caller, GECCO + DeepBGC only": 1, "Single caller": 2}
    table = pd.DataFrame(rows)
    table["_group"] = table["Category"].map(order)
    return table.sort_values(["_group", "Contig", "Start"]).drop(columns="_group").reset_index(drop=True)


def summary_lines(table):
    single = table[table["Category"] == "Single caller"]
    lines = ["{0} loci".format(len(table))]
    for name, group in table.groupby("Category", sort=False):
        lines.append("  {0}: {1} ({2} with a biosynthetic core gene)".format(
            name, len(group), int((group["Biosynthetic core gene"] == "yes").sum())))
    for caller, group in single.groupby("Callers"):
        lines.append("    single {0}: {1}".format(caller, len(group)))
    return lines


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    table = build_table(args.results, args.sample)
    if args.out.endswith(".xlsx"):
        table.to_excel(args.out, index=False, sheet_name=args.sample[:31])
    else:
        table.to_csv(args.out, sep="\t", index=False)
    print("\n".join(summary_lines(table)))


if __name__ == "__main__":
    main()
