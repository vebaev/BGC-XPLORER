#!/usr/bin/env python3
"""Cross-check the recovery counts with BGC-QUAST, the tool Reviewer 1 names.

BGC-QUAST's compare-to-reference mode scores the BGCs of an assembly against
the BGCs of a reference genome. Here the "reference BGCs" are the MIBiG truth
loci, the "assembly" is the reference genome itself (so the QUAST alignment is
the identity), and each prediction set is scored in turn. Every set is written
in the GECCO table layout, which BGC-QUAST reads for any tool.

Usage:
  bgcquast_crosscheck.py --workdir work --variants variants --truth-dir truth \
      --samples S_coelicolor_A3_2,... --out bgcquast
"""
import argparse
import glob
import json
import os
import subprocess
import sys

import pandas as pd

BIN = os.path.dirname(sys.executable)
HERE = os.path.dirname(os.path.abspath(__file__))

from validate_consensus import contig_map, load_arms

ARMS = ["antismash", "gecco", "deepbgc", "side_by_side", "containment:agreed",
        "voting2:union", "voting2:agreed"]
MIBIG_TO_GECCO = {"NRPS": "nrp", "PKS": "polyketide", "ribosomal": "ripp",
                  "terpene": "terpene", "saccharide": "saccharide"}
PROBABILITIES = ["alkaloid_probability", "nrp_probability", "polyketide_probability",
                 "ripp_probability", "saccharide_probability", "terpene_probability"]


def gecco_table(frame, ids, types):
    table = pd.DataFrame({"sequence_id": frame["contig"].values, "cluster_id": ids,
                          "start": frame["start"].values, "end": frame["end"].values,
                          "type": types})
    for column in PROBABILITIES:
        table[column] = 0.0
    return table


def truth_types(classes):
    mapped = sorted({MIBIG_TO_GECCO[c] for c in str(classes).split(";") if c in MIBIG_TO_GECCO})
    return ";".join(mapped) if mapped else "Unknown"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--variants", required=True)
    ap.add_argument("--truth-dir", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rows, per_locus = [], []
    for sample in args.samples.split(","):
        base = os.path.join(args.out, sample)
        genome = os.path.join(args.workdir, "data", "fasta", sample + ".fasta")
        quast = os.path.join(base, "quast")
        if not glob.glob(os.path.join(quast, "contigs_reports", "minimap_output", sample + ".coords")):
            subprocess.run([os.path.join(BIN, "quast.py"), genome, "-r", genome, "-o", quast, "-l", sample,
                            "-t", str(args.threads)], check=True, stdout=subprocess.DEVNULL)
        truth = pd.read_csv(os.path.join(args.truth_dir, sample + ".tsv"), sep="\t")
        reference = os.path.join(base, "reference.clusters.tsv")
        gecco_table(truth, truth["mibig_id"].values, truth["class"].map(truth_types).values).to_csv(
            reference, sep="\t", index=False)

        arms = load_arms(args.workdir, args.variants, sample, contig_map(args.workdir, sample))
        for arm in ARMS:
            if arm not in arms:
                continue
            folder = os.path.join(base, arm.replace(":", "_"))
            os.makedirs(folder, exist_ok=True)
            preds = arms[arm].reset_index(drop=True)
            ids = ["{0}_{1}".format(arm.replace(":", "_"), i + 1) for i in range(len(preds))]
            gecco_table(preds, ids, ["Unknown"] * len(preds)).to_csv(
                os.path.join(folder, sample + ".clusters.tsv"), sep="\t", index=False)
            report_dir = os.path.join(folder, "report")
            statuses_path = os.path.join(folder, "per_reference_status.json")
            subprocess.run([sys.executable, os.path.join(HERE, "bgcquast_percluster.py"), statuses_path,
                            os.path.join(folder, sample + ".clusters.tsv"),
                            "-r", reference, "-R", genome, "-q", quast,
                            "--mode", "compare-to-reference", "-o", report_dir],
                           check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            for mibig_id, status in json.load(open(statuses_path)).items():
                per_locus.append({"sample": sample, "arm": arm, "mibig_id": mibig_id,
                                  "bgcquast_status": status})
            report = pd.read_csv(os.path.join(report_dir, "report.tsv"), sep="\t", index_col=0)
            values = report[sample]
            rows.append({"sample": sample, "arm": arm,
                         "bgcs": int(values["# BGCs"]),
                         "full": int(values["# full ref. BGCs"]),
                         "partial": int(values["# partial ref. BGCs"]),
                         "missed": int(values["# missed ref. BGCs"]),
                         "unmapped_to_ref": int(values["# unmapped BGCs to ref."])})

    table = pd.DataFrame(rows)
    table.to_csv(os.path.join(args.out, "bgcquast_summary.tsv"), sep="\t", index=False)
    pd.DataFrame(per_locus).to_csv(os.path.join(args.out, "bgcquast_per_locus.tsv"), sep="\t", index=False)
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
