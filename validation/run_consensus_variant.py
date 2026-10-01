#!/usr/bin/env python3
"""Run the pipeline's own scripts/build_consensus.py under another grouping rule.

The caller outputs of a finished run are reused; only the grouping changes.
The script is executed unmodified, with a stand-in for the `snakemake` object
Snakemake would inject, so every variant is produced by the released code
rather than by a re-implementation of it.

Usage:
  run_consensus_variant.py --workdir work --sample S_coelicolor_A3_2 \
      --mode containment --out variants/S_coelicolor_A3_2/containment
  run_consensus_variant.py ... --mode voting --min-callers 3 --out .../voting3
"""
import argparse
import os
import runpy
import sys
from types import SimpleNamespace

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True, help="directory holding results/ and data/")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--mode", choices=["voting", "containment"], required=True)
    ap.add_argument("--min-callers", type=int, default=2)
    ap.add_argument("--min-containment", type=float, default=0.80)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    summary = os.path.join(args.workdir, "results", args.sample, "summary")
    os.makedirs(args.out, exist_ok=True)
    snakemake = SimpleNamespace(
        input=SimpleNamespace(
            antismash=os.path.join(summary, "antismash.bgc.tsv"),
            gecco=os.path.join(summary, "gecco.bgc.tsv"),
            deepbgc=os.path.join(summary, "deepbgc.bgc.tsv"),
            arts=os.path.join(summary, "arts.hits.tsv"),
            bakta=os.path.join(args.workdir, "results", args.sample, "bakta",
                               args.sample + ".features.tsv"),
            fna=os.path.join(args.workdir, "data", "bakta", args.sample, args.sample + ".fna"),
        ),
        output=SimpleNamespace(
            consensus=os.path.join(args.out, "consensus_bgcs.tsv"),
            overlap=os.path.join(args.out, "tool_overlap.tsv"),
            gene_support=os.path.join(args.out, "gene_caller_support.tsv"),
        ),
        config={"consensus": {
            "mode": args.mode,
            "min_callers_per_gene": args.min_callers,
            "min_containment": args.min_containment,
        }},
    )
    sys.path.insert(0, os.path.join(ROOT, "scripts"))
    runpy.run_path(os.path.join(ROOT, "scripts", "build_consensus.py"),
                   init_globals={"snakemake": snakemake}, run_name="__main__")


if __name__ == "__main__":
    main()
