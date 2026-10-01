#!/usr/bin/env python3
"""Re-parse finished DeepBGC runs with the current scripts/parse_deepbgc.py.

Used once, when v1.2.0 fixed the 1 bp offset in DeepBGC start coordinates
(DeepBGC writes 0-based starts; the parser now converts them). The raw DeepBGC
output is unchanged, so only the parsed table and the consensus built from it
need regenerating; this runs the pipeline's own parser and grouping on the raw
output of the original workflow runs.

Usage:
  reparse_deepbgc.py --runs ../benchmark/work --workdir work \
      --samples S_coelicolor_A3_2,S_albidoflavus_J1074,...
"""
import argparse
import os
import runpy
import sys
from types import SimpleNamespace

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", required=True, help="workflow work directory holding results/<sample>/deepbgc")
    ap.add_argument("--workdir", required=True, help="validation work directory to update")
    ap.add_argument("--samples", required=True)
    args = ap.parse_args()
    sys.path.insert(0, os.path.join(ROOT, "scripts"))
    workdir = os.path.abspath(args.workdir)
    here = os.getcwd()
    for sample in args.samples.split(","):
        out = os.path.join(workdir, "results", sample, "summary", "deepbgc.bgc.tsv")
        # Paths relative to the run directory, as Snakemake would pass them,
        # so the source_file column matches the workflow's.
        os.chdir(args.runs)
        snakemake = SimpleNamespace(
            wildcards=SimpleNamespace(sample=sample),
            input=SimpleNamespace(done=os.path.join("results", sample, "deepbgc", ".done")),
            output=[out],
        )
        runpy.run_path(os.path.join(ROOT, "scripts", "parse_deepbgc.py"),
                       init_globals={"snakemake": snakemake}, run_name="__main__")
        os.chdir(here)
        print("re-parsed", sample)


if __name__ == "__main__":
    main()
