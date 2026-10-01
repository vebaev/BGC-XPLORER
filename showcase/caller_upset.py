#!/usr/bin/env python3
"""UpSet plot of which callers support each locus in a finished run.

Bars count loci per exact caller combination; the dots below name the
combination. Drawn with matplotlib alone so it runs in the validation
environment without an UpSet package.

Usage:
  caller_upset.py --results results --sample Soil_1 --out showcase/Soil_1_upset.png
  (a .pdf path writes a vector figure)
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

CALLERS = [("antismash", "antiSMASH"), ("gecco", "GECCO"), ("deepbgc", "DeepBGC")]
COLOURS = {3: "#2b6cb0", 2: "#4a9fd8", 1: "#a0aec0"}
# GECCO and DeepBGC are both Pfam-based, so their agreement alone is weaker
# evidence than two callers suggest; the bar is set apart.
PFAM_ONLY = frozenset({"gecco", "deepbgc"})
PFAM_ONLY_COLOUR = "#dd8452"


def combinations(results, sample):
    loci = pd.read_csv(os.path.join(results, sample, "summary", "region_evidence.tsv"), sep="\t")
    counts = {}
    for value in loci["support_tools"].fillna(""):
        tools = frozenset(tool.strip().lower() for tool in str(value).split(",") if tool.strip())
        if tools:
            counts[tools] = counts.get(tools, 0) + 1
    # Most callers first, then most loci, so agreement reads left to right.
    return sorted(counts.items(), key=lambda item: (-len(item[0]), -item[1])), len(loci)


def draw(combos, total, sample, out):
    fig, (bars, dots) = plt.subplots(2, 1, figsize=(1.0 + 0.62 * len(combos), 4.2),
                                     gridspec_kw={"height_ratios": [2.4, 1]}, sharex=True)
    xs = range(len(combos))
    heights = [count for _, count in combos]
    colours = [PFAM_ONLY_COLOUR if tools == PFAM_ONLY else COLOURS[len(tools)] for tools, _ in combos]
    bars.bar(xs, heights, color=colours, width=0.6)
    for x, height in zip(xs, heights):
        bars.text(x, height + 0.3, str(height), ha="center", va="bottom", fontsize=8)
    bars.set_ylabel("Loci", fontsize=9)
    bars.set_ylim(0, max(heights) * 1.18)
    bars.set_title("{0}: caller support of {1} loci".format(sample, total), fontsize=10)
    for spine in ("top", "right"):
        bars.spines[spine].set_visible(False)
    bars.tick_params(axis="x", length=0)
    if any(tools == PFAM_ONLY for tools, _ in combos):
        bars.text(0.99, 0.97, "orange: GECCO + DeepBGC only (both Pfam-based)", transform=bars.transAxes,
                  ha="right", va="top", fontsize=7.5, color="#9c4a1a")

    for row, (key, _) in enumerate(CALLERS):
        y = len(CALLERS) - 1 - row
        for x, (tools, _) in enumerate(combos):
            dots.scatter(x, y, s=46, color="#2d3748" if key in tools else "#e2e8f0", zorder=3)
    for x, (tools, _) in enumerate(combos):
        ys = [len(CALLERS) - 1 - row for row, (key, _) in enumerate(CALLERS) if key in tools]
        if len(ys) > 1:
            dots.plot([x, x], [min(ys), max(ys)], color="#2d3748", linewidth=1.6, zorder=2)
    dots.set_yticks(range(len(CALLERS)))
    dots.set_yticklabels([label for _, label in reversed(CALLERS)], fontsize=8.5)
    dots.set_ylim(-0.6, len(CALLERS) - 0.4)
    dots.set_xticks([])
    for spine in dots.spines.values():
        spine.set_visible(False)
    dots.tick_params(length=0)
    fig.tight_layout(h_pad=0.2)
    fig.savefig(out, dpi=300)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    combos, total = combinations(args.results, args.sample)
    draw(combos, total, args.sample, args.out)
    for tools, count in combos:
        print("{0}\t{1}".format(" + ".join(label for key, label in CALLERS if key in tools), count))


if __name__ == "__main__":
    main()
