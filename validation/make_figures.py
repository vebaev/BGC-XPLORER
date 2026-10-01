#!/usr/bin/env python3
"""Figures for the benchmark, from metrics/per_cluster.tsv.

  recovery.png   recovery categories per prediction set, all truth loci pooled
  jaccard.png    per-locus Jaccard with the truth locus, paired across sets
  boundaries.png signed extension beyond the truth locus at each end
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# The union span is the gene-map window, not the reported locus; it stays in
# the tables but not in the figures. "side_by_side:any" is the same rows as
# side_by_side, scored by the mean over the rows overlapping the locus - what a
# user gets without knowing which row is right. The plain side_by_side score
# picks the best row with knowledge of the truth, an unattainable upper bound.
ARMS = [
    ("antismash", "antiSMASH"),
    ("gecco", "GECCO"),
    ("deepbgc", "DeepBGC"),
    ("side_by_side:any", "Side-by-side,\nany row"),
    ("side_by_side", "Side-by-side,\nbest row\n(oracle)"),
    ("containment:agreed", "Containment\n(previous)"),
    ("voting2:multi", "Voting 2/3,\n≥2-caller\ntier only"),
    ("voting2:agreed", "Voting 2/3\n(default)"),
    ("voting3:agreed", "Voting 3/3"),
]
RECOVERY_ARMS = [a for a in ARMS if a[0] != "side_by_side:any"]
CATEGORIES = [("full", "Full (one prediction ≥95 %)", "#2b6cb0"),
              ("full_multi", "Full, several predictions", "#90cdf4"),
              ("partial", "Partial (10–95 %)", "#f6ad55"),
              ("missed", "Missed (<10 %)", "#e2e8f0")]


def recovery(df, out):
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    present = [(arm, label.replace("\n", " ")) for arm, label in RECOVERY_ARMS if arm in df.arm.unique()]
    arms = [arm for arm, _ in present]
    labels = [label for _, label in present]
    left = np.zeros(len(arms))
    for key, name, colour in CATEGORIES:
        counts = np.array([(df[df.arm == a].category == key).sum() for a in arms])
        ax.barh(labels, counts, left=left, color=colour, label=name, edgecolor="white")
        left += counts
    delimited = [df[df.arm == a].delimited.sum() for a in arms]
    for y, d in enumerate(delimited):
        ax.text(left[y] + 0.5, y, "{0} delimited".format(d), va="center", fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("MIBiG truth loci (n = {0}, four genomes)".format(len(df[df.arm == arms[0]])))
    ax.set_xlim(0, left.max() * 1.3)
    ax.tick_params(axis="y", labelsize=8)
    ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.45, -0.14), ncol=4, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=300)


def jaccard(df, out):
    any_row = df[df.arm == "side_by_side"].assign(arm="side_by_side:any", jaccard=lambda d: d.jaccard_any_fragment)
    df = pd.concat([df, any_row])
    arms = [arm for arm, _ in ARMS if arm in df.arm.unique()]
    labels = [label for arm, label in ARMS if arm in df.arm.unique()]
    wide = df.pivot_table(index=["sample", "mibig_id"], columns="arm", values="jaccard")[arms]
    fig, ax = plt.subplots(figsize=(9, 4.4))
    for _, row in wide.iterrows():
        ax.plot(range(len(arms)), row.values, color="#cbd5e0", lw=0.6, zorder=1)
    ax.boxplot([wide[a].values for a in arms], positions=range(len(arms)), widths=0.45,
               showfliers=False, medianprops={"color": "#c53030", "lw": 2}, zorder=2)
    ax.set_xticks(range(len(arms)))
    ax.set_xticklabels(labels, fontsize=7.5)
    ax.set_ylabel("Jaccard with the MIBiG locus (bp)")
    ax.set_ylim(-0.02, 1.02)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=300)


def boundaries(df, out):
    arms = [arm for arm, _ in RECOVERY_ARMS if arm in df.arm.unique()]
    labels = [label for arm, label in RECOVERY_ARMS if arm in df.arm.unique()]
    fig, ax = plt.subplots(figsize=(9, 4.4))
    data = []
    for a in arms:
        sub = df[(df.arm == a) & df.detected]
        data.append(np.concatenate([sub.left_extension.values, sub.right_extension.values]) / 1000.0)
    ax.axhline(0, color="#718096", lw=0.8)
    ax.boxplot(data, positions=range(len(arms)), widths=0.45, showfliers=True,
               flierprops={"markersize": 2}, medianprops={"color": "#c53030", "lw": 2})
    ax.set_yscale("symlog", linthresh=1)
    ax.set_xticks(range(len(arms)))
    ax.set_xticklabels(labels, fontsize=7.5)
    ax.set_ylabel("Extension beyond the MIBiG locus, kb\n(each end; negative = truncated)")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=300)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metrics", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    df = pd.read_csv(os.path.join(args.metrics, "per_cluster.tsv"), sep="\t")
    recovery(df, os.path.join(args.out, "recovery.png"))
    jaccard(df, os.path.join(args.out, "jaccard.png"))
    boundaries(df, os.path.join(args.out, "boundaries.png"))


if __name__ == "__main__":
    main()
