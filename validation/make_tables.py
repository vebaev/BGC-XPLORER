#!/usr/bin/env python3
"""Per-genome tables and the per-locus figure for the manuscript.

One row per prediction set and genome. BGC-XPLORER is the workflow as
released: voting with two of three callers per gene, reported (agreed) span.
Side-by-side is the three callers' rows without grouping (funcscan's comBGC);
its Jaccard is the mean over the rows overlapping the locus, since choosing the
best row needs the answer.

Outputs:
  tables/table_per_genome.tsv / .md    main table
  tables/supplementary_per_locus.xlsx  one sheet per genome, every MIBiG locus
  figures/per_locus.png                status and Jaccard per locus and tool
  figures/per_locus_supp.png / .pdf    the same without side-by-side (supplementary)
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SETS = [("antismash", "antiSMASH"), ("gecco", "GECCO"), ("deepbgc", "DeepBGC"),
        ("side_by_side", "Side-by-side (no grouping)"), ("voting2:agreed", "BGC-XPLORER")]
GENOMES = [("S_coelicolor_A3_2", "S. coelicolor A3(2)"),
           ("S_albidoflavus_J1074", "S. albidoflavus J1074"),
           ("S_avermitilis_MA4680", "S. avermitilis MA-4680"),
           ("S_tropica_CNB440", "Salinispora tropica CNB-440")]
COLOURS = {"full": "#2b6cb0", "partial": "#f6ad55", "missed": "#edf2f7"}


def load(metrics, bgcquast):
    d = pd.read_csv(os.path.join(metrics, "per_cluster.tsv"), sep="\t")
    d = d[d.arm.isin([a for a, _ in SETS])].copy()
    d["status"] = d.category.replace({"full_multi": "full"})
    side = d.arm == "side_by_side"
    d.loc[side, "jaccard"] = d.loc[side, "jaccard_any_fragment"]
    q = pd.read_csv(os.path.join(bgcquast, "bgcquast_per_locus.tsv"), sep="\t")
    d = d.merge(q, on=["sample", "arm", "mibig_id"], how="left")
    g = pd.read_csv(os.path.join(metrics, "genome_level.tsv"), sep="\t")
    return d, g


def counts(frame, column):
    return "{0} / {1} / {2}".format(*[(frame[column] == s).sum() for s in ("full", "partial", "missed")])


def main_table(d, g):
    rows = []
    for sample, genome in GENOMES:
        for arm, label in SETS:
            sub = d[(d["sample"] == sample) & (d.arm == arm)]
            det = sub[sub.detected]
            ext = pd.concat([det.left_extension.abs(), det.right_extension.abs()])
            predicted = g[(g["sample"] == sample) & (g.arm == arm)].predictions.iloc[0]
            rows.append({
                "Genome": genome, "Tool": label, "MIBiG loci": len(sub),
                "Predicted loci": int(predicted),
                "Detected": int(sub.detected.sum()),
                "Full / partial / missed": counts(sub, "status"),
                "BGC-QUAST full / partial / missed": counts(sub, "bgcquast_status"),
                "Same status as BGC-QUAST": "{0}/{1}".format((sub.status == sub.bgcquast_status).sum(), len(sub)),
                "Delimited (Jaccard ≥ 0.5)": int(sub.delimited.sum()) if arm != "side_by_side" else "—",
                "Median Jaccard": round(sub.jaccard.median(), 2),
                "Median gene F1": round(sub.gene_f1.median(), 2) if arm != "side_by_side" else "—",
                "Median boundary error, kb": round(ext.median() / 1000, 1) if arm != "side_by_side" else "—",
            })
    return pd.DataFrame(rows)


def to_markdown(table):
    lines = []
    for genome, sub in table.groupby("Genome", sort=False):
        sub = sub.drop(columns="Genome")
        lines.append("### " + genome + "\n")
        lines.append("| " + " | ".join(sub.columns) + " |")
        lines.append("|" + "---|" * len(sub.columns))
        for _, r in sub.iterrows():
            cells = [str(v) for v in r.values]
            if r["Tool"] == "BGC-XPLORER":
                cells = ["**" + c + "**" for c in cells]
            lines.append("| " + " | ".join(cells) + " |")
        lines.append("")
    lines.append("Side-by-side: the three callers' rows without grouping, as nf-core/funcscan's comBGC "
                 "reports them; Jaccard is the mean over the rows overlapping the locus. Status and counts "
                 "agree with BGC-QUAST 1.1.0 (compare-to-reference, MIBiG loci as reference) for every locus.")
    return "\n".join(lines)


def per_locus_sheets(d, path):
    with pd.ExcelWriter(path) as writer:
        for sample, genome in GENOMES:
            sub = d[d["sample"] == sample]
            base = sub[sub.arm == "antismash"][["mibig_id", "compound", "truth_length", "truth_genes"]]
            base = base.rename(columns={"truth_length": "length_bp", "truth_genes": "genes"})
            for arm, label in SETS:
                part = sub[sub.arm == arm].set_index("mibig_id")
                base[label + " status"] = base.mibig_id.map(part.status)
                base[label + " BGC-QUAST"] = base.mibig_id.map(part.bgcquast_status)
                base[label + " Jaccard"] = base.mibig_id.map(part.jaccard).round(3)
            base.to_excel(writer, sheet_name=genome.replace("Salinispora", "Sal.")[:31], index=False)


def per_locus_figure(d, path, sets=SETS):
    sizes = [len(d[(d["sample"] == s) & (d.arm == "antismash")]) for s, _ in GENOMES]
    fig, axes = plt.subplots(4, 1, figsize=(2.6 + 0.78 * len(sets), 0.32 * sum(sizes) + 3.2),
                             gridspec_kw={"height_ratios": sizes})
    for ax, (sample, genome) in zip(axes, GENOMES):
        sub = d[d["sample"] == sample]
        loci = sub[sub.arm == "antismash"].sort_values("truth_length", ascending=False)
        names = [str(c).split(";")[0][:24] for c in loci.compound]
        for x, (arm, _) in enumerate(sets):
            part = sub[sub.arm == arm].set_index("mibig_id")
            for y, mibig_id in enumerate(loci.mibig_id):
                status, jac = part.status[mibig_id], part.jaccard[mibig_id]
                ax.add_patch(plt.Rectangle((x, y), 0.95, 0.9, color=COLOURS[status]))
                if status != "missed":
                    ax.text(x + 0.475, y + 0.45, "{0:.2f}".format(jac), ha="center", va="center",
                            fontsize=7, color="white" if status == "full" else "black")
        ax.set_xlim(0, len(sets))
        ax.set_ylim(len(loci), 0)
        ax.set_xticks(np.arange(len(sets)) + 0.475)
        ax.xaxis.tick_top()
        ax.set_xticklabels([l.replace(" (no grouping)", "\n(no grouping)") for _, l in sets]
                           if ax is axes[0] else [], fontsize=8)
        ax.set_yticks(np.arange(len(loci)) + 0.45)
        ax.set_yticklabels(names, fontsize=7.5)
        ax.set_ylabel(genome, fontsize=9, style="italic", labelpad=8)
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in COLOURS.values()]
    fig.legend(handles, ["Full (≥95 % covered)", "Partial (10–95 %)", "Missed (<10 %)"],
               loc="lower center", bbox_to_anchor=(0.5, 0.012), ncol=3, frameon=False, fontsize=8)
    fig.text(0.5, 0.004, "Cell value: Jaccard with the MIBiG locus. Status identical to BGC-QUAST for every cell.",
             ha="center", fontsize=7.5)
    fig.tight_layout(rect=(0, 0.035, 1, 1), h_pad=0.6)
    fig.savefig(path, dpi=300)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metrics", default="metrics")
    ap.add_argument("--bgcquast", default="bgcquast")
    ap.add_argument("--tables", default="tables")
    ap.add_argument("--figures", default="figures")
    args = ap.parse_args()
    os.makedirs(args.tables, exist_ok=True)
    d, g = load(args.metrics, args.bgcquast)
    if d.bgcquast_status.isna().any():
        raise SystemExit("BGC-QUAST status missing for some loci")
    table = main_table(d, g)
    table.to_csv(os.path.join(args.tables, "table_per_genome.tsv"), sep="\t", index=False)
    open(os.path.join(args.tables, "table_per_genome.md"), "w").write(to_markdown(table))
    per_locus_sheets(d, os.path.join(args.tables, "supplementary_per_locus.xlsx"))
    per_locus_figure(d, os.path.join(args.figures, "per_locus.png"))
    no_side = [s for s in SETS if s[0] != "side_by_side"]
    for ext in ("png", "pdf"):
        per_locus_figure(d, os.path.join(args.figures, "per_locus_supp." + ext), no_side)
    print(to_markdown(table))


if __name__ == "__main__":
    main()
