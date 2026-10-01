#!/usr/bin/env python3
"""Score BGC prediction sets against MIBiG-derived truth sets.

Prediction sets ("arms"), all from one run of the workflow per genome:
  antismash, gecco, deepbgc   each caller alone
  side_by_side                the three callers' rows concatenated without any
                              grouping - what nf-core/funcscan's comBGC reports
  <variant>:agreed            a consensus variant, interval = agreed gene span
  <variant>:union             the same loci, interval = union of the members
Consensus variants come from scripts/build_consensus.py itself (see
run_consensus_variant.py): voting with 1, 2 (the default) or 3 callers per
gene, and the previous containment rule.

Definitions, fixed before any result was inspected. T is a truth locus, P a
prediction on the same contig, |.| a length in bp. A gene belongs to an
interval when at least half of the CDS lies inside it.
  recovery category (BGC-QUAST thresholds):
    full        one prediction covers >= 95 % of T
    full_multi  several predictions together cover >= 95 % of T, none alone
    partial     predictions cover 10-95 % of T
    missed      predictions cover < 10 % of T
  detected      not missed
  best match    the prediction with the highest bp Jaccard with T
  delimited     best-match Jaccard >= 0.5 (a one-to-one recovery)
  split         no prediction covers >= 50 % of T, several together do
  fused         the best match also covers >= 50 % of another truth locus
                that is neither nested in T nor contains it
  left/right extension  signed bp by which the best match extends beyond T at
                each end (negative = truncated)
  gene precision / recall / F1  of the best match's genes against T's genes
  fragments     number of predictions sharing at least one gene with T

Precision at the locus level is not reported: MIBiG does not annotate every
cluster in these genomes, so a prediction without a MIBiG counterpart is an
uncharacterised candidate, not a false positive. Gene precision inside a
recovered locus is well defined and is reported instead.

Usage:
  validate_consensus.py --workdir work --variants variants --truth-dir truth \
      --samples S_coelicolor_A3_2,S_albidoflavus_J1074 --out metrics
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
from scipy import stats

CALLERS = ["antismash", "gecco", "deepbgc"]
VARIANTS = ["voting2", "voting1", "voting3", "containment"]
REFERENCE_ARM = "voting2:agreed"
BOOTSTRAP = 10000
SEED = 20260930


def to_num(series):
    return pd.to_numeric(series, errors="coerce")


def contig_map(workdir, sample):
    record = json.load(open(os.path.join(workdir, "data", "bakta", sample, sample + ".json")))
    return {s["id"]: s.get("orig_id", s["id"]) for s in record["sequences"]}


def load_genes(workdir, sample, cmap):
    genes = pd.read_csv(os.path.join(workdir, "results", sample, "bakta", sample + ".features.tsv"), sep="\t")
    genes = genes[genes["type"].astype(str) == "cds"].copy()
    genes["contig"] = genes["contig"].map(cmap)
    genes["start"], genes["end"] = to_num(genes["start"]), to_num(genes["end"])
    return genes[["contig", "start", "end", "locus_tag"]].reset_index(drop=True)


def load_arms(workdir, variants_dir, sample, cmap):
    summary = os.path.join(workdir, "results", sample, "summary")
    arms = {}
    frames = []
    for tool in CALLERS:
        df = pd.read_csv(os.path.join(summary, tool + ".bgc.tsv"), sep="\t")
        df = pd.DataFrame({"contig": df["contig"].astype(str).map(cmap),
                           "start": to_num(df["start"]), "end": to_num(df["end"]), "tool": tool})
        arms[tool] = df.dropna(subset=["start", "end"])
        frames.append(arms[tool])
    arms["side_by_side"] = pd.concat(frames, ignore_index=True)
    for variant in VARIANTS:
        path = os.path.join(variants_dir, sample, variant, "consensus_bgcs.tsv")
        if not os.path.exists(path):
            continue
        df = pd.read_csv(path, sep="\t")
        contig = df["contig"].astype(str).map(cmap)
        support = df["support_tools"].astype(str)
        arms[variant + ":agreed"] = pd.DataFrame({
            "contig": contig, "start": to_num(df["start"]), "end": to_num(df["end"]),
            "tool": support}).dropna(subset=["start", "end"])
        arms[variant + ":union"] = pd.DataFrame({
            "contig": contig, "start": to_num(df["union_start"]), "end": to_num(df["union_end"]),
            "tool": support}).dropna(subset=["start", "end"])
        # Post hoc: the multi-caller tier alone, agreed span.
        multi = to_num(df["support_count"]) >= 2
        arms[variant + ":multi"] = arms[variant + ":agreed"][multi.loc[arms[variant + ":agreed"].index]]
    for name, df in arms.items():
        if df["contig"].isna().any():
            raise SystemExit("{0}/{1}: predictions on contigs absent from the Bakta contig map".format(sample, name))
        df["start"], df["end"] = df["start"].astype(int), df["end"].astype(int)
    return arms


def overlap(a1, a2, b1, b2):
    return max(0, min(a2, b2) - max(a1, b1) + 1)


def genes_in(genes, contig, start, end):
    """Locus tags of CDSs with at least half their length inside [start, end]."""
    local = genes[(genes.contig == contig) & (genes.start <= end) & (genes.end >= start)]
    inside = (np.minimum(local.end, end) - np.maximum(local.start, start) + 1) / (local.end - local.start + 1)
    return set(local.locus_tag[inside >= 0.5])


def union_cover(intervals, start, end):
    clipped = sorted((max(s, start), min(e, end)) for s, e in intervals if e >= start and s <= end)
    total, current = 0, None
    for s, e in clipped:
        if current is None or s > current[1] + 1:
            if current:
                total += current[1] - current[0] + 1
            current = [s, e]
        else:
            current[1] = max(current[1], e)
    if current:
        total += current[1] - current[0] + 1
    return total


def related(a, b):
    """True when one truth locus is nested in the other."""
    return a["nested_in"] == b["mibig_id"] or b["nested_in"] == a["mibig_id"]


def score_arm(preds, truth, genes):
    rows = []
    for _, t in truth.iterrows():
        tlen = t.end - t.start + 1
        tgenes = genes_in(genes, t.contig, t.start, t.end)
        local = preds[(preds.contig == t.contig) & (preds.start <= t.end) & (preds.end >= t.start)]
        cover_any = union_cover(list(zip(local.start, local.end)), t.start, t.end) / tlen
        best, cover_single, fragments, fragment_jaccards = None, 0.0, 0, []
        for _, p in local.iterrows():
            ov = overlap(t.start, t.end, p.start, p.end)
            cover_single = max(cover_single, ov / tlen)
            pgenes = genes_in(genes, p.contig, p.start, p.end)
            jac = ov / (max(t.end, p.end) - min(t.start, p.start) + 1)
            if pgenes & tgenes:
                fragments += 1
                fragment_jaccards.append(jac)
            if best is None or jac > best["jaccard"]:
                best = {"jaccard": jac, "start": p.start, "end": p.end, "genes": pgenes}
        if cover_single >= 0.95:
            category = "full"
        elif cover_any >= 0.95:
            category = "full_multi"
        elif cover_any >= 0.10:
            category = "partial"
        else:
            category = "missed"
        row = {"mibig_id": t.mibig_id, "compound": t.compound, "truth_length": tlen,
               "truth_genes": len(tgenes), "nested_in": t.nested_in,
               "category": category, "detected": category != "missed",
               "cover_single": round(cover_single, 3), "cover_any": round(cover_any, 3),
               "split": cover_single < 0.5 and cover_any >= 0.5 and len(local) >= 2,
               "fragments": fragments,
               # Post hoc: what a user gets by taking any one of the overlapping rows.
               "jaccard_any_fragment": round(float(np.mean(fragment_jaccards)), 4) if fragment_jaccards else 0.0}
        if best is not None and category != "missed":
            shared = len(best["genes"] & tgenes)
            precision = shared / len(best["genes"]) if best["genes"] else 0.0
            recall = shared / len(tgenes) if tgenes else 0.0
            fused = any(
                overlap(o.start, o.end, best["start"], best["end"]) / (o.end - o.start + 1) >= 0.5
                for _, o in truth[(truth.contig == t.contig) & (truth.mibig_id != t.mibig_id)].iterrows()
                if not related(t, o))
            row.update({
                "jaccard": round(best["jaccard"], 4),
                "delimited": best["jaccard"] >= 0.5,
                "pred_length": best["end"] - best["start"] + 1,
                "length_ratio": round((best["end"] - best["start"] + 1) / tlen, 3),
                "left_extension": int(t.start - best["start"]),
                "right_extension": int(best["end"] - t.end),
                "gene_precision": round(precision, 4),
                "gene_recall": round(recall, 4),
                "gene_f1": round(2 * precision * recall / (precision + recall), 4) if precision + recall else 0.0,
                "fused": fused,
            })
        else:
            row.update({"jaccard": 0.0, "delimited": False, "pred_length": 0, "length_ratio": np.nan,
                        "left_extension": np.nan, "right_extension": np.nan, "gene_precision": np.nan,
                        "gene_recall": 0.0, "gene_f1": 0.0, "fused": False})
        rows.append(row)
    return pd.DataFrame(rows)


def genome_level(preds, truth, genome_length):
    covered = sum(union_cover(list(zip(g.start, g.end)), 1, 10 ** 10) for _, g in preds.groupby("contig"))
    touches = 0
    for _, p in preds.iterrows():
        local = truth[truth.contig == p.contig]
        if any(overlap(p.start, p.end, t.start, t.end) > 0 for _, t in local.iterrows()):
            touches += 1
    return {"predictions": len(preds), "bp_flagged": covered,
            "genome_fraction": round(covered / genome_length, 4),
            "predictions_touching_truth": touches,
            "predictions_without_mibig": len(preds) - touches}


def summarise(df):
    det = df[df.detected]
    abs_ext = pd.concat([det.left_extension.abs(), det.right_extension.abs()])
    return {
        "truth_loci": len(df),
        "detected": int(df.detected.sum()),
        "full": int((df.category == "full").sum()),
        "full_multi": int((df.category == "full_multi").sum()),
        "partial": int((df.category == "partial").sum()),
        "missed": int((df.category == "missed").sum()),
        "delimited": int(df.delimited.sum()),
        "split": int(df.split.sum()),
        "fused": int(df.fused.sum()),
        "median_jaccard": round(df.jaccard.median(), 3),
        "mean_jaccard": round(df.jaccard.mean(), 3),
        "median_gene_f1": round(df.gene_f1.median(), 3),
        "median_gene_precision": round(det.gene_precision.median(), 3) if len(det) else np.nan,
        "median_gene_recall": round(df.gene_recall.median(), 3),
        "median_length_ratio": round(det.length_ratio.median(), 2) if len(det) else np.nan,
        "median_abs_boundary_error_bp": int(abs_ext.median()) if len(abs_ext) else np.nan,
        "median_fragments_detected": float(det.fragments.median()) if len(det) else np.nan,
        "median_jaccard_any_fragment": round(df.jaccard_any_fragment.median(), 3),
    }


def paired_tests(per_cluster, metric):
    """Reference arm against every other arm, paired by truth locus, pooled over genomes."""
    rows = []
    rng = np.random.default_rng(SEED)
    wide = per_cluster.pivot_table(index=["sample", "mibig_id"], columns="arm", values=metric)
    if REFERENCE_ARM not in wide:
        return pd.DataFrame()
    for arm in wide.columns:
        if arm == REFERENCE_ARM:
            continue
        pair = wide[[REFERENCE_ARM, arm]].dropna()
        diff = (pair[REFERENCE_ARM] - pair[arm]).to_numpy()
        samples = pair.index.get_level_values("sample").to_numpy()
        boot = []
        for _ in range(BOOTSTRAP):
            idx = np.concatenate([rng.choice(np.flatnonzero(samples == s), size=(samples == s).sum())
                                  for s in np.unique(samples)])
            boot.append(np.mean(diff[idx]))
        nonzero = diff[diff != 0]
        p = stats.wilcoxon(nonzero).pvalue if len(nonzero) >= 5 else np.nan
        rows.append({"metric": metric, "reference": REFERENCE_ARM, "other": arm, "n": len(diff),
                     "reference_better": int((diff > 0).sum()), "other_better": int((diff < 0).sum()),
                     "ties": int((diff == 0).sum()), "mean_difference": round(diff.mean(), 4),
                     "ci95_low": round(np.percentile(boot, 2.5), 4),
                     "ci95_high": round(np.percentile(boot, 97.5), 4),
                     "wilcoxon_p": p})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--variants", required=True)
    ap.add_argument("--truth-dir", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    per_cluster, genome_rows = [], []
    for sample in args.samples.split(","):
        cmap = contig_map(args.workdir, sample)
        truth = pd.read_csv(os.path.join(args.truth_dir, sample + ".tsv"), sep="\t").fillna({"nested_in": ""})
        genes = load_genes(args.workdir, sample, cmap)
        genome_length = sum(s["length"] for s in json.load(open(os.path.join(
            args.workdir, "data", "bakta", sample, sample + ".json")))["sequences"])
        for arm, preds in load_arms(args.workdir, args.variants, sample, cmap).items():
            df = score_arm(preds, truth, genes)
            df.insert(0, "arm", arm)
            df.insert(0, "sample", sample)
            per_cluster.append(df)
            genome_rows.append(dict(sample=sample, arm=arm, **genome_level(preds, truth, genome_length)))

    per_cluster = pd.concat(per_cluster, ignore_index=True)
    per_cluster.to_csv(os.path.join(args.out, "per_cluster.tsv"), sep="\t", index=False)
    pd.DataFrame(genome_rows).to_csv(os.path.join(args.out, "genome_level.tsv"), sep="\t", index=False)

    by_genome = [dict(sample=s, arm=a, **summarise(g)) for (s, a), g in per_cluster.groupby(["sample", "arm"])]
    pd.DataFrame(by_genome).to_csv(os.path.join(args.out, "summary_by_genome.tsv"), sep="\t", index=False)

    subsets = {
        "all": per_cluster,
        "non_nested_ge5kb": per_cluster[(per_cluster.nested_in == "") & (per_cluster.truth_length >= 5000)],
    }
    pooled, tests = [], []
    for subset, frame in subsets.items():
        for arm, g in frame.groupby("arm"):
            pooled.append(dict(subset=subset, arm=arm, **summarise(g)))
        for metric in ("jaccard", "gene_f1", "jaccard_any_fragment"):
            t = paired_tests(frame, metric)
            if len(t):
                t.insert(0, "subset", subset)
                tests.append(t)
    pooled = pd.DataFrame(pooled)
    pooled.to_csv(os.path.join(args.out, "summary_pooled.tsv"), sep="\t", index=False)
    if tests:
        pd.concat(tests).to_csv(os.path.join(args.out, "paired_tests.tsv"), sep="\t", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(pooled[pooled.subset == "all"].drop(columns="subset").to_string(index=False))


if __name__ == "__main__":
    main()
