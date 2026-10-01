#!/usr/bin/env python3
"""Build a MIBiG-derived truth set for one reference genome.

Coordinates are not copied from MIBiG. Every active MIBiG 4.0 entry is aligned
to the genome with megablast, and an entry is accepted as present when its
locus sequence aligns at >= MIN_IDENTITY over >= MIN_COVERAGE of its length in
one collinear chain. The truth interval is the span of that chain on the
genome. This recovers entries deposited under another accession (INSDC vs
RefSeq, a cosmid, a sister strain) and never trusts a coordinate that belongs
to a superseded assembly.

Entries that land on the same place are collapsed into one truth locus
(reciprocal overlap >= DUPLICATE_OVERLAP, or one inside the other with a
compound in common - the same cluster deposited twice). The representative is the entry
deposited on this genome's own record, then the higher MIBiG quality, then the
complete one, then the lowest accession; the others are kept as aliases.
Entries lying inside a larger truth locus are kept and flagged, because MIBiG
records them as separate clusters.

Usage:
  build_truth.py --genome work/data/fasta/S_coelicolor_A3_2.fasta \
                 --accessions NC_003888.3,AL645882.2 \
                 --mibig-json mibig/mibig_json_4.0 --mibig-gbk mibig/mibig_gbk_4.0 \
                 --out truth/S_coelicolor_A3_2.tsv
"""
import argparse
import glob
import json
import os
import subprocess
import tempfile

import pandas as pd
from Bio import SeqIO

MIN_IDENTITY = 95.0
MIN_COVERAGE = 0.90
MAX_SPAN_FACTOR = 1.5
DUPLICATE_OVERLAP = 0.80
NESTED_FRACTION = 0.90
QUALITY_RANK = {"high": 0, "medium": 1, "questionable": 2}

BLAST_FIELDS = ["qseqid", "sseqid", "pident", "length", "qstart", "qend",
                "sstart", "send", "evalue", "bitscore", "qlen", "slen"]


def load_mibig(json_dir):
    """Active entries only; retired and pending records are not curated truth."""
    entries = {}
    for path in glob.glob(os.path.join(json_dir, "*.json")):
        record = json.load(open(path))
        if record.get("status") != "active":
            continue
        compounds = [c.get("name", "") for c in record.get("compounds", []) if c.get("name")]
        classes = sorted({c.get("class", "") for c in record.get("biosynthesis", {}).get("classes", [])})
        loci = record.get("loci", [])
        entries[record["accession"]] = {
            "mibig_id": record["accession"],
            "compound": ";".join(compounds),
            "class": ";".join(c for c in classes if c),
            "quality": record.get("quality", ""),
            "completeness": record.get("completeness", ""),
            "source_accession": loci[0].get("accession", "") if loci else "",
            "organism": record.get("taxonomy", {}).get("name", ""),
        }
    return entries


def write_queries(entries, gbk_dir, path):
    written = 0
    with open(path, "w") as handle:
        for mibig_id in sorted(entries):
            gbk = os.path.join(gbk_dir, mibig_id + ".gbk")
            if not os.path.exists(gbk):
                continue
            record = next(SeqIO.parse(gbk, "genbank"))
            if len(record.seq) == 0:
                continue
            handle.write(">{0}\n{1}\n".format(mibig_id, str(record.seq)))
            written += 1
    return written


def run_blast(queries, genome, threads, workdir):
    db = os.path.join(workdir, "genome")
    subprocess.run(["makeblastdb", "-in", genome, "-dbtype", "nucl", "-out", db],
                   check=True, stdout=subprocess.DEVNULL)
    out = os.path.join(workdir, "hits.tsv")
    subprocess.run([
        "blastn", "-task", "megablast", "-query", queries, "-db", db,
        "-perc_identity", "90", "-evalue", "1e-20", "-max_target_seqs", "10",
        "-max_hsps", "500", "-num_threads", str(threads),
        "-outfmt", "6 " + " ".join(BLAST_FIELDS), "-out", out,
    ], check=True)
    return pd.read_csv(out, sep="\t", names=BLAST_FIELDS)


def union_length(intervals):
    total, current = 0, None
    for start, end in sorted(intervals):
        if current is None or start > current[1] + 1:
            if current:
                total += current[1] - current[0] + 1
            current = [start, end]
        else:
            current[1] = max(current[1], end)
    if current:
        total += current[1] - current[0] + 1
    return total


def best_chain(hsps, qlen):
    """The densest set of HSPs on one contig and strand within 1.5x the query length."""
    best = None
    for (contig, strand), group in hsps.groupby(["sseqid", "strand"]):
        group = group.sort_values("s_lo").reset_index(drop=True)
        for i in range(len(group)):
            window = group[(group.s_lo >= group.s_lo[i])
                           & (group.s_hi <= group.s_lo[i] + MAX_SPAN_FACTOR * qlen)]
            covered = union_length(list(zip(window.qstart, window.qend)))
            identity = (window.pident * window.length).sum() / window.length.sum()
            candidate = (covered / qlen, identity, contig, int(window.s_lo.min()),
                         int(window.s_hi.max()), strand)
            if best is None or candidate[:2] > best[:2]:
                best = candidate
    return best


def reciprocal_overlap(a, b):
    ov = max(0, min(a["end"], b["end"]) - max(a["start"], b["start"]) + 1)
    return min(ov / a["length"], ov / b["length"]), ov / a["length"]


def compounds(row):
    return {c.strip().lower() for c in row["compound"].split(";") if c.strip()}


def collapse(rows, own_accessions):
    rows = sorted(rows, key=lambda r: (
        r["source_accession"].split(".")[0] not in own_accessions,
        QUALITY_RANK.get(r["quality"], 3),
        r["completeness"] != "complete",
        r["mibig_id"],
    ))
    kept = []
    for row in rows:
        duplicate_of = None
        for keeper in kept:
            if keeper["contig"] != row["contig"]:
                continue
            reciprocal, row_inside = reciprocal_overlap(row, keeper)
            keeper_inside = reciprocal_overlap(keeper, row)[1]
            same_compound = bool(compounds(row) & compounds(keeper))
            if reciprocal >= DUPLICATE_OVERLAP or (
                    same_compound and max(row_inside, keeper_inside) >= NESTED_FRACTION):
                duplicate_of = keeper
                break
        if duplicate_of is not None:
            duplicate_of["aliases"].append(row["mibig_id"])
        else:
            row["aliases"] = []
            kept.append(row)
    for row in kept:
        row["nested_in"] = ""
        for other in kept:
            if other is row or other["contig"] != row["contig"] or other["length"] <= row["length"]:
                continue
            if reciprocal_overlap(row, other)[1] >= NESTED_FRACTION:
                row["nested_in"] = other["mibig_id"]
                break
        row["aliases"] = ",".join(row["aliases"])
    return sorted(kept, key=lambda r: (r["contig"], r["start"]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--accessions", default="",
                    help="comma-separated accessions of this genome's own records (INSDC and RefSeq)")
    ap.add_argument("--mibig-json", required=True)
    ap.add_argument("--mibig-gbk", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    own = {a.split(".")[0] for a in args.accessions.split(",") if a}
    entries = load_mibig(args.mibig_json)
    with tempfile.TemporaryDirectory() as workdir:
        queries = os.path.join(workdir, "mibig.fna")
        n = write_queries(entries, args.mibig_gbk, queries)
        hits = run_blast(queries, args.genome, args.threads, workdir)

    hits["strand"] = (hits.send >= hits.sstart).map({True: "+", False: "-"})
    hits["s_lo"] = hits[["sstart", "send"]].min(axis=1)
    hits["s_hi"] = hits[["sstart", "send"]].max(axis=1)
    hits = hits[hits.pident >= MIN_IDENTITY - 5]

    accepted, near_misses = [], []
    for mibig_id, hsps in hits.groupby("qseqid"):
        qlen = int(hsps.qlen.iloc[0])
        coverage, identity, contig, start, end, strand = best_chain(hsps, qlen)
        row = dict(entries[mibig_id], contig=contig, start=start, end=end,
                   length=end - start + 1, identity=round(identity, 2),
                   coverage=round(coverage, 3), mibig_length=qlen)
        row["source_match"] = "own_record" if row["source_accession"].split(".")[0] in own else "other_record"
        if coverage >= MIN_COVERAGE and identity >= MIN_IDENTITY:
            accepted.append(row)
        elif coverage >= 0.5:
            near_misses.append(row)

    truth = pd.DataFrame(collapse(accepted, own))
    columns = ["mibig_id", "compound", "class", "contig", "start", "end", "length",
               "identity", "coverage", "mibig_length", "source_accession", "source_match",
               "quality", "completeness", "organism", "aliases", "nested_in"]
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    truth[columns].to_csv(args.out, sep="\t", index=False)
    if near_misses:
        pd.DataFrame(near_misses).sort_values("coverage", ascending=False).to_csv(
            args.out.replace(".tsv", ".near_misses.tsv"), sep="\t", index=False)
    print("{0}: {1} MIBiG loci searched, {2} accepted, {3} truth loci after collapsing, "
          "{4} near misses".format(os.path.basename(args.genome), n, len(accepted),
                                   len(truth), len(near_misses)))


if __name__ == "__main__":
    main()
