import runpy
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

BGC_COLUMNS = ["sample", "tool", "contig", "start", "end", "strand", "bgc_id", "bgc_type",
               "product", "score", "confidence", "core_gene_records", "source_file"]
# Six 100 bp genes, 100 bp apart.
GENES = [(1 + 200 * i, 100 + 200 * i) for i in range(6)]


def prediction(tool, start, end):
    return {"sample": "s", "tool": tool, "contig": "contig_1", "start": start, "end": end,
            "strand": ".", "bgc_id": tool + "_1", "bgc_type": "NRPS", "product": "NRPS",
            "score": "", "confidence": "", "core_gene_records": "", "source_file": ""}


def run_consensus(tmp, predictions, **consensus):
    tmp = Path(tmp)
    by_tool = {tool: [p for p in predictions if p["tool"] == tool] for tool in ("antismash", "gecco", "deepbgc")}
    for tool, rows in by_tool.items():
        pd.DataFrame(rows, columns=BGC_COLUMNS).to_csv(tmp / f"{tool}.tsv", sep="\t", index=False)
    pd.DataFrame(columns=["sample", "tool", "contig", "start", "end", "feature", "score", "evidence",
                          "source_file"]).to_csv(tmp / "arts.tsv", sep="\t", index=False)
    pd.DataFrame([
        {"contig": "contig_1", "type": "cds", "start": start, "end": end, "strand": "+",
         "locus_tag": f"G{number}", "gene": "", "product": "", "dbxrefs": ""}
        for number, (start, end) in enumerate(GENES, 1)
    ]).to_csv(tmp / "bakta.tsv", sep="\t", index=False)
    (tmp / "genome.fna").write_text(">contig_1\n" + "A" * 1200 + "\n")
    snakemake = SimpleNamespace(
        input=SimpleNamespace(antismash=str(tmp / "antismash.tsv"), gecco=str(tmp / "gecco.tsv"),
                              deepbgc=str(tmp / "deepbgc.tsv"), arts=str(tmp / "arts.tsv"),
                              bakta=str(tmp / "bakta.tsv"), fna=str(tmp / "genome.fna")),
        output=SimpleNamespace(consensus=str(tmp / "consensus.tsv"), overlap=str(tmp / "overlap.tsv"),
                               gene_support=str(tmp / "support.tsv")),
        config={"consensus": {"mode": "voting", "min_callers_per_gene": 2, **consensus}},
    )
    runpy.run_path(str(ROOT / "scripts" / "build_consensus.py"),
                   init_globals={"snakemake": snakemake}, run_name="__main__")
    consensus_table = pd.read_csv(tmp / "consensus.tsv", sep="\t")
    support = pd.read_csv(tmp / "support.tsv", sep="\t")
    return consensus_table, dict(zip(support.locus_tag, support.callers))


# antiSMASH covers all six genes; GECCO covers 51 % of gene 1 and genes 2-4 in full.
PREDICTIONS = [prediction("antismash", 1, 1100), prediction("gecco", 50, 700)]


def multi_caller(table):
    rows = table[table.support_count >= 2]
    return [(int(row.start), int(row.end)) for _, row in rows.iterrows()]


class VotingRuleTests(unittest.TestCase):
    def test_defaults_count_any_overlap_and_accept_any_run_length(self):
        with tempfile.TemporaryDirectory() as tmp:
            table, support = run_consensus(tmp, PREDICTIONS)
        self.assertEqual(multi_caller(table), [(1, 700)])
        self.assertEqual(support["G1"], 2)

    def test_overlap_fraction_drops_a_partial_vote(self):
        with tempfile.TemporaryDirectory() as tmp:
            table, support = run_consensus(tmp, PREDICTIONS, min_gene_overlap_fraction=0.6)
        self.assertEqual(multi_caller(table), [(201, 700)])
        self.assertEqual(support["G1"], 1)

    def test_runs_shorter_than_the_minimum_keep_predictions_whole(self):
        with tempfile.TemporaryDirectory() as tmp:
            table, _ = run_consensus(tmp, PREDICTIONS, min_genes_per_locus=5)
        self.assertEqual(multi_caller(table), [])
        self.assertEqual(sorted(zip(table.start, table.end)), [(1, 1100), (50, 700)])

    def test_a_gene_below_the_threshold_splits_the_run(self):
        predictions = [prediction("antismash", 1, 1100), prediction("gecco", 1, 300),
                       prediction("deepbgc", 601, 1100)]
        with tempfile.TemporaryDirectory() as tmp:
            table, _ = run_consensus(tmp, predictions)
        self.assertEqual(multi_caller(table), [(1, 300), (601, 1100)])


if __name__ == "__main__":
    unittest.main()
