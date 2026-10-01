import runpy
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))


class ParseDeepBGCCoordinates(unittest.TestCase):
    """DeepBGC is 0-based and end-exclusive; the parsed table must be 1-based and end-inclusive."""

    def test_starts_shift_by_one_and_ends_do_not(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            # One cluster over a gene that Prodigal reports as 1144..3813 (1-based, inclusive).
            pd.DataFrame([{
                "sequence_id": "contig_1", "nucl_start": 1143, "nucl_end": 3813,
                "bgc_candidate_id": "contig_1_1143-3813.1", "product_class": "NRP",
                "product_activity": "antibacterial", "detector_score": 0.9, "classifier_score": 0.8,
                "bio_pfam_ids": "PF00109",
            }]).to_csv(tmp / "deepbgc.bgc.tsv", sep="\t", index=False)
            pd.DataFrame([{
                "sequence_id": "contig_1", "protein_id": "p1", "gene_start": 1143, "gene_end": 3813,
                "gene_strand": 1, "pfam_id": "PF00109", "deepbgc_score": 0.9, "in_cluster": 1,
            }]).to_csv(tmp / "deepbgc.pfam.tsv", sep="\t", index=False)
            (tmp / "done").write_text("")
            snakemake = SimpleNamespace(
                wildcards=SimpleNamespace(sample="s"),
                input=SimpleNamespace(done=str(tmp / "done")),
                output=[str(tmp / "out.tsv")],
            )
            runpy.run_path(str(ROOT / "scripts" / "parse_deepbgc.py"),
                           init_globals={"snakemake": snakemake}, run_name="__main__")
            row = pd.read_csv(tmp / "out.tsv", sep="\t").iloc[0]
            self.assertEqual((row["start"], row["end"]), (1144, 3813))
            self.assertIn("1144", str(row["core_gene_records"]))
            self.assertNotIn("1143", str(row["core_gene_records"]))


if __name__ == "__main__":
    unittest.main()
