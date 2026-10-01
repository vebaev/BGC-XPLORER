import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from render_cluster_maps import annotation_conflict, render_svg, tooltip_for


CORE_GT1 = {"locus_tag": "G1", "predictor_core_evidence": "gecco: PF00501 (AMP-binding enzyme)",
            "dbcan_diamond": "GT1", "dbcan_hmm": "", "dbcan_subfamily": "", "dbcan_recommendation": ""}


class AnnotationConflictTests(unittest.TestCase):
    def test_cazyme_family_on_a_core_gene_is_a_conflict_naming_its_source(self):
        self.assertEqual(
            annotation_conflict(CORE_GT1),
            "dbCAN GT1 (DIAMOND) on a biosynthetic core gene (gecco: PF00501 (AMP-binding enzyme))",
        )

    def test_no_conflict_without_both_annotations(self):
        self.assertEqual(annotation_conflict({**CORE_GT1, "dbcan_diamond": ""}), "")
        self.assertEqual(annotation_conflict({**CORE_GT1, "predictor_core_evidence": ""}), "")

    def test_conflict_reaches_the_tooltip_and_the_gene_map(self):
        row = {**CORE_GT1, "contig": "c", "gene_start": 1, "gene_end": 900, "strand": "+",
               "bakta_product": "AMP-binding protein", "gene_category": "biosynthetic_core",
               "display_label": "G1"}
        row["annotation_conflict"] = annotation_conflict(row)
        row["tooltip_text"] = tooltip_for(row)
        self.assertIn("Annotation conflict: dbCAN GT1", row["tooltip_text"])
        cluster = {"consensus_id": "s_consensus_1", "cluster_label": "s_1", "contig": "c",
                   "start": 1, "end": 900, "union_start": 1, "union_end": 900}
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "map.svg"
            render_svg(cluster, pd.DataFrame([row]), out)
            svg = out.read_text()
        self.assertIn("<g class='gene conflict'>", svg)
        self.assertIn("Annotation conflict</text>", svg)


if __name__ == "__main__":
    unittest.main()
