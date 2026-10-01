import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from render_cluster_maps import annotation_conflict, classify_gene, display_label, render_svg, tooltip_for


# A family recommended by dbCAN (two methods agree) on a core gene.
CORE_GT1 = {"locus_tag": "G1", "predictor_core_evidence": "gecco: PF00501 (AMP-binding enzyme)",
            "dbcan_diamond": "GT1", "dbcan_hmm": "GT1(10-300)", "dbcan_subfamily": "", "dbcan_recommendation": "GT1"}
# The same family from DIAMOND alone, which dbCAN does not recommend.
DIAMOND_ONLY = {**CORE_GT1, "dbcan_hmm": "", "dbcan_recommendation": ""}


class AnnotationConflictTests(unittest.TestCase):
    def test_cazyme_family_on_a_core_gene_is_a_conflict_naming_its_source(self):
        self.assertEqual(
            annotation_conflict(CORE_GT1),
            "dbCAN GT1 (HMMER, DIAMOND) on a biosynthetic core gene (gecco: PF00501 (AMP-binding enzyme))",
        )

    def test_single_method_hit_is_not_a_conflict_but_is_shown(self):
        self.assertEqual(annotation_conflict(DIAMOND_ONLY), "")
        self.assertIn("single method DIAMOND GT1; not recommended by dbCAN", tooltip_for(DIAMOND_ONLY))

    def test_single_method_hit_does_not_make_a_cazyme_or_a_label(self):
        plain = {**DIAMOND_ONLY, "predictor_core_evidence": "", "bakta_product": "Catalase"}
        self.assertNotEqual(classify_gene(plain), "cazyme")
        self.assertNotIn("GT1", display_label(plain))
        self.assertEqual(classify_gene({**plain, "dbcan_hmm": "GT1(1-9)", "dbcan_recommendation": "GT1"}), "cazyme")

    def test_no_conflict_without_both_annotations(self):
        self.assertEqual(annotation_conflict({**CORE_GT1, "dbcan_recommendation": ""}), "")
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
