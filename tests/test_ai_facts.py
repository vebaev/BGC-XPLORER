import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from ai_facts import build_fact_card, result, template_summary, verify  # noqa: E402

PAYLOAD = {
    "region": {"contig": "contig_1", "start": "1000", "end": "5000", "length_bp": "4001",
               "support_tools": "antismash,gecco", "bgc_types": "NRPS,terpene",
               "core_gene_support": "antismash:ABC_00002; gecco:ABC_00002,ABC_00003",
               "arts_known_hits": "0", "arts_duf_hits": "2", "nearest_contig_edge_bp": "999"},
    "genes": [{"locus_tag": "ABC_00001", "gene_start": "500", "gene_end": "900", "bakta_product": "Transporter"},
              {"locus_tag": "ABC_00002", "gene_start": "1000", "gene_end": "3000", "bakta_product": "NRPS"},
              {"locus_tag": "ABC_00003", "gene_start": "3100", "gene_end": "4900", "bakta_product": "Thioesterase"}],
    "mibig_dereplication": {"best_mibig_id": "BGC0000001.1", "best_mibig_product": "examplomycin/examplomycin B",
                            "score_metric": "KnownClusterBlast empirical score", "match_score": "100.0",
                            "matched_genes": "12.0", "core_gene_hits": "2.0"},
}
FACTS = build_fact_card(PAYLOAD, conflicts=["ABC_00002"])
GOOD = ("The locus spans contig_1:1000-5000 (4001 bp, 2 genes) and is supported by antiSMASH and GECCO "
        "(2 of 3 callers); its MIBiG compound class is NRP and Terpene. ABC_00002 (NRPS) is a core gene named by "
        "antiSMASH and GECCO, and ABC_00003 (Thioesterase) is a core gene named by GECCO. The closest MIBiG entry is "
        "BGC0000001.1 (examplomycin), KnownClusterBlast empirical score 100 on 12 matched genes, 2 of them core genes. "
        "ARTS reports 0 known resistance hits and 2 DUF hits. 1 gene carries an annotation conflict: ABC_00002. "
        "The locus lies 999 bp from a contig edge.")


class FactCard(unittest.TestCase):

    def test_card_contents(self):
        keys = [f["key"] for f in FACTS]
        self.assertEqual(keys, ["location", "callers", "class", "core:ABC_00002", "core:ABC_00003",
                                "mibig", "arts", "conflict", "edge"])
        self.assertIn("2 genes", FACTS[0]["sentence"])  # genes with their midpoint in the span

    def test_template_always_verifies(self):
        self.assertTrue(verify(template_summary(FACTS), FACTS)["verified"])

    def test_faithful_paragraph_verifies(self):
        check = verify(GOOD, FACTS)
        self.assertTrue(check["verified"], check["reasons"])
        self.assertEqual(check["facts_kept"], len(FACTS))

    def test_dropped_fact_fails(self):
        check = verify(GOOD.replace(" ARTS reports 0 known resistance hits and 2 DUF hits.", ""), FACTS)
        self.assertFalse(check["verified"])
        self.assertTrue(any("arts" in r for r in check["reasons"]))

    def test_added_number_fails(self):
        self.assertFalse(verify(GOOD + " It holds 17 genes.", FACTS)["verified"])

    def test_added_identifier_fails(self):
        self.assertFalse(verify(GOOD + " ABC_00009 is nearby.", FACTS)["verified"])

    def test_caller_not_in_facts_fails(self):
        self.assertFalse(verify(GOOD.replace("named by GECCO.", "named by GECCO and DeepBGC."), FACTS)["verified"])

    def test_core_gene_must_share_a_sentence_with_its_callers(self):
        moved = GOOD.replace("ABC_00003 (Thioesterase) is a core gene named by GECCO.",
                             "ABC_00003 (Thioesterase) is a core gene. GECCO names it.")
        self.assertFalse(verify(moved, FACTS)["verified"])

    def test_grouped_genes_share_their_callers(self):
        facts = build_fact_card(dict(PAYLOAD, region=dict(PAYLOAD["region"], core_gene_support="gecco:ABC_00002,ABC_00003")))
        text = template_summary([f for f in facts if not f["key"].startswith("core:")]) + (
            " ABC_00002 (NRPS) and ABC_00003 (Thioesterase) are core genes named by GECCO.")
        self.assertTrue(verify(text, facts)["verified"], verify(text, facts)["reasons"])
        wrong = text.replace("named by GECCO.", "named by GECCO and antiSMASH.")
        self.assertFalse(verify(wrong, facts)["verified"])

    def test_numbers_are_not_matched_inside_locus_tags(self):
        dropped = GOOD.replace(" 1 gene carries an annotation conflict: ABC_00002.", " A gene carries an annotation conflict.")
        self.assertFalse(verify(dropped, FACTS)["verified"])

    def test_number_words_count_as_numbers(self):
        self.assertTrue(verify(GOOD.replace("1 gene carries", "One gene carries"), FACTS)["verified"])
        self.assertFalse(verify(GOOD.replace("1 gene carries", "Two genes carry"), FACTS)["verified"])

    def test_interpretation_fails(self):
        self.assertFalse(verify(GOOD + " The locus likely produces a peptide.", FACTS)["verified"])

    def test_failed_text_is_replaced_by_template(self):
        shown = result({"summary": GOOD + " It probably makes an antibiotic."}, FACTS)
        self.assertEqual(shown["shown_source"], "template")
        self.assertEqual(shown["shown_summary"], template_summary(FACTS))
        self.assertEqual(result({"summary": GOOD}, FACTS)["shown_source"], "model")


if __name__ == "__main__":
    unittest.main()
