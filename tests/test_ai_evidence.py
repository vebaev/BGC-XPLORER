import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from ai_evidence import check_statement, coverage, normalize, verify  # noqa: E402

PAYLOAD = {
    "region": {"contig": "contig_1", "start": "1332869", "end": "1390082", "length_bp": "57214",
               "support_tools": "antismash,deepbgc,gecco",
               "core_gene_support": "antismash:ABC_01147; gecco:ABC_01147,ABC_01154"},
    "gene_count": 36,
    "gene_category_counts": {"transporter": 11},
    "representative_genes": [
        {"locus_tag": "ABC_01147", "bakta_product": "D-alanine--poly(phosphoribitol) ligase subunit DltA",
         "pfams": "AMP-binding,AMP-binding_C,Condensation", "gene_category": "biosynthetic_core"},
        {"locus_tag": "ABC_01154", "bakta_product": "Thioesterase", "gene_category": "biosynthetic_core"},
    ],
    "tool_predictions": {"gecco": [{"start": "1336178", "end": "1389372"}]},
    "arts_hits": [{"feature": "PF09594.5", "evidence": "duf_hit"}],
    "mibig_dereplication": {"best_mibig_id": "BGC0002358.3", "best_mibig_product": "cyclofaulknamycin",
                            "match_score": "100.0", "matched_genes": "24.0"},
}


def statement(text, *paths):
    return {"text": text, "evidence": list(paths)}


class EvidenceCheck(unittest.TestCase):

    def test_restated_values_verify(self):
        ok, reasons = check_statement(statement(
            "ABC_01147 is annotated as D-alanine--poly(phosphoribitol) ligase subunit DltA with Pfam domains AMP-binding, AMP-binding_C and Condensation.",
            "ABC_01147.bakta_product",
            "ABC_01147.pfams"), PAYLOAD)
        self.assertTrue(ok, reasons)

    def test_paths_and_numbers_verify(self):
        ok, reasons = check_statement(statement(
            "KnownClusterBlast matches BGC0002358.3 (cyclofaulknamycin) with score 100.0 on 24.0 genes.",
            "mibig_dereplication.best_mibig_id",
            "mibig_dereplication.best_mibig_product",
            "mibig_dereplication.match_score",
            "mibig_dereplication.matched_genes"), PAYLOAD)
        self.assertTrue(ok, reasons)

    def test_list_index_path(self):
        ok, reasons = check_statement(statement(
            "GECCO reports contig_1 from 1,336,178 to 1,389,372.",
            "tool_predictions.gecco[0].start",
            "tool_predictions.gecco[0].end"), PAYLOAD)
        self.assertTrue(ok, reasons)

    def test_wrong_value_fails(self):
        ok, reasons = check_statement(statement(
            "The MIBiG match score is 80.", "mibig_dereplication.match_score"), PAYLOAD)
        self.assertFalse(ok)
        self.assertTrue(any("80" in r for r in reasons))

    def test_number_not_in_evidence_fails(self):
        ok, reasons = check_statement(statement(
            "The locus holds 37 genes.", "gene_count"), PAYLOAD)
        self.assertFalse(ok)
        self.assertTrue(any("37" in r for r in reasons))

    def test_inference_wording_fails(self):
        ok, reasons = check_statement(statement(
            "ABC_01154 Thioesterase likely releases the peptide.",
            "ABC_01154.bakta_product"), PAYLOAD)
        self.assertFalse(ok)
        self.assertIn("interpretive wording", reasons[0])

    def test_unknown_locus_tag_fails(self):
        ok, reasons = check_statement(statement(
            "ABC_09999 is a Thioesterase.",
            "ABC_09999.bakta_product"), PAYLOAD)
        self.assertFalse(ok)

    def test_identifier_must_be_cited(self):
        ok, reasons = check_statement(statement(
            "ABC_01154 sits next to ABC_01147, a Thioesterase.",
            "ABC_01154.bakta_product"), PAYLOAD)
        self.assertFalse(ok)
        self.assertTrue(any("ABC_01147" in r for r in reasons))

    def test_cited_value_must_be_restated(self):
        ok, reasons = check_statement(statement(
            "ABC_01154 is a hydrolase.",
            "ABC_01154.bakta_product"), PAYLOAD)
        self.assertFalse(ok)
        self.assertTrue(any("does not restate" in r for r in reasons))

    def test_names_with_digits_are_not_numbers(self):
        ok, reasons = check_statement(statement(
            "ARTS reports PF09594.5 as a duf_hit on contig_1.",
            "arts_hits[0].feature",
            "arts_hits[0].evidence"), PAYLOAD)
        self.assertTrue(ok, reasons)

    def test_coverage_counts_callers_core_genes_mibig_arts(self):
        analysis = normalize({"overview": [statement("antiSMASH, DeepBGC and GECCO support the locus.",
                                                      "region.support_tools")],
                              "core_genes": [statement("ABC_01147 is a core gene.",
                                                        "ABC_01147.gene_category")]})
        result = coverage(analysis, PAYLOAD)
        self.assertEqual(result["callers"], {"expected": 3, "covered": 3})
        self.assertEqual(result["core_genes"], {"expected": 2, "covered": 1})
        self.assertEqual(result["mibig"]["covered"], 0)
        self.assertEqual(result["arts"]["covered"], 0)

    def test_core_support_callers_checked(self):
        ok, _ = check_statement(statement(
            "ABC_01154 is annotated as Thioesterase and is named as a core gene by gecco.",
            "ABC_01154.bakta_product", "region.core_gene_support"), PAYLOAD)
        self.assertTrue(ok)
        ok, reasons = check_statement(statement(
            "ABC_01147 is named as a core gene by antismash only.", "region.core_gene_support"), PAYLOAD)
        self.assertFalse(ok)
        self.assertIn("antismash, gecco", reasons[0])

    def test_core_support_negation_checked(self):
        ok, _ = check_statement(statement(
            "antismash does not name ABC_01154 as a core gene.", "region.core_gene_support"), PAYLOAD)
        self.assertTrue(ok)
        ok, reasons = check_statement(statement(
            "gecco does not name ABC_01154 as a core gene.", "region.core_gene_support"), PAYLOAD)
        self.assertFalse(ok)

    def test_empty_value_needs_negation(self):
        payload = dict(PAYLOAD, dbcan_cgc=[])
        self.assertTrue(check_statement(statement("No dbCAN CGC is reported.", "dbcan_cgc"), payload)[0])
        self.assertFalse(check_statement(statement("A dbCAN CGC is reported.", "dbcan_cgc"), payload)[0])

    def test_region_fields_found_without_prefix(self):
        payload = dict(PAYLOAD, region=dict(PAYLOAD["region"], arts_known_hits="0"))
        self.assertTrue(check_statement(statement("ARTS reports 0 known hits.", "arts_known_hits"), payload)[0])

    def test_whole_record_is_not_evidence(self):
        ok, reasons = check_statement(statement("There is a MIBiG comparison.", "mibig_dereplication"), PAYLOAD)
        self.assertFalse(ok)
        self.assertIn("whole records", reasons[0])
        ok, reasons = check_statement(statement("The MIBiG match score is 100.0.", "mibig_dereplication",
                                                "mibig_dereplication.match_score"), PAYLOAD)
        self.assertTrue(ok, reasons)

    def test_long_value_must_be_quoted_in_part(self):
        payload = dict(PAYLOAD, representative_genes=PAYLOAD["representative_genes"] + [
            {"locus_tag": "ABC_01143", "eggnog_description": "ABC-type antimicrobial peptide transport system, permease component"}])
        self.assertTrue(check_statement(statement(
            "ABC_01143 is described as an ABC-type antimicrobial peptide transport system permease.",
            "ABC_01143.eggnog_description"), payload)[0])
        self.assertFalse(check_statement(statement(
            "ABC_01143 exports the product.", "ABC_01143.eggnog_description"), payload)[0])

    def test_verify_totals(self):
        analysis = normalize({"overview": [statement("The locus holds 36 genes.", "gene_count"),
                                           statement("The locus holds 40 genes.", "gene_count")]})
        result = verify(analysis, PAYLOAD)
        self.assertEqual((result["statements"], result["verified"]), (2, 1))


if __name__ == "__main__":
    unittest.main()
