import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from core_gene_evidence import (
    antismash_core_records,
    bakta_core_match,
    domain_core_record,
    parse_records,
)


class CoreGeneEvidenceTests(unittest.TestCase):
    def test_antismash_accepts_only_explicit_biosynthetic_roles(self):
        features = [
            {"type": "CDS", "location": "[99:300](+)", "qualifiers": {
                "locus_tag": ["core_1"], "gene_kind": ["biosynthetic"],
                "gene_functions": ["biosynthetic (rule-based-clusters) terpene: terpene_synth"],
            }},
            {"type": "CDS", "location": "[400:600](+)", "qualifiers": {
                "locus_tag": ["reg_1"], "gene_kind": ["regulatory"],
                "product": ["transcriptional regulator"],
            }},
        ]
        records = antismash_core_records(features)
        self.assertEqual(len(records), 1)
        parsed = parse_records(records[0])[0]
        self.assertEqual(parsed["identifier"], "core_1")
        self.assertEqual((parsed["start"], parsed["end"]), (100, 300))
        self.assertEqual(parsed["source"], "antismash")

    def test_domain_evidence_excludes_generic_tailoring_domains(self):
        record = domain_core_record("protein_1", 10, 90, "gecco", ["PF00501", "PF00005"])
        parsed = parse_records(record)[0]
        self.assertIn("PF00501", parsed["evidence"])
        self.assertNotIn("PF00005", parsed["evidence"])
        self.assertEqual(domain_core_record("protein_2", 10, 90, "gecco", ["PF00005"]), "")

    def test_bakta_fallback_is_strict(self):
        self.assertEqual(
            bakta_core_match({"product": "Type I polyketide synthase"}),
            "polyketide synthase",
        )
        self.assertEqual(bakta_core_match({"product": "putative oxidase"}), "")
        self.assertEqual(bakta_core_match({"product": "SAM-dependent methyltransferase"}), "")
        self.assertEqual(bakta_core_match({"product": "DNA repair cyclase-like protein"}), "")


if __name__ == "__main__":
    unittest.main()
