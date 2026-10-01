import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from ai_cluster_server import trace_identifiers

PAYLOAD = {"genes": [{"locus_tag": "BFJFPI_01155"}, {"locus_tag": "BFJFPI_01157"}],
           "mibig": {"best_mibig_id": "BGC0002358.3"}}


class TraceabilityTests(unittest.TestCase):
    def test_identifiers_absent_from_the_payload_are_flagged(self):
        analysis = {
            "summary": "BFJFPI_01155 encodes the adenylation domain. It resembles BGC0002358. "
                       "BFJFPI_09999 adds a methyl group.",
            "key_genes": ["BFJFPI_01157 (NRPS): elongation", "BFJFPI_00001 (P450): oxidation"],
            "recommended_followup": ["Knock out the cluster."],
        }
        trace = trace_identifiers(analysis, PAYLOAD)
        self.assertEqual(trace["statements"], 6)
        self.assertEqual(trace["with_identifiers"], 5)
        self.assertEqual([item["missing"] for item in trace["untraceable"]], [["BFJFPI_09999"], ["BFJFPI_00001"]])
        self.assertEqual(trace["untraceable"][0]["field"], "summary")

    def test_a_versioned_mibig_id_matches_the_unversioned_payload_entry(self):
        trace = trace_identifiers({"summary": "Similar to BGC0002358.3."}, {"mibig": "BGC0002358"})
        self.assertEqual(trace["untraceable"], [])


if __name__ == "__main__":
    unittest.main()
