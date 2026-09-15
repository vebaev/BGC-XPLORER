import unittest

from scripts.mibig_hits import extract_structured_mibig_hits, mibig_method_priority


class StructuredMibigHitTests(unittest.TestCase):
    def test_knownclusterblast_keeps_core_supported_top_hit_and_metric(self):
        payload = {"records": [{"modules": {
            "antismash.modules.clusterblast": {"knowncluster": {"results": [
                {"region_number": 2, "ranking": [
                    [{"accession": "BGC0000001.1", "description": "unrelated transporter"},
                     {"hits": 1, "core_gene_hits": 0, "similarity": 20}],
                    [{"accession": "BGC0000002.1", "description": "known product",
                      "cluster_type": "NRPS"},
                     {"hits": 5, "core_gene_hits": 2, "similarity": 17}],
                ]},
            ]}},
        }}]}

        hits = extract_structured_mibig_hits(payload)

        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0]["region_id"], "2")
        self.assertEqual(hits[0]["best_mibig_id"], "BGC0000002.1")
        self.assertEqual(hits[0]["matched_genes"], 5)
        self.assertEqual(hits[0]["core_gene_hits"], 2)
        self.assertEqual(hits[0]["match_score"], 17)
        self.assertEqual(hits[0]["score_metric"], "KnownClusterBlast empirical score")
        self.assertIsNone(hits[0]["mibig_similarity"])

    def test_clustercompare_score_is_not_percent_identity(self):
        payload = {"records": [{"modules": {
            "antismash.modules.cluster_compare": {"db_results": {"MIBiG": {
                "by_region": {"3": {"RegionToRegion_RiQ": {
                    "scores_by_region": {"BGC0000003.1: 0-1000": 0.72,
                                         "BGC0000004.1: 0-1000": 0.81},
                }}}
            }}},
        }}]}

        hits = extract_structured_mibig_hits(payload)

        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0]["best_mibig_id"], "BGC0000004.1")
        self.assertEqual(hits[0]["match_score"], 0.81)
        self.assertEqual(hits[0]["score_metric"], "ClusterCompare score (0-1)")
        self.assertIsNone(hits[0]["mibig_similarity"])

    def test_strong_peptide_match_is_retained_over_generic_cluster_comparison(self):
        self.assertLess(
            mibig_method_priority("comparippson_html", 81.8),
            mibig_method_priority("knownclusterblast", None),
        )
        self.assertGreater(
            mibig_method_priority("comparippson_html", 14.9),
            mibig_method_priority("knownclusterblast", None),
        )


if __name__ == "__main__":
    unittest.main()
