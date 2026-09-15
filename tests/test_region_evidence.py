import sys
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from aggregate_region_evidence import aggregate_evidence
from report_evidence_table import render_evidence_table


class RegionEvidenceTests(unittest.TestCase):
    def test_art_records_are_separate_deduplicated_and_contig_specific(self):
        consensus = pd.DataFrame([
            {"consensus_id": "sample_consensus_1", "consensus_label": "sample_1", "contig": "a",
             "start": 100, "end": 200, "support_tools": "antismash,gecco", "support_count": 2,
             "distance_to_left_edge": 99, "distance_to_right_edge": 800},
            {"consensus_id": "sample_consensus_2", "consensus_label": "sample_2", "contig": "b",
             "start": 100, "end": 200, "support_tools": "deepbgc", "support_count": 1,
             "distance_to_left_edge": 99, "distance_to_right_edge": 800},
        ])
        arts = pd.DataFrame([
            {"contig": "a", "start": 150, "end": 160, "feature": "resistance", "evidence": "known_hit"},
            {"contig": "a", "start": 150, "end": 160, "feature": "resistance", "evidence": "known_hit"},
            {"contig": "a", "start": 170, "end": 180, "feature": "PF00000", "evidence": "duf_hit"},
            {"contig": "b", "start": 150, "end": 160, "feature": "other", "evidence": "known_hit"},
        ])
        mibig = pd.DataFrame([
            {"consensus_id": "sample_consensus_1", "best_mibig_id": "BGC0000001.1",
             "match_score": 75, "score_metric": "KnownClusterBlast empirical score"}
        ])
        result = aggregate_evidence(consensus, arts, mibig).set_index("consensus_id")
        first = result.loc["sample_consensus_1"]
        second = result.loc["sample_consensus_2"]
        self.assertEqual((first.arts_known_hits, first.arts_duf_hits, first.arts_hits), (1, 1, 2))
        self.assertEqual((second.arts_known_hits, second.arts_duf_hits), (1, 0))
        self.assertEqual(first.best_mibig_id, "BGC0000001.1")
        self.assertEqual(first.nearest_contig_edge_bp, 99)
        self.assertNotIn("priority_score", result.columns)
        self.assertNotIn("confidence_category", result.columns)
        self.assertNotIn("interest_category", result.columns)

    def test_single_table_includes_all_rows_sortable_markers_and_safe_text(self):
        frame = pd.DataFrame([
            {"consensus_id": "s_consensus_1", "consensus_label": "s_1", "contig": "a",
             "start": 100, "end": 200, "length_bp": 101, "support_tools": "antismash,gecco",
             "support_count": 2, "arts_known_hits": 1, "arts_duf_hits": 0,
             "best_mibig_id": "BGC0000001.1", "evidence_source": "knownclusterblast",
             "match_score": 75, "score_metric": "KnownClusterBlast empirical score",
             "bgc_types": "terpene", "products": "<script>alert(1)</script>",
             "nearest_contig_edge_bp": 99},
            {"consensus_id": "s_consensus_2", "consensus_label": "s_2", "contig": "b",
             "start": 300, "end": 400, "length_bp": 101, "support_tools": "deepbgc",
             "support_count": 1, "arts_known_hits": 0, "arts_duf_hits": 2,
             "best_mibig_id": "", "bgc_types": "unknown", "products": "unknown",
             "nearest_contig_edge_bp": 299},
        ])
        html = render_evidence_table(
            frame, {}, lambda row, index: "<button>Map</button>", overlap_fraction=0.45
        )
        self.assertEqual(html.count("<tr data-search="), 2)
        self.assertIn("data-column='7' data-type='number'", html)
        self.assertIn("id='evidence-mibig'", html)
        self.assertIn("id='evidence-known'", html)
        self.assertIn("id='evidence-duf'", html)
        self.assertIn("KnownClusterBlast empirical score", html)
        self.assertIn("antiSMASH", html)
        self.assertIn("GECCO", html)
        self.assertIn("DeepBGC", html)
        self.assertIn("reciprocal overlap of at least 0.45", html)
        self.assertNotIn("<script>alert(1)</script>", html)
        self.assertNotIn("High-confidence", html)
        self.assertNotIn("Potentially Novel", html)


if __name__ == "__main__":
    unittest.main()
