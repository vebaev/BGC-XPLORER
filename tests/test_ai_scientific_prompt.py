import unittest
from types import SimpleNamespace
from scripts.ai_cluster_server import (ANALYSIS_JSON_SCHEMA, compact_gene_for_model, normalize_analysis, analysis_fingerprint)

class ScientificPromptTests(unittest.TestCase):
    def test_preserves_evidence(self):
        gene = {'locus_tag': 'test1', 'gene_start': 10, 'strand': '+', 'ec': '1.2.3.4', 'kegg_ko': 'K00001', 'bakta_product': 'long annotation ' * 30}
        result = compact_gene_for_model(gene)
        for key, value in gene.items():
            self.assertEqual(result[key], str(value).strip())

    def test_six_sections(self):
        result = normalize_analysis({'summary': 'Observed', 'confidence': 'high', 'caveats': 'x', 'novelty_assessment': 'x'})
        self.assertEqual(set(result), set(ANALYSIS_JSON_SCHEMA['required']))
        self.assertEqual(len(result), 6)

    def test_cache_changes_with_model_evidence_and_parameters(self):
        server = SimpleNamespace(model='a', base_url='https://example.test', temperature=.1, top_p=.95, max_tokens=1800, reasoning_effort='none', seed=42, use_guided_json=True)
        initial = analysis_fingerprint({'genes': []}, server)
        self.assertNotEqual(initial, analysis_fingerprint({'genes': ['new']}, server))
        server.model = 'b'
        self.assertNotEqual(initial, analysis_fingerprint({'genes': []}, server))
