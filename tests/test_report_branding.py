import tempfile
import unittest
from pathlib import Path

from scripts.report_branding import image_data_uri, report_home_link
from scripts.report_design import DONUT_COLORS, PRIMARY_GLANCE_LABELS, reproducibility_panel


class ReportBrandingTests(unittest.TestCase):
    def test_at_a_glance_uses_four_primary_metrics(self):
        self.assertEqual(
            PRIMARY_GLANCE_LABELS,
            ("Consensus", "Multi-tool", "High-confidence", "High-interest"),
        )

    def test_reproducibility_panel_uses_live_model_placeholder_without_download_link(self):
        panel = reproducibility_panel({
            "application": {"version": "1.0", "git_commit": "abc"},
            "arts_reference": "actinobacteria",
            "ai": {"model": "nvidia/model-at-generation"},
        })

        self.assertIn("id='active-ai-model'", panel)
        self.assertIn("data-generation-model='nvidia/model-at-generation'", panel)
        self.assertIn("Checking active model", panel)
        self.assertNotIn("Download provenance", panel)
        self.assertEqual(panel.count("<p"), 1)
        self.assertIn("ARTS reference: <strong>actinobacteria</strong> · ", panel)
        self.assertIn("<span id='ai-model-label'>Active AI model:</span>", panel)

    def test_donut_palette_is_pastel(self):
        self.assertEqual(
            DONUT_COLORS,
            ("#b8a9e8", "#a8d8c7", "#a9cce8", "#f3c6a8", "#e8b4c4", "#9fd6d2"),
        )

    def test_embeds_first_existing_image_with_correct_media_type(self):
        with tempfile.TemporaryDirectory() as directory:
            logo = Path(directory) / "logo.jpg"
            logo.write_bytes(b"jpeg-bytes")

            uri = image_data_uri([Path(directory) / "missing.png", logo])

            self.assertEqual(uri, "data:image/jpeg;base64,anBlZy1ieXRlcw==")

    def test_home_link_returns_to_upload_page_from_iframe_or_new_tab(self):
        link = report_home_link()

        self.assertIn("href='/'", link)
        self.assertIn("target='_top'", link)
        self.assertIn("Back to BGC-XPLORER", link)


if __name__ == "__main__":
    unittest.main()
