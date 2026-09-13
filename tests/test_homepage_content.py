import unittest

from scripts.homepage_content import HERO_BODY, HERO_TITLE, RESULT_FEATURES, WORKFLOW_STEPS, section_header


class HomepageContentTests(unittest.TestCase):
    def test_section_header_keeps_accent_with_title_and_description_below(self):
        markup = section_header("Start a new analysis", "Description")

        self.assertIn("<div class='section-head'><h2>Start a new analysis</h2><span class='section-accent'></span></div>", markup)
        self.assertTrue(markup.endswith("<p class='section-description muted'>Description</p>"))

    def test_hero_explains_scientific_purpose(self):
        text = "{0} {1}".format(HERO_TITLE, HERO_BODY).lower()

        self.assertIn("biosynthetic gene clusters", text)
        self.assertIn("priorit", text)
        self.assertNotIn("upload page", text)

    def test_workflow_covers_discovery_context_and_integrated_report(self):
        self.assertEqual(len(WORKFLOW_STEPS), 5)
        rendered = " ".join(" ".join(step) for step in WORKFLOW_STEPS)

        for name in ("antiSMASH", "GECCO", "DeepBGC", "eggNOG", "dbCAN", "ARTS", "MIBiG"):
            self.assertIn(name, rendered)
        self.assertEqual(sum("Bakta" in " ".join(step) for step in WORKFLOW_STEPS), 1)

    def test_results_describe_the_primary_scientific_outputs(self):
        rendered = " ".join(" ".join(feature) for feature in RESULT_FEATURES).lower()

        for term in ("consensus", "priorit", "resistance", "gene maps", "provenance", "ai"):
            self.assertIn(term, rendered)


if __name__ == "__main__":
    unittest.main()
