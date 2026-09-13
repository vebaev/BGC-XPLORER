import unittest
from pathlib import Path


class EggnogRerunTests(unittest.TestCase):
    def test_incomplete_outputs_are_overridden_in_both_execution_modes(self):
        rule = (Path(__file__).resolve().parents[1] / "rules" / "eggnog.smk").read_text()
        self.assertEqual(rule.count("--override"), 2)


if __name__ == "__main__":
    unittest.main()
