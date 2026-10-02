import unittest
from pathlib import Path

RULE = (Path(__file__).resolve().parents[1] / "rules" / "arts.smk").read_text()


class ArtsRule(unittest.TestCase):

    def test_arts_helpers_are_on_path(self):
        # ARTS runs mafft and FastTree by name for its phylogeny check; they are installed next to the
        # ARTS Python, so that directory must come first on PATH when ARTS is called.
        self.assertIn('PATH="$(dirname "$PYTHON_BIN"):', RULE)


if __name__ == "__main__":
    unittest.main()
