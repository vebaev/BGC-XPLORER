import tempfile
import unittest
from pathlib import Path

from scripts.common import load_table_if_exists


class CommonTableTests(unittest.TestCase):
    def test_empty_existing_tsv_is_a_valid_empty_table(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "tool_overlap.tsv"
            path.write_text("\n", encoding="utf-8")

            table = load_table_if_exists(path, ["sample", "overlap_count"])

        self.assertTrue(table.empty)
        self.assertEqual(table.columns.tolist(), ["sample", "overlap_count"])


if __name__ == "__main__":
    unittest.main()
