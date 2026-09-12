import tempfile
import unittest
from pathlib import Path

import yaml

from scripts.check_databases import inspect_databases, missing_required_names


class InspectDatabasesTests(unittest.TestCase):
    def test_maps_manifest_paths_under_external_database_root(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "antismash" / "pfam").mkdir(parents=True)
            (root / "antismash" / "pfam" / "Pfam-A.hmm").write_text("ready")
            manifest = {
                "antismash": {
                    "required": True,
                    "host_path": "db/antismash",
                    "required_files": ["pfam/Pfam-A.hmm"],
                }
            }

            status = inspect_databases(manifest, root)

            self.assertTrue(status["antismash"]["ready"])
            self.assertEqual(status["antismash"]["path"], str(root / "antismash"))

    def test_reports_only_missing_required_resources(self):
        with tempfile.TemporaryDirectory() as directory:
            manifest = {
                "required_db": {
                    "required": True,
                    "host_path": "db/required_db",
                    "required_files": ["sentinel.dat"],
                },
                "optional_db": {
                    "required": False,
                    "host_path": "db/optional_db",
                    "required_files": ["sentinel.dat"],
                },
            }

            status = inspect_databases(manifest, Path(directory))

            self.assertFalse(status["required_db"]["ready"])
            self.assertTrue(status["optional_db"]["ready"])
            self.assertEqual(missing_required_names(status), ["required_db"])


if __name__ == "__main__":
    unittest.main()
