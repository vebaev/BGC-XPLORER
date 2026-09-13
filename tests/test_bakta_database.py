import json
import tempfile
import unittest
from pathlib import Path

from scripts.check_bakta_database import inspect_bakta_database


def write_database(root, db_type):
    name = "db-light" if db_type == "light" else "db"
    path = root / "bakta" / name
    path.mkdir(parents=True)
    (path / "version.json").write_text(json.dumps({"major": 6, "minor": 0, "date": "2025-02-24", "type": db_type == "light"}))
    (path / "bakta.db").write_text("database")
    (path / ("pscc.dmnd" if db_type == "light" else "psc.dmnd")).write_text("diamond")
    (path / "amrfinderplus-db").mkdir()
    (path / "amrfinderplus-db" / "database_format_version.txt").write_text("1")
    return path


class BaktaDatabaseTests(unittest.TestCase):
    def test_accepts_selected_light_database(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = write_database(root, "light")
            status = inspect_bakta_database(root, "light")
            self.assertTrue(status["ready"])
            self.assertEqual(status["path"], str(path))

    def test_full_selection_does_not_reuse_light_database(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_database(root, "light")
            status = inspect_bakta_database(root, "full")
            self.assertFalse(status["ready"])
            self.assertTrue(status["path"].endswith("/bakta/db"))

    def test_rejects_database_without_amrfinder_setup(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = write_database(root, "full")
            for child in (path / "amrfinderplus-db").iterdir():
                child.unlink()
            status = inspect_bakta_database(root, "full")
            self.assertFalse(status["ready"])
            self.assertIn("amrfinderplus-db", status["missing"])


if __name__ == "__main__":
    unittest.main()
