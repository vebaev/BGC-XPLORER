import unittest

from scripts.check_public_tree import find_public_tree_violations


class PublicTreeTests(unittest.TestCase):
    def test_rejects_private_and_generated_paths(self):
        paths = [
            ".env",
            "config/local_ai.env",
            "db/eggnog/eggnog.db",
            "results/sample/report.html",
            "tools/arts/.git/config",
            "scripts/__pycache__/common.pyc",
        ]

        self.assertEqual(find_public_tree_violations(paths), paths)

    def test_accepts_public_configuration_and_source(self):
        paths = [
            ".env.example",
            "db/manifest.yaml",
            "scripts/check_databases.py",
            "tools/arts/artspipeline1.py",
        ]

        self.assertEqual(find_public_tree_violations(paths), [])


if __name__ == "__main__":
    unittest.main()
