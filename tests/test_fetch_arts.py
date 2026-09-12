import os
import subprocess
import tempfile
import unittest
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class FetchArtsTests(unittest.TestCase):
    def make_archives(self, source):
        with zipfile.ZipFile(source / "actinobacteria.zip", "w") as archive:
            archive.writestr("actinobacteria/coremodels.hmm", "core")
            archive.writestr("actinobacteria/model_metadata.json", "{}")
            archive.writestr("actinobacteria/genematrix.txt", "genes")
        with zipfile.ZipFile(source / "hmm_models.zip", "w") as archive:
            archive.writestr("knownresistance.hmm", "known")
            archive.writestr("dufmodels.hmm", "duf")

    def run_fetch(self, database_root, source):
        environment = os.environ.copy()
        environment.update({
            "ARTS_REFERENCE": "actinobacteria",
            "ARTS_SOURCE_DIR": str(source),
            "BGC_DB_ROOT": str(database_root),
        })
        return subprocess.run(
            ["bash", str(ROOT / "scripts" / "fetch_arts.sh")],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            env=environment,
        )

    def test_installs_actinobacteria_and_skips_valid_second_run(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source"
            database_root = root / "db"
            source.mkdir()
            self.make_archives(source)

            first = self.run_fetch(database_root, source)
            second = self.run_fetch(database_root, source)

            self.assertTrue((database_root / "arts/actinobacteria/coremodels.hmm").is_file())
            self.assertTrue((database_root / "arts/knownresistance.hmm").is_file())
            self.assertTrue((database_root / "arts/dufmodels.hmm").is_file())
            self.assertTrue((database_root / "arts/.bgc-xplorer-actinobacteria.ready").is_file())
            self.assertIn("installed successfully", first.stdout)
            self.assertIn("already valid", second.stdout)

    def test_rejects_unknown_reference(self):
        with tempfile.TemporaryDirectory() as directory:
            environment = os.environ.copy()
            environment.update({
                "ARTS_REFERENCE": "unknown",
                "BGC_DB_ROOT": directory,
            })
            result = subprocess.run(
                ["bash", str(ROOT / "scripts" / "fetch_arts.sh")],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                env=environment,
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Unsupported ARTS reference", result.stderr)


if __name__ == "__main__":
    unittest.main()
