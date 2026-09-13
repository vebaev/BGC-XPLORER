import os
import subprocess
import tempfile
import unittest
from pathlib import Path


class FetchBaktaDatabaseTests(unittest.TestCase):
    def test_bakta_helpers_are_available_to_downloader(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bin_dir = root / "bin"
            bin_dir.mkdir()
            (bin_dir / "amrfinder").write_text("#!/bin/sh\nexit 0\n")
            (bin_dir / "bakta_db").write_text(
                "#!/bin/sh\n"
                "command -v amrfinder >/dev/null || exit 42\n"
                "mkdir -p \"$3/db-light/amrfinderplus-db\"\n"
                "printf '{\"type\": true, \"major\": 6, \"minor\": 0}\\n' > \"$3/db-light/version.json\"\n"
                "touch \"$3/db-light/bakta.db\" \"$3/db-light/pscc.dmnd\"\n"
                "touch \"$3/db-light/amrfinderplus-db/database_format_version.txt\"\n"
            )
            for executable in bin_dir.iterdir():
                executable.chmod(0o755)

            environment = os.environ.copy()
            environment.update({
                "BGC_DB_ROOT": str(root / "database"),
                "BAKTA_DB_TYPE": "light",
                "BAKTA_ENV_BIN": str(bin_dir),
                "BAKTA_DB_BIN": str(bin_dir / "bakta_db"),
                "BGC_PYTHON": os.environ.get("PYTHON", "python3"),
            })
            result = subprocess.run(
                ["bash", "scripts/fetch_bakta_db.sh"],
                cwd=Path(__file__).resolve().parents[1],
                env=environment,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                universal_newlines=True,
            )

            self.assertEqual(result.returncode, 0, result.stdout)


if __name__ == "__main__":
    unittest.main()
