import tempfile
import unittest
from pathlib import Path

from scripts.startup_checks import missing_executables


class StartupChecksTests(unittest.TestCase):
    def test_accepts_executable_paths(self):
        with tempfile.TemporaryDirectory() as directory:
            executable = Path(directory) / "tool"
            executable.write_text("#!/bin/sh\n")
            executable.chmod(0o755)

            self.assertEqual(missing_executables([str(executable)]), [])

    def test_reports_missing_and_non_executable_paths(self):
        with tempfile.TemporaryDirectory() as directory:
            non_executable = Path(directory) / "tool"
            non_executable.write_text("data")

            missing = missing_executables([str(non_executable), "/missing/tool"])

            self.assertEqual(missing, [str(non_executable), "/missing/tool"])


if __name__ == "__main__":
    unittest.main()
