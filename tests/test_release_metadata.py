import unittest
from pathlib import Path

import yaml


class ReleaseMetadataTests(unittest.TestCase):
    def test_release_workflow_uses_limited_permissions_and_amd64(self):
        workflow = yaml.safe_load(Path(".github/workflows/release-container.yml").read_text())

        self.assertEqual(workflow["permissions"], {"contents": "read", "packages": "write"})
        job = workflow["jobs"]["publish"]
        rendered = yaml.safe_dump(job)
        self.assertIn("linux/amd64", rendered)
        self.assertIn("ghcr.io", rendered)
        self.assertIn("docker/build-push-action@v6", rendered)

    def test_ci_has_read_only_permissions(self):
        workflow = yaml.safe_load(Path(".github/workflows/ci.yml").read_text())

        self.assertEqual(workflow["permissions"], {"contents": "read"})


if __name__ == "__main__":
    unittest.main()
