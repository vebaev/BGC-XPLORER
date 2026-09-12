import json
import tempfile
import unittest
from pathlib import Path

import yaml

from scripts.build_provenance import build_provenance, sha256_file


class BuildProvenanceTests(unittest.TestCase):
    def test_records_reproducibility_inputs_without_secrets(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_path = root / "sample.fna"
            artifact_path = root / "consensus.tsv"
            config_path = root / "config.yaml"
            manifest_path = root / "manifest.yaml"
            database_file = root / "db/arts/actinobacteria/coremodels.hmm"
            input_path.write_text(">contig\nACGT\n")
            artifact_path.write_text("candidate\nBGC1\n")
            database_file.parent.mkdir(parents=True)
            database_file.write_text("model")
            config_path.write_text(yaml.safe_dump({
                "tools": {"arts": {"reference_set": "actinobacteria"}},
                "api_key": "config-secret-must-not-appear",
            }))
            manifest_path.write_text(yaml.safe_dump({
                "arts": {
                    "required": True,
                    "host_path": "db/arts",
                    "source": "https://example.test/arts",
                    "required_files": ["actinobacteria/coremodels.hmm"],
                }
            }))

            provenance = build_provenance(
                sample="sample",
                input_paths=[str(input_path)],
                artifact_paths=[str(artifact_path)],
                config_path=str(config_path),
                database_manifest_path=str(manifest_path),
                database_root=str(root / "db"),
                environment={
                    "BGC_XPLORER_VERSION": "0.1.0",
                    "BGC_XPLORER_COMMIT": "abc123",
                    "BGC_IMAGE_REFERENCE": "ghcr.io/example/bgc@sha256:123",
                    "NVIDIA_MODEL": "nvidia/test-model",
                    "NVIDIA_API_KEY": "must-not-appear",
                    "ARTS_REFERENCE": "actinobacteria",
                    "BGC_TOOL_VERSIONS_JSON": '{"antismash":"8.0.4"}',
                },
            )

            encoded = json.dumps(provenance)
            self.assertNotIn("must-not-appear", encoded)
            self.assertNotIn("config-secret-must-not-appear", encoded)
            self.assertEqual(provenance["configuration"]["effective"]["api_key"], "[REDACTED]")
            self.assertEqual(provenance["application"]["version"], "0.1.0")
            self.assertEqual(provenance["application"]["git_commit"], "abc123")
            self.assertEqual(provenance["ai"]["model"], "nvidia/test-model")
            self.assertEqual(provenance["arts_reference"], "actinobacteria")
            self.assertEqual(provenance["inputs"][0]["sha256"], sha256_file(str(input_path)))
            self.assertEqual(provenance["databases"]["arts"]["files"][0]["sha256"], sha256_file(str(database_file)))
            self.assertEqual(provenance["tools"]["antismash"], "8.0.4")

    def test_missing_optional_database_file_is_recorded(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config_path = root / "config.yaml"
            manifest_path = root / "manifest.yaml"
            config_path.write_text("tools: {}\n")
            manifest_path.write_text(yaml.safe_dump({
                "optional": {
                    "required": False,
                    "host_path": "db/optional",
                    "required_files": ["missing.dat"],
                }
            }))

            provenance = build_provenance(
                sample="sample",
                input_paths=[],
                artifact_paths=[],
                config_path=str(config_path),
                database_manifest_path=str(manifest_path),
                database_root=str(root / "db"),
                environment={},
            )

            record = provenance["databases"]["optional"]["files"][0]
            self.assertFalse(record["exists"])
            self.assertIsNone(record["sha256"])


if __name__ == "__main__":
    unittest.main()
