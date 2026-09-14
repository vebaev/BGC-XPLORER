import unittest
from pathlib import Path

from scripts.ai_cluster_server import DEFAULT_MODEL, api_key_state


EXPECTED_MODEL = "nvidia/nemotron-3.5-lightning-30b-a3b"


class AIRuntimeConfigTests(unittest.TestCase):
    def test_default_model_is_consistent_across_runtime_entrypoints(self):
        self.assertEqual(DEFAULT_MODEL, EXPECTED_MODEL)
        root = Path(__file__).resolve().parents[1]
        for relative in (
            "docker-compose.yml",
            "scripts/start_ai_cluster_server.sh",
            "scripts/save_ai_key.sh",
            ".env.example",
        ):
            self.assertIn(EXPECTED_MODEL, (root / relative).read_text(encoding="utf-8"), relative)
            self.assertNotIn("deepseek-ai/deepseek-v4-pro", (root / relative).read_text(encoding="utf-8"), relative)

    def test_api_key_state_rejects_empty_and_known_placeholder_values(self):
        self.assertEqual(api_key_state(""), "missing")
        self.assertEqual(api_key_state("validation-only"), "placeholder")
        self.assertEqual(api_key_state("replace-with-your-nvidia-api-key"), "placeholder")
        self.assertEqual(api_key_state("configured-test-key"), "configured")


if __name__ == "__main__":
    unittest.main()
