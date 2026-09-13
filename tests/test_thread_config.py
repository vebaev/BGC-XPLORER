import os
import unittest
from unittest.mock import patch

from scripts.thread_config import configured_threads


class ThreadConfigTests(unittest.TestCase):
    def test_defaults_to_four_threads(self):
        with patch.dict(os.environ, {}, clear=True):
            self.assertEqual(configured_threads(), 4)

    def test_reads_positive_docker_thread_count(self):
        with patch.dict(os.environ, {"BGC_THREADS": "8"}, clear=True):
            self.assertEqual(configured_threads(), 8)

    def test_rejects_zero_non_integer_and_whitespace(self):
        for value in ("0", "-2", "1.5", "", " 8 "):
            with self.subTest(value=value):
                with patch.dict(os.environ, {"BGC_THREADS": value}, clear=True):
                    with self.assertRaisesRegex(ValueError, "positive integer"):
                        configured_threads()


if __name__ == "__main__":
    unittest.main()
