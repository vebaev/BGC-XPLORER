import unittest

from scripts.fasta_input import MAX_FASTA_BYTES, fasta_suffix, normalize_sample_name


class FastaInputTests(unittest.TestCase):
    def test_upload_limit_is_30_decimal_megabytes(self):
        self.assertEqual(MAX_FASTA_BYTES, 30_000_000)

    def test_accepts_supported_fasta_names(self):
        for name in ("genome.fa", "genome.fasta", "genome.fna", "GENOME.FASTA"):
            with self.subTest(name=name):
                self.assertTrue(fasta_suffix(name))

    def test_rejects_non_fasta_uploads(self):
        for name in ("genome.gbff", "genome.faa", "genome.fasta.gz", "genome"):
            with self.subTest(name=name):
                self.assertIsNone(fasta_suffix(name))

    def test_normalizes_safe_sample_name(self):
        self.assertEqual(normalize_sample_name(" Soil sample 1 "), "Soil_sample_1")
        self.assertEqual(normalize_sample_name("../../sample"), "sample")
        with self.assertRaises(ValueError):
            normalize_sample_name("---")


if __name__ == "__main__":
    unittest.main()
