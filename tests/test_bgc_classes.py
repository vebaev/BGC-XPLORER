import sys
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from bgc_classes import activity_labels, count_per_locus, harmonised_classes, mibig_class


class BgcClassTests(unittest.TestCase):
    def test_caller_names_for_one_class_map_to_the_same_mibig_class(self):
        for label in ("NRPS", "NRP", "NRPS-like", "thioamide-NRP"):
            self.assertEqual(mibig_class(label), "NRP")
        for label in ("T1PKS", "T3PKS", "transAT-PKS", "Polyketide", "hglE-KS"):
            self.assertEqual(mibig_class(label), "Polyketide")
        for label in ("terpene", "Terpene", "terpene-precursor"):
            self.assertEqual(mibig_class(label), "Terpene")
        for label in ("RiPP", "RiPP-like", "lanthipeptide-class-iii", "lassopeptide", "RRE-containing"):
            self.assertEqual(mibig_class(label), "RiPP")
        self.assertEqual(mibig_class("ectoine"), "Other")
        self.assertEqual(mibig_class("Unknown"), "Unknown")
        self.assertIsNone(mibig_class("antibacterial"))

    def test_hybrids_count_for_both_classes_and_each_class_once(self):
        self.assertEqual(harmonised_classes("NRP;Polyketide,NRPS,T1PKS"), ["NRP", "Polyketide"])
        self.assertEqual(harmonised_classes("NRP-Polyketide,terpene,Terpene"), ["NRP", "Polyketide", "Terpene"])
        self.assertEqual(harmonised_classes(float("nan")), [])

    def test_activity_is_kept_apart_from_class(self):
        self.assertEqual(activity_labels("antibacterial,NRPS,T1PKS"), ["antibacterial"])
        self.assertEqual(activity_labels("antibacterial-cytotoxic,cytotoxic"), ["antibacterial", "cytotoxic"])
        self.assertEqual(activity_labels("terpene"), [])

    def test_a_locus_counts_once_per_label(self):
        series = pd.Series(["NRPS,NRP,NRPS-like", "terpene,Terpene", "NRP;Polyketide", None])
        self.assertEqual(count_per_locus(series, harmonised_classes),
                         [("NRP", 2), ("Polyketide", 1), ("Terpene", 1)])


if __name__ == "__main__":
    unittest.main()
