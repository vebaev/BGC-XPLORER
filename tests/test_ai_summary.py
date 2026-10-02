import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from ai_summary import _vocabulary, add_definitions, check_hypothesis, check_sentence, fact, result  # noqa: E402

FACTS = [
    fact("O1", "overview", "Soil_1_13 is a nonribosomal peptide (NRPS) cluster, found by antiSMASH, DeepBGC and GECCO "
         "(all three detection tools).", values=["antiSMASH", "DeepBGC", "GECCO", 3], names=["NRPS"]),
    fact("G2", "genes", "CJLEIP_01147 (DltA) is named as a core gene by antiSMASH, DeepBGC and GECCO; antiSMASH finds "
         "3 modules in it, and predicts phenylalanine, no substrate and valine for its adenylation domains.",
         genes=["CJLEIP_01147"], names=["phenylalanine", "valine"]),
    fact("G3", "genes", "CJLEIP_01149 (DltA) is named as a core gene by antiSMASH, DeepBGC and GECCO; antiSMASH finds "
         "3 modules in it, and predicts no substrate, threonine and threonine for its adenylation domains.",
         genes=["CJLEIP_01149"], names=["threonine"]),
    fact("S1", "similar", "The closest MIBiG entry is BGC0002358 (cyclofaulknamycin), from Streptomyces albidoflavus.",
         names=["cyclofaulknamycin", "Streptomyces albidoflavus"]),
]
VOCAB = _vocabulary(FACTS, known_compounds={"surugamide", "cyclofaulknamycin"})


def check(text, ids):
    return check_sentence({"text": text, "facts": ids}, FACTS, VOCAB)


class SentenceCheck(unittest.TestCase):

    def test_faithful_sentences_pass(self):
        self.assertEqual(check("Soil_1_13 is an NRPS cluster that all three detection tools found.", ["O1"]), [])
        self.assertEqual(check("CJLEIP_01147 has three modules, predicted to load phenylalanine, an unknown residue "
                               "and valine.", ["G2"]), [])
        self.assertEqual(check("Its closest known relative is the cyclofaulknamycin cluster of Streptomyces "
                               "albidoflavus (BGC0002358).", ["S1"]), [])

    def test_sentence_must_cite_a_fact(self):
        self.assertTrue(check("Soil_1_13 is an NRPS cluster.", []))
        self.assertTrue(check("Soil_1_13 is an NRPS cluster.", ["X9"]))

    def test_added_number_fails(self):
        self.assertTrue(check("CJLEIP_01147 has four modules.", ["G2"]))

    def test_identifier_not_in_cited_facts_fails(self):
        self.assertTrue(check("CJLEIP_01149 has three modules.", ["G2"]))
        self.assertTrue(check("It resembles BGC0000001.", ["S1"]))

    def test_substrate_from_another_gene_fails(self):
        self.assertTrue(check("CJLEIP_01147 is predicted to load threonine.", ["G2"]))
        self.assertTrue(check("CJLEIP_01147 is predicted to load threonine.", ["G2", "G3"]))
        self.assertTrue(check("CJLEIP_01147 is predicted to load Thr.", ["G2"]))

    def test_compound_not_in_cited_facts_fails(self):
        self.assertTrue(check("The cluster resembles the surugamide cluster.", ["S1"]))

    def test_caller_not_in_cited_facts_fails(self):
        self.assertTrue(check("Its closest MIBiG entry, BGC0002358, was found by GECCO.", ["S1"]))

    def test_interpretation_fails(self):
        self.assertTrue(check("CJLEIP_01147 likely makes a peptide with phenylalanine.", ["G2"]))
        self.assertTrue(check("Soil_1_13 probably produces an antibiotic.", ["O1"]))


class AuditCases(unittest.TestCase):
    """Sentences that passed an earlier version of the check although they were wrong."""

    COUNTS = [fact("G1", "genes", "By the annotation categories of BGC-XPLORER, it holds 2 regulatory, 2 transport, "
                   "2 tailoring-enzyme and 3 resistance genes.")]
    CORE = [fact("G3", "genes", "CJLEIP_06060, annotated by Bakta as Dimethylallyltransferase, is named as a core gene by "
                 "antiSMASH.", genes=["CJLEIP_06060"], callers=["antismash"]),
            fact("G4", "genes", "CJLEIP_06067, annotated by Bakta as Polyketide synthase, is named as a core gene by "
                 "antiSMASH, DeepBGC and GECCO.", genes=["CJLEIP_06067"], callers=["antismash", "deepbgc", "gecco"])]
    MATCH = [fact("S2", "similar", "antiSMASH's ClusterCompare gives this entry a score of 0.47 on a scale of 0 to 1. "
                  "Overall this is a distant match.")]

    def run_check(self, text, facts, section="genes"):
        return check_sentence({"text": text, "facts": [f["id"] for f in facts]}, facts, _vocabulary(facts, set()), section)

    def test_number_counting_another_category_fails(self):
        self.assertTrue(self.run_check("The cluster contains 3 biosynthetic genes, all tailoring enzymes.", self.COUNTS))
        self.assertTrue(self.run_check("It also holds two regulatory, two transport and two resistance genes.", self.COUNTS))
        self.assertEqual(self.run_check("It holds three resistance genes and two regulatory genes.", self.COUNTS), [])

    def test_callers_moved_to_another_gene_fail(self):
        self.assertTrue(self.run_check("The core genes CJLEIP_06060 and CJLEIP_06067 are named by antiSMASH, DeepBGC "
                                       "and GECCO.", self.CORE))
        self.assertEqual(self.run_check("CJLEIP_06067 is named by antiSMASH, DeepBGC and GECCO, and CJLEIP_06060 by "
                                        "antiSMASH.", self.CORE), [])

    def test_total_given_to_one_gene_fails(self):
        facts = [fact("G3", "genes", "CJLEIP_01147 is named as a core gene by antiSMASH; antiSMASH finds 3 modules in it; "
                      "1 of its modules carries an epimerization domain.", genes=["CJLEIP_01147"], callers=["antismash"]),
                 fact("G6", "genes", "Together CJLEIP_01147 and CJLEIP_01149 hold 6 modules; 3 modules carry an "
                      "epimerization domain.", genes=["CJLEIP_01147", "CJLEIP_01149"])]
        self.assertTrue(self.run_check("CJLEIP_01147 has three modules, and three modules carry an epimerization domain.", facts))
        self.assertEqual(self.run_check("CJLEIP_01147 has three modules, one of which carries an epimerization domain.", facts), [])
        self.assertEqual(self.run_check("CJLEIP_01147 and CJLEIP_01149 hold six modules, three with an epimerization "
                                        "domain.", facts), [])

    def test_closeness_judged_by_the_model_fails(self):
        self.assertTrue(self.run_check("The match is moderate.", self.MATCH, "similar"))
        self.assertEqual(self.run_check("ClusterCompare scores it 0.47 out of 1, a distant match.", self.MATCH, "similar"), [])


class Hypothesis(unittest.TestCase):

    def test_hedged_hypothesis_passes(self):
        self.assertEqual(check_hypothesis({"text": "The cluster may make a peptide related to cyclofaulknamycin.",
                                           "facts": ["S1"]}, FACTS), [])

    def test_unhedged_or_invented_hypothesis_fails(self):
        self.assertTrue(check_hypothesis({"text": "The cluster makes cyclofaulknamycin.", "facts": []}, FACTS))
        self.assertTrue(check_hypothesis({"text": "It may resemble BGC0000001 with 12 genes.", "facts": []}, FACTS))


class Result(unittest.TestCase):

    def test_failed_sentences_are_removed_and_empty_parts_filled_from_facts(self):
        analysis = {"overview": [{"text": "Soil_1_13 is an NRPS cluster that all three detection tools found.",
                                  "facts": ["O1"]}],
                    "genes": [{"text": "CJLEIP_01147 has four modules.", "facts": ["G2"]}],
                    "similar": [], "hypothesis": {"text": "", "facts": []}}
        shown = result(analysis, FACTS)
        sections = {s["key"]: s["sentences"] for s in shown["summary_sections"]}
        self.assertEqual(sections["overview"][0]["source"], "model")
        self.assertEqual({s["source"] for s in sections["genes"]}, {"template"})
        self.assertEqual(shown["sentences_removed"], 1)
        self.assertFalse(shown["hypothesis"]["shown"])

    def test_definitions_are_added_once(self):
        items = add_definitions([{"text": "Each gene has three modules."}, {"text": "The modules differ."}])
        self.assertIn("(the unit of", items[0]["shown_text"])
        self.assertNotIn("(the unit of", items[1]["shown_text"])


if __name__ == "__main__":
    unittest.main()
