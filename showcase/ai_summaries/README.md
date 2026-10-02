# AI summaries for ten showcase loci (Supplementary S4)

Summaries returned by the optional AI layer (`AI summary`, mode `summary`) for
ten loci of the showcase genome (STR-S001, GenBank CP158056.1, sample
`Soil_1`), with the facts the model was given, the software's check of every
sentence, the separate AI hypothesis, and a manual reading of every sentence.

Model `nvidia/nemotron-3-super-120b-a12b`, prompt version 5.6-summary,
temperature 0.2, top-p 0.95, at most 2,500 tokens, seed 42, reasoning off,
structured JSON output; generated on 2 October 2026, one request per locus,
10 s apart. The SHA-256 request fingerprint stored in each file is recomputed
from the facts and the prompt in `scripts/ai_summary.py` and matches.

## How a summary is made

1. Code computes numbered facts for the locus from the run's own outputs:
   location and callers; gene roles, NRPS/PKS modules, substrate predictions
   and KnownClusterBlast gene pairs from antiSMASH's JSON; the closest
   cluster's MIBiG 4.0 entry (compound, organism, class, recorded activity,
   evidence, reference). Words that judge a MIBiG match ("close",
   "moderate", "distant") are set by code from the share of genes matched
   and their identity.
2. The model writes three short parts (what the cluster is, what it
   contains, what it resembles) and cites, for every sentence, the facts it
   uses.
3. Each sentence is checked on its own: every number, locus tag, accession,
   caller, substrate and compound name (including every compound name in
   MIBiG) must come from the facts it cites; a number that counts a category
   must count the same category there; a gene's attributes and the callers
   that name it must come from a fact about that gene; no interpretive
   wording. A sentence that fails is removed; a part left empty is filled
   with the facts as written by code. Definitions of technical terms come
   from a fixed glossary, added by code.
4. A hypothesis about what the cluster may make or do is kept apart, shown
   in its own box labelled "Hypothesis (AI-generated, not verified)", and
   checked only for invented identifiers and numbers and for being worded as
   a hypothesis.

## Loci

| Locus | Kind | Sentences shown / removed |
|---|---|---|
| Soil_1_13 | NRPS, all three callers, close MIBiG match (cyclofaulknamycin) | 9 / 0 |
| Soil_1_32 | NRPS, all three callers, close MIBiG match (surugamide) | 9 / 1 |
| Soil_1_48 | 272 kb multi-class locus, all three callers, close MIBiG match (candicidin) | 18 / 0 |
| Soil_1_15 | terpene, all three callers, moderate MIBiG match (geosmin) | 6 / 0 |
| Soil_1_6 | RiPP-like, antiSMASH only, distant MIBiG match (colicin V) | 6 / 0 |
| Soil_1_33 | NRPS, DeepBGC and GECCO, no MIBiG comparison | 5 / 1 |
| Soil_1_22 | no class, DeepBGC and GECCO, no core gene, no MIBiG comparison | 4 / 0 |
| Soil_1_10 | DeepBGC only, no MIBiG comparison | 6 / 0 |
| Soil_1_12 | GECCO only, no core gene, no MIBiG comparison | 6 / 0 |
| Soil_1_38 | DeepBGC and GECCO, core genes only outside the agreed span | 6 / 0 |

Manual reading of all 77 model sentences: 75 shown, of which 74 correct and
1 imprecisely worded (Soil_1_15, "identified as such"); 2 removed by the
check, both of them correct but citing a number not in their cited facts.
No shown sentence was wrong. All ten hypotheses were shown; they are not
judged as facts. The manual reading is the authors' and is recorded per
sentence in `S4_sentence_check.tsv`.

| File | Content |
|---|---|
| `Soil_1_consensus_<n>.json` | the stored result: model, date, prompt version, parameters, fingerprint, facts, the model's answer, every sentence check, the shown text, the hypothesis |
| `prompt_5.6-summary.md` | system prompt, an example user prompt, output schema, glossary |
| `S4_sentence_check.tsv` | every sentence and hypothesis with the check's verdict and the manual reading |
| `v2_interpretation/` | the earlier free-interpretation design (prompt 2.0), kept for the record: of its 92 statements 11 were partly incorrect or overstated and 3 incorrect, which is why it was replaced |
