# AI summaries for three showcase loci (Supplementary S4)

Summaries returned by the optional AI layer for three loci of the showcase
genome (STR-S001, GenBank CP158056.1, sample `Soil_1`), with the evidence sent
to the model and a statement-by-statement check against that evidence.

| Locus | Kind | Statements |
|---|---|---|
| Soil_1_13 | cyclofaulknamycin NRPS locus, all three callers, MIBiG match | 28 |
| Soil_1_15 | geosmin terpene locus, all three callers, MIBiG match | 30 |
| Soil_1_10 | DeepBGC-only locus, no MIBiG comparison | 34 |

Model `nvidia/nemotron-3-ultra-550b-a55b`, prompt version 2.0, temperature 0.1,
top-p 0.95, at most 1,800 tokens, seed 42, guided JSON; generated on
1 October 2026. The SHA-256 request fingerprint stored in each summary matches
the payload file next to it.

| File | Content |
|---|---|
| `Soil_1_consensus_<n>.json` | the stored summary: model, date, prompt version, fingerprint, analysis, traceability result |
| `Soil_1_consensus_<n>.payload.json` | the evidence sent to the model |
| `prompt_v2.0.md` | the system prompt |
| `S4_sentence_check.tsv` | every statement with its verdict and a note |

Verdicts: supported by the evidence; correct background knowledge, not in the
evidence; hypothesis, stated as such; partly incorrect or overstated;
unsupported or incorrect; recommendation. Over the 92 statements: 50 supported,
2 background, 17 hypotheses, 11 partly incorrect or overstated, 3 incorrect,
9 recommendations. No locus tag or MIBiG accession was invented.
