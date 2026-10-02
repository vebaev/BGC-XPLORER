# AI summary prompt, version 5.6-summary

Model `nvidia/nemotron-3-super-120b-a12b`; temperature 0.2, top-p 0.95, at most 2500 tokens, seed 42, reasoning off (`enable_thinking: false`), structured JSON output (`response_format: json_schema`).

## System prompt

```text
You describe one biosynthetic gene cluster for a biologist, from numbered facts computed by software.
Write three short parts:
- "overview": 1-2 sentences, what the cluster is (class, which tools found it, size).
- "genes": 2-5 sentences, what it contains: how many biosynthetic genes there are and of which kinds, the core genes and
  their make-up (modules, domains, predicted substrates), and the most informative additional genes. Do not list every
  gene; choose. Prefer words a biologist reads easily (e.g. "phenylalanine" over "Phe").
- "similar": 1-3 sentences, the closest known cluster in MIBiG and how close it is, or that there is none.
Rules for these three parts:
- Write for a reader, in natural, varied prose, not a list. Sentences of 15-30 words; one idea per sentence.
- Start "genes" with the overall picture (how many biosynthetic genes, of which kinds), then the core genes (describe at
  most three one by one and use the "Together" fact for the rest), then at most four additional genes, grouped by kind (e.g. "two oxidising enzymes, a cytochrome P450 and an L-asparagine oxygenase").
- Say modules and substrates in words ("three modules, predicted to load phenylalanine, an unidentified residue and
  valine"); give a domain order only if it adds something. Write "an unidentified residue" for "no substrate".
- In "similar", say how close the match is in words first, then the numbers.
- After every sentence give the ids of the facts it uses. Use only values from those facts: keep every number, locus
  tag, accession, organism, compound, substrate and caller name exactly as in the facts (numbers up to twenty may be
  written as words). Do not add numbers, genes, names or comparisons of your own.
- Do not interpret: do not say what the cluster makes, what a gene does in this cluster, or whether it is new.
- Use technical terms such as module, adenylation domain, epimerization domain or thioesterase without defining them;
  the software adds the definitions.
Then a separate "hypothesis": 1-3 sentences on what the cluster might make or do, worded as a hypothesis (may, could,
possibly) and based on the facts and on the gene products listed as context; consider the whole set of genes, not
only the core genes, and say so if the genes point to primary metabolism rather than a specialised metabolite. Cite
the fact ids it rests on. It will be shown labelled as an unverified AI
hypothesis. Do not invent identifiers or numbers in it.
Return ONLY JSON matching the schema.
```

## User prompt

The numbered facts computed by `scripts/ai_summary.py` for one locus, by part, then the locus's gene products as context for the hypothesis only. Example (Soil_1_13):

```text
overview facts:
  [O1] Soil_1_13 is a nonribosomal peptide (NRPS) cluster, found by antiSMASH, DeepBGC and GECCO (all three detection tools).
  [O2] It spans 57 kb on contig_1 (1,332,869–1,390,082) and contains 30 genes.
genes facts:
  [G1] antiSMASH marks 11 genes as biosynthetic (2 core and 9 additional), 4 as regulatory and 3 as transport.
  [G2] 3 core genes are named in the locus, each by at least one of antiSMASH, DeepBGC and GECCO.
  [G3] CJLEIP_01147, annotated by Bakta as D-alanine--poly(phosphoribitol) ligase subunit DltA, is named as a core gene by antiSMASH, DeepBGC and GECCO; antiSMASH finds 3 modules in it, and predicts phenylalanine, no substrate and valine for its adenylation domains; 1 of its modules carries an epimerization domain.
  [G4] CJLEIP_01149, annotated by Bakta as D-alanine--poly(phosphoribitol) ligase subunit DltA, is named as a core gene by antiSMASH, DeepBGC and GECCO; antiSMASH finds 3 modules in it, and predicts no substrate for 1 and threonine for 2 for its adenylation domains; 2 of its modules carry an epimerization domain.
  [G5] CJLEIP_01154, annotated by Bakta as Thioesterase, is named as a core gene by DeepBGC and GECCO.
  [G6] Together CJLEIP_01147 and CJLEIP_01149 hold 6 modules; antiSMASH predicts a substrate for 4 modules (phenylalanine for 1, valine for 1 and threonine for 2) and none for 2 modules; 3 modules carry an epimerization domain.
  [G7] CJLEIP_01139 encodes Glycosyltransferase family 1 protein and is marked by antiSMASH as an additional biosynthetic gene.
  [G8] CJLEIP_01148 encodes MbtH family protein and is marked by antiSMASH as an additional biosynthetic gene.
  [G9] CJLEIP_01152 encodes Cytochrome P450 and is marked by antiSMASH as an additional biosynthetic gene.
  [G10] CJLEIP_01161 encodes 23S rRNA (uracil(1939)-C(5))-methyltransferase RlmD and is marked by antiSMASH as an additional biosynthetic gene.
  [G11] CJLEIP_01136 encodes 2-deoxy-scyllo-inosamine dehydrogenase and is marked by antiSMASH as an additional biosynthetic gene.
similar facts:
  [S1] The closest MIBiG entry is BGC0002358 (cyclofaulknamycin), from Streptomyces albidoflavus.
  [S2] 21 genes of the locus have a counterpart in BGC0002358, at 96–100 % identity; antiSMASH finds counterparts for 100 % of that cluster's genes. Overall this is a close match.
  [S3] CJLEIP_01147 matches XNR_0983 (Peptide synthase) of BGC0002358 at 99 % identity.
  [S4] CJLEIP_01149 matches XNR_0985 (Peptide synthase) of BGC0002358 at 99 % identity.
  [S5] CJLEIP_01154 matches XNR_0990 (putative thioesterase involved in non-ribosomalpeptide biosynthesis) of BGC0002358 at 99 % identity.
  [S6] MIBiG lists the compound class of BGC0002358 as PKS.
  [S7] MIBiG gives knock-out studies as evidence for BGC0002358 (PubMed 34442689).
  [S8] MIBiG rates the annotation quality of BGC0002358 as questionable.
context for the hypothesis only (do not use it in the three parts):
  [H1] Gene products in the locus, in genome order: CJLEIP_01135 Carbohydrate ABC transporter permease; CJLEIP_01136 2-deoxy-scyllo-inosamine dehydrogenase; CJLEIP_01137 DUF2029 domain-containing protein; CJLEIP_01138 Methyltransferase; CJLEIP_01139 Glycosyltransferase family 1 protein; CJLEIP_01140 DNA-binding response regulator; CJLEIP_01141 histidine kinase; CJLEIP_01142 ABC transporter ATP-binding protein; CJLEIP_01143 ABC3 transporter permease C-terminal domain-containing protein; CJLEIP_01144 DNA-binding response regulator; CJLEIP_01145 histidine kinase; CJLEIP_01146 Class A beta-lactamase-related serine hydrolase; CJLEIP_01147 D-alanine--poly(phosphoribitol) ligase subunit DltA; CJLEIP_01148 MbtH family protein; CJLEIP_01149 D-alanine--poly(phosphoribitol) ligase subunit DltA; CJLEIP_01150 L-asparagine oxygenase; CJLEIP_01151 MFS transporter; CJLEIP_01152 Cytochrome P450; CJLEIP_01153 HTH arsR-type domain-containing protein; CJLEIP_01154 Thioesterase; CJLEIP_01155 Erythromycin esterase; CJLEIP_01156 ATP-binding protein; CJLEIP_01157 Cation efflux family protein; CJLEIP_01158 ATP/GTP-binding protein; CJLEIP_01159 IS5 family IS112 transposase; CJLEIP_01160 PE-PGRS family protein; CJLEIP_01161 23S rRNA (uracil(1939)-C(5))-methyltransferase RlmD; CJLEIP_01162 APC family permease; CJLEIP_01163 Trk system potassium uptake protein TrkA; CJLEIP_01164 TrkA family potassium uptake protein.
```

## Output schema

```json
{
  "type": "object",
  "additionalProperties": false,
  "required": [
    "overview",
    "genes",
    "similar",
    "hypothesis"
  ],
  "properties": {
    "overview": {
      "type": "array",
      "minItems": 1,
      "items": {
        "type": "object",
        "additionalProperties": false,
        "required": [
          "text",
          "facts"
        ],
        "properties": {
          "text": {
            "type": "string",
            "minLength": 10
          },
          "facts": {
            "type": "array",
            "minItems": 1,
            "items": {
              "type": "string"
            }
          }
        }
      }
    },
    "genes": {
      "type": "array",
      "minItems": 1,
      "items": {
        "type": "object",
        "additionalProperties": false,
        "required": [
          "text",
          "facts"
        ],
        "properties": {
          "text": {
            "type": "string",
            "minLength": 10
          },
          "facts": {
            "type": "array",
            "minItems": 1,
            "items": {
              "type": "string"
            }
          }
        }
      }
    },
    "similar": {
      "type": "array",
      "minItems": 1,
      "items": {
        "type": "object",
        "additionalProperties": false,
        "required": [
          "text",
          "facts"
        ],
        "properties": {
          "text": {
            "type": "string",
            "minLength": 10
          },
          "facts": {
            "type": "array",
            "minItems": 1,
            "items": {
              "type": "string"
            }
          }
        }
      }
    },
    "hypothesis": {
      "type": "object",
      "additionalProperties": false,
      "required": [
        "text",
        "facts"
      ],
      "properties": {
        "text": {
          "type": "string"
        },
        "facts": {
          "type": "array",
          "items": {
            "type": "string"
          }
        }
      }
    }
  }
}
```

## Glossary added by code

The first use of each term in the shown text gets this definition in brackets; the model is told not to define terms.

| Term | Definition |
|---|---|
| module | the unit of an assembly-line enzyme that adds one building block to the growing chain |
| adenylation domain | selects and activates the amino acid a module adds |
| epimerization domain | converts an amino acid to its D-form |
| condensation domain | joins the amino acid of its module to the growing chain |
| thioesterase | an enzyme type that releases the finished chain from NRPS and PKS assembly lines |
| MbtH-like protein | a small helper protein of NRPS adenylation domains |
| MbtH family protein | a small helper protein of NRPS adenylation domains |
| cytochrome P450 | a family of oxidising enzymes |
| ketosynthase domain | extends a polyketide chain by one unit |
| acyltransferase domain | selects the extender unit a polyketide module adds |
| glycosyltransferase | an enzyme type that attaches sugars |
| methyltransferase | an enzyme type that adds methyl groups |
| terpene synthase | an enzyme type that folds a linear precursor into a terpene skeleton |
