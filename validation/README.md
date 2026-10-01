# Validation of the caller grouping against MIBiG

Benchmark of the gene-level grouping of antiSMASH, GECCO and DeepBGC
predictions on four reference genomes from two genera, against MIBiG 4.0 loci.
Everything needed to re-score the results is in this directory; re-running the
workflow itself needs the genomes and the reference databases.

## Design

| | |
|---|---|
| Genomes | *Streptomyces coelicolor* A3(2) (NC_003888.3, plasmids SCP1 NC_003903.1 and SCP2 NC_003904.1), *S. albidoflavus* J1074 (CP004370.1), *S. avermitilis* MA-4680 (NC_003155.5, plasmid SAP1 NC_004719.1), *Salinispora tropica* CNB-440 (CP000667.1) |
| Workflow | commit `81d608f` (voting rule) mounted over image `ghcr.io/vebaev/bgc-xplorer:1.1.1`; tool versions are those of v1.1.1 (Bakta 1.12.0 light, antiSMASH 8.0.4, GECCO 0.10.3, DeepBGC 0.1.31). The grouping code is unchanged in v1.2.0 with default settings: every variant below is reproduced byte for byte. |
| Truth | 42 MIBiG 4.0 loci (18 / 7 / 11 / 6). Every active MIBiG entry was aligned to each genome with megablast and accepted at ≥95 % identity over ≥90 % of its length; coordinates come from the alignment, not from MIBiG. Entries at the same place (reciprocal overlap ≥80 %, or nested with a shared compound) were collapsed. `build_truth.py` |
| Prediction sets | each caller alone; **side-by-side** = the three callers' rows concatenated without grouping, as nf-core/funcscan's comBGC reports them; the earlier **containment** rule; **voting** with 1, 2 (default) or 3 callers per gene, agreed span and union span. Variants are produced by the unmodified `scripts/build_consensus.py` (`run_consensus_variant.py`). |
| Metrics | fixed before scoring, in `validate_consensus.py`: BGC-QUAST recovery categories (full ≥95 %, partial 10–95 %, missed <10 %); delimited = Jaccard ≥0.5; bp Jaccard of the best-matching prediction; signed boundary extension; gene-level precision, recall and F1 inside the recovered locus; split, fused, fragments. Locus-level precision is not reported, because MIBiG does not annotate every cluster in these genomes. |
| Statistics | paired by truth locus, pooled over genomes: Wilcoxon signed-rank test and genome-stratified bootstrap 95 % CI (10,000 resamples) of the mean difference |
| Checks | unit tests of the metric definitions on synthetic loci (`test_validate_consensus.py`); BGC-QUAST 1.1.0 (compare-to-reference mode, MIBiG loci as reference) assigns the same status as `validate_consensus.py` to all 294 locus × set pairs; metrics identical in two independently built environments |

Added after the first scoring, and labelled as such in the tables: the
≥2-caller tier as a separate set, and the "any row" score for side-by-side
reporting (mean Jaccard over the rows overlapping the locus, i.e. what a user
gets without knowing which row is right).

## Results

Per genome: [`tables/table_per_genome.md`](tables/table_per_genome.md); per
locus: [`tables/supplementary_per_locus.xlsx`](tables/supplementary_per_locus.xlsx)
and [`figures/per_locus.png`](figures/per_locus.png)
([`figures/per_locus_supp.png`](figures/per_locus_supp.png) without the
side-by-side column).

Pooled over the four genomes (42 MIBiG loci):

| Prediction set | Loci reported | Genome flagged | Detected | Full / partial / missed | Delimited | Median Jaccard | Median gene F1 | Median gene precision |
|---|---|---|---|---|---|---|---|---|
| antiSMASH | 112 | 15.8 % | 40 | 36 / 4 / 2 | 19 | 0.458 | 0.515 | 0.362 |
| GECCO | 147 | 14.0 % | 38 | 34 / 4 / 4 | 23 | 0.582 | 0.667 | 0.614 |
| DeepBGC | 141 | 30.7 % | 32 | 28 / 4 / 10 | 7 | 0.163 | 0.220 | 0.213 |
| Side-by-side, any row | 400 | 37.2 % | 41 | 41 / 0 / 1 | — | 0.445 | — | — |
| Side-by-side, best row (oracle) | 400 | 37.2 % | 41 | 41 / 0 / 1 | 26 | 0.650 | 0.766 | 0.636 |
| Containment (earlier rule) | 247 | 36.9 % | 41 | 40 / 1 / 1 | 12 | 0.269 | 0.360 | 0.222 |
| Voting 2/3, ≥2-caller tier only | 129 | 14.8 % | 38 | 34 / 4 / 4 | 20 | 0.482 | 0.554 | 0.456 |
| **Voting 2/3 (default), agreed span** | **226** | **23.0 %** | **41** | **37 / 4 / 1** | **22** | **0.524** | **0.626** | **0.478** |
| Voting 3/3, agreed span | 264 | 23.7 % | 41 | 34 / 7 / 1 | 27 | 0.660 | 0.782 | 0.667 |
| Voting 2/3, union span (gene-map window) | 226 | 37.2 % | 41 | 41 / 0 / 1 | 8 | 0.225 | 0.289 | 0.170 |

Paired comparisons of voting 2/3 (agreed span) with each other set, Jaccard:

| Other set | Loci better / worse / tied | Mean Δ Jaccard [95 % CI] | Wilcoxon p |
|---|---|---|---|
| Containment (earlier rule) | 32 / 6 / 4 | +0.168 [0.106, 0.234] | 9 × 10⁻⁷ |
| DeepBGC | 34 / 7 / 1 | +0.236 [0.163, 0.312] | 4 × 10⁻⁷ |
| antiSMASH | 22 / 17 / 3 | +0.079 [0.020, 0.145] | 0.12 |
| Side-by-side, any row | 29 / 9 / 4 | +0.054 [0.027, 0.081] | 3 × 10⁻⁴ |
| GECCO | 9 / 26 / 7 | −0.028 [−0.084, 0.036] | 0.02 |
| Side-by-side, best row (oracle) | 4 / 32 / 6 | −0.087 [−0.121, −0.057] | 1 × 10⁻⁷ |
| Voting 3/3 | 8 / 29 / 5 | −0.097 [−0.150, −0.050] | 2 × 10⁻⁴ |

Gene F1 gives the same pattern, and so does the subset of 34 non-nested loci
of at least 5 kb (`metrics/paired_tests.tsv`, `metrics/summary_pooled.tsv`).
The one locus no set detects (ochronotic pigment, *S. avermitilis*, 5 kb) is
missed by all three callers.

## Reproducing

Re-score from the files here (needs the environment in `environment.yaml`):

```bash
cd validation
mamba env create -p .env -f environment.yaml
.env/bin/python -m pytest -q test_validate_consensus.py
.env/bin/python validate_consensus.py --workdir work --variants variants --truth-dir truth \
    --samples S_coelicolor_A3_2,S_albidoflavus_J1074,S_avermitilis_MA4680,S_tropica_CNB440 --out metrics
.env/bin/python make_figures.py --metrics metrics --out figures
.env/bin/python make_tables.py
```

`work/` holds only what scoring reads: the callers' tables and the consensus
table from each workflow run, the Bakta feature tables, the run provenance,
and the contig names and lengths from the Bakta JSON.

From scratch:

1. Download the genomes above and MIBiG 4.0 (`mibig_json_4.0`, `mibig_gbk_4.0`) and build the truth sets with `build_truth.py`.
2. Run the workflow on the genomes with `run_benchmark_container.sh` (set `DB_ROOT` and, if they live elsewhere, `BAKTA_DB_DIR`, `EGGNOG_DB_DIR`, `DBCAN_DB_DIR`; `IMAGE` and `MOUNT_CODE=0` to use a later release).
3. Run `run_validation.sh`: grouping variants, the check that the default variant equals the workflow's own output, scoring and figures.
4. The BGC-QUAST cross-check uses BGC-QUAST 1.1.0 (commit `45f50f2d82cd423d275ef38f9a8c902278b04da9`) in the environment from `environment.bgcquast.yaml`: `bgcquast_crosscheck.py`, with per-locus status captured by `bgcquast_percluster.py`.

## Files

| Path | Content |
|---|---|
| `build_truth.py` | MIBiG truth set per genome |
| `run_benchmark_container.sh` | workflow run on the benchmark genomes |
| `run_consensus_variant.py` | runs `scripts/build_consensus.py` under another grouping rule |
| `run_validation.sh` | variants, scoring and figures |
| `validate_consensus.py`, `test_validate_consensus.py` | metric definitions, scoring and their tests |
| `make_figures.py`, `make_tables.py` | figures and per-genome tables |
| `bgcquast_crosscheck.py`, `bgcquast_percluster.py` | BGC-QUAST cross-check |
| `truth/` | truth sets (and near misses for *S. tropica*) |
| `variants/` | consensus tables for every grouping variant |
| `metrics/` | per-locus, per-genome, pooled and paired results |
| `tables/`, `figures/` | per-genome tables, per-locus workbook, figures |
| `bgcquast/` | BGC-QUAST status per locus and set |
