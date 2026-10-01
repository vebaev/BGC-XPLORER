# Test dataset

`STR-S001_CP158056.1_1200001-1520000.fasta` is a 320 kb fragment of the
*Streptomyces albidoflavus* STR-S001 chromosome (GenBank
[CP158056.1](https://www.ncbi.nlm.nih.gov/nuccore/CP158056.1), positions
1,200,001–1,520,000). It contains two complete biosynthetic loci, far enough
from the fragment ends that neither is cut, so a full run finishes in minutes
rather than hours and still exercises every caller, the grouping and the
report.

It is meant to check an installation, not as an example of a draft assembly:
BGC-XPLORER is developed and tested for whole bacterial genomes.

## Running it

Upload the FASTA in the web interface, or copy it into `data/fasta/` as
`STR_S001_test.fasta`, add `STR_S001_test	actinobacteria` to
`config/samples.tsv`, and run the workflow.

## Expected result

A complete run gives three loci (`expected_loci.tsv`, written by
`showcase/loci_table.py`; coordinates are on the fragment):

| Locus | Agreed span | Callers | MIBiG representative |
|---|---|---|---|
| STR_S001_test_2 | 132,869–190,082 (57,214 bp) | antiSMASH, GECCO, DeepBGC | BGC0002358.3, cyclofaulknamycin |
| STR_S001_test_3 | 271,431–282,301 (10,871 bp) | antiSMASH, GECCO | BGC0002470.2, synechobactins |
| STR_S001_test_1 | 1–16,869 (16,869 bp) | DeepBGC | none |

The first two are the full-genome loci at the same place. The third, a
DeepBGC call against the start of the fragment, is not found in the
full-genome run: it comes from cutting the chromosome, and it is why
fragmented assemblies are not supported in this version.

In the cyclofaulknamycin locus dbCAN's DIAMOND search alone assigns
glycosyltransferase family GT1 to the two adenylation-domain genes. dbCAN does
not recommend that family, because only one of its three methods finds it, so
the report shows the hit in the gene tooltip as a single-method hit and neither
labels the genes as CAZymes nor flags an annotation conflict.

On a 16-thread run (Intel Xeon Gold 6148) the fragment took 73 minutes, 70 of
them in eggNOG-mapper, whose DIAMOND search scans the whole eggNOG protein
database regardless of input size.
