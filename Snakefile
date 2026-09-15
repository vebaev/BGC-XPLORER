configfile: "config/config.yaml"

import os
import csv

from scripts.thread_config import configured_threads


BGC_THREADS = configured_threads()


with open(config["samples"], "r", newline="") as samples_handle:
    SAMPLES_TABLE = [
        {key: value or "" for key, value in row.items()}
        for row in csv.DictReader(samples_handle, delimiter="\t")
    ]

SAMPLES = [row["sample"] for row in SAMPLES_TABLE]
SAMPLE_BY_NAME = {row["sample"]: row for row in SAMPLES_TABLE}


include: "rules/bakta_input.smk"
include: "rules/runtime.smk"
include: "rules/antismash.smk"
include: "rules/gecco.smk"
include: "rules/deepbgc.smk"
include: "rules/arts.smk"
include: "rules/eggnog.smk"
include: "rules/dbcan.smk"
include: "rules/normalize.smk"
include: "rules/consensus.smk"
include: "rules/dereplication.smk"
include: "rules/evidence.smk"
include: "rules/catalog.smk"
include: "rules/cluster_maps.smk"
include: "rules/provenance.smk"
include: "rules/report.smk"


rule all:
    input:
        expand("results/{sample}/report/{sample}.html", sample=SAMPLES),
        expand("results/{sample}/summary/region_evidence.tsv", sample=SAMPLES),
        expand("results/{sample}/summary/consensus_bgcs.tsv", sample=SAMPLES),
        expand("results/{sample}/summary/mibig_dereplication.tsv", sample=SAMPLES),
        expand("results/{sample}/summary/eggnog.annotations.tsv", sample=SAMPLES),
        expand("results/{sample}/summary/dbcan.cgc.tsv", sample=SAMPLES),
        expand("results/{sample}/summary/cluster_genes.tsv", sample=SAMPLES),
        expand("results/{sample}/summary/cluster_gene_maps.tsv", sample=SAMPLES),
        expand("results/{sample}/summary/cluster_catalog.tsv", sample=SAMPLES),
        expand("results/{sample}/qc/bakta_input_check.json", sample=SAMPLES),
        expand("results/{sample}/provenance.json", sample=SAMPLES)
