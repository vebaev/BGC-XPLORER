rule build_cluster_catalog:
    input:
        consensus="results/{sample}/summary/consensus_bgcs.tsv",
        prioritized="results/{sample}/summary/prioritized_bgcs.tsv",
        antismash="results/{sample}/summary/antismash.bgc.tsv",
        gecco="results/{sample}/summary/gecco.bgc.tsv",
        deepbgc="results/{sample}/summary/deepbgc.bgc.tsv"
    output:
        tsv="results/{sample}/summary/cluster_catalog.tsv",
        html="results/{sample}/summary/cluster_catalog.html"
    conda:
        "../envs/report.yaml"
    script:
        "../scripts/build_cluster_catalog.py"
