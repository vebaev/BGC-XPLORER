rule build_consensus:
    input:
        antismash="results/{sample}/summary/antismash.bgc.tsv",
        gecco="results/{sample}/summary/gecco.bgc.tsv",
        deepbgc="results/{sample}/summary/deepbgc.bgc.tsv",
        arts="results/{sample}/summary/arts.hits.tsv",
        bakta="results/{sample}/bakta/{sample}.features.tsv",
        fna=lambda wc: bakta_file(wc, "fna")
    output:
        consensus="results/{sample}/summary/consensus_bgcs.tsv",
        overlap="results/{sample}/summary/tool_overlap.tsv",
        gene_support="results/{sample}/summary/gene_caller_support.tsv"
    conda:
        "../envs/report.yaml"
    script:
        "../scripts/build_consensus.py"
