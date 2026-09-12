rule prioritize_regions:
    input:
        consensus="results/{sample}/summary/consensus_bgcs.tsv",
        overlap="results/{sample}/summary/tool_overlap.tsv",
        arts="results/{sample}/summary/arts.hits.tsv",
        mibig="results/{sample}/summary/mibig_dereplication.tsv"
    output:
        prioritized="results/{sample}/summary/prioritized_bgcs.tsv",
        prioritized_regions="results/{sample}/summary/prioritized_regions.tsv"
    conda:
        "../envs/report.yaml"
    script:
        "../scripts/prioritize_regions.py"
