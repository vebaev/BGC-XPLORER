rule aggregate_region_evidence:
    input:
        consensus="results/{sample}/summary/consensus_bgcs.tsv",
        overlap="results/{sample}/summary/tool_overlap.tsv",
        arts="results/{sample}/summary/arts.hits.tsv",
        mibig="results/{sample}/summary/mibig_dereplication.tsv"
    output:
        evidence="results/{sample}/summary/region_evidence.tsv"
    conda:
        "../envs/report.yaml"
    script:
        "../scripts/aggregate_region_evidence.py"
