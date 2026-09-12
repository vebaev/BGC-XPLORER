rule parse_mibig_dereplication:
    input:
        antismash="results/{sample}/summary/antismash.bgc.tsv",
        consensus="results/{sample}/summary/consensus_bgcs.tsv",
        index_html="results/{sample}/antismash/index.html",
        regions_js="results/{sample}/antismash/regions.js"
    output:
        "results/{sample}/summary/mibig_dereplication.tsv"
    conda:
        "../envs/report.yaml"
    script:
        "../scripts/parse_mibig.py"
