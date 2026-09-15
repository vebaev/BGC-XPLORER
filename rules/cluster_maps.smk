rule render_cluster_maps:
    input:
        evidence="results/{sample}/summary/region_evidence.tsv",
        dbcan_summary="results/{sample}/summary/dbcan.cgc.tsv",
        bakta="results/{sample}/bakta/{sample}.features.tsv",
        eggnog="results/{sample}/summary/eggnog.annotations.tsv",
        arts="results/{sample}/summary/arts.hits.tsv",
        dbcan_done="results/{sample}/dbcan/.done"
    output:
        cluster_genes="results/{sample}/summary/cluster_genes.tsv",
        cluster_maps="results/{sample}/summary/cluster_gene_maps.tsv",
        map_assets=directory("results/{sample}/report/assets/cluster_maps")
    conda:
        "../envs/report.yaml"
    script:
        "../scripts/render_cluster_maps.py"
