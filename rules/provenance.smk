rule build_provenance:
    input:
        raw=lambda wc: [os.path.join(config.get("fasta_dir", "data/fasta"), wc.sample + ".fasta")] + [bakta_file(wc, key) for key in ("gbff", "fna", "faa", "gff3", "json", "tsv")],
        artifacts=[
            "results/{sample}/qc/bakta_input_check.json",
            "results/{sample}/bakta/{sample}.metadata.json",
            "results/{sample}/summary/antismash.bgc.tsv",
            "results/{sample}/summary/gecco.bgc.tsv",
            "results/{sample}/summary/deepbgc.bgc.tsv",
            "results/{sample}/summary/arts.hits.tsv",
            "results/{sample}/summary/dbcan.cgc.tsv",
            "results/{sample}/summary/consensus_bgcs.tsv",
            "results/{sample}/summary/region_evidence.tsv",
            "results/{sample}/summary/tool_overlap.tsv",
        ],
        config="config/config.yaml",
        database_manifest=database_manifest_path
    output:
        "results/{sample}/provenance.json"
    params:
        database_root=lambda wc: config["execution"]["host_db_root"]
    conda:
        "../envs/report.yaml"
    script:
        "../scripts/build_provenance.py"
