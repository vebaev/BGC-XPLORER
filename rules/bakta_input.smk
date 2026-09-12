import os


def bakta_file(wildcards, ext_key):
    base = os.path.join(config["bakta_dir"], wildcards.sample, wildcards.sample)
    return f"{base}{config['bakta_extensions'][ext_key]}"


rule validate_bakta_input:
    input:
        gbff=lambda wc: bakta_file(wc, "gbff"),
        fna=lambda wc: bakta_file(wc, "fna"),
        faa=lambda wc: bakta_file(wc, "faa"),
        gff3=lambda wc: bakta_file(wc, "gff3"),
        json=lambda wc: bakta_file(wc, "json"),
        tsv=lambda wc: bakta_file(wc, "tsv"),
    output:
        "results/{sample}/qc/bakta_input_check.json"
    params:
        sample="{sample}"
    conda:
        "../envs/report.yaml"
    script:
        "../scripts/check_bakta_output.py"


rule extract_bakta_metadata:
    input:
        gbff=lambda wc: bakta_file(wc, "gbff"),
        faa=lambda wc: bakta_file(wc, "faa"),
        gff3=lambda wc: bakta_file(wc, "gff3"),
        json=lambda wc: bakta_file(wc, "json"),
        tsv=lambda wc: bakta_file(wc, "tsv"),
        qc="results/{sample}/qc/bakta_input_check.json"
    output:
        features="results/{sample}/bakta/{sample}.features.tsv",
        metadata="results/{sample}/bakta/{sample}.metadata.json"
    params:
        sample="{sample}"
    conda:
        "../envs/report.yaml"
    script:
        "../scripts/parse_bakta.py"
