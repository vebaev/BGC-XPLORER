import os


def bakta_file(wildcards, ext_key):
    base = os.path.join(config["bakta_dir"], wildcards.sample, wildcards.sample)
    return f"{base}{config['bakta_extensions'][ext_key]}"


def bakta_database_path():
    db_type = os.environ.get("BAKTA_DB_TYPE", "light")
    if db_type not in {"light", "full"}:
        raise ValueError("BAKTA_DB_TYPE must be 'light' or 'full'")
    return os.path.join(config["tools"]["bakta"].get("db_root", "/db/bakta"), "db-light" if db_type == "light" else "db")


rule run_bakta:
    input:
        fasta=lambda wc: os.path.join(config.get("fasta_dir", "data/fasta"), wc.sample + ".fasta")
    output:
        gbff=config["bakta_dir"] + "/{sample}/{sample}" + config["bakta_extensions"]["gbff"],
        fna=config["bakta_dir"] + "/{sample}/{sample}" + config["bakta_extensions"]["fna"],
        faa=config["bakta_dir"] + "/{sample}/{sample}" + config["bakta_extensions"]["faa"],
        gff3=config["bakta_dir"] + "/{sample}/{sample}" + config["bakta_extensions"]["gff3"],
        json=config["bakta_dir"] + "/{sample}/{sample}" + config["bakta_extensions"]["json"],
        tsv=config["bakta_dir"] + "/{sample}/{sample}" + config["bakta_extensions"]["tsv"],
    params:
        outdir=lambda wc: os.path.join(config["bakta_dir"], wc.sample),
        executable=lambda wc: config["tools"]["bakta"]["executable"],
        database=lambda wc: bakta_database_path(),
        extra=lambda wc: config["tools"]["bakta"].get("extra_args", "")
    threads: lambda wc: int(config["tools"]["bakta"].get("threads", 4))
    shell:
        """
        mkdir -p {params.outdir}
        {params.executable} --db {params.database} --output {params.outdir} \
          --prefix {wildcards.sample} --threads {threads} --force {params.extra} {input.fasta}
        """


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
