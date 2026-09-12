rule parse_antismash:
    input:
        done="results/{sample}/antismash/.done"
    output:
        "results/{sample}/summary/antismash.bgc.tsv"
    conda:
        "../envs/report.yaml"
    script:
        "../scripts/parse_antismash.py"


rule parse_gecco:
    input:
        done="results/{sample}/gecco/.done"
    output:
        "results/{sample}/summary/gecco.bgc.tsv"
    conda:
        "../envs/report.yaml"
    script:
        "../scripts/parse_gecco.py"


rule parse_deepbgc:
    input:
        done="results/{sample}/deepbgc/.done"
    output:
        "results/{sample}/summary/deepbgc.bgc.tsv"
    conda:
        "../envs/report.yaml"
    script:
        "../scripts/parse_deepbgc.py"


rule parse_arts:
    input:
        done="results/{sample}/arts/.done"
    output:
        "results/{sample}/summary/arts.hits.tsv"
    conda:
        "../envs/report.yaml"
    script:
        "../scripts/parse_arts.py"


rule parse_eggnog:
    input:
        done="results/{sample}/eggnog/.done",
        bakta="results/{sample}/bakta/{sample}.features.tsv"
    output:
        "results/{sample}/summary/eggnog.annotations.tsv"
    conda:
        "../envs/report.yaml"
    script:
        "../scripts/parse_eggnog.py"


rule parse_dbcan:
    input:
        done="results/{sample}/dbcan/.done"
    output:
        "results/{sample}/summary/dbcan.cgc.tsv"
    conda:
        "../envs/report.yaml"
    script:
        "../scripts/parse_dbcan.py"
