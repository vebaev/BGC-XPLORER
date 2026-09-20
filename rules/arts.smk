def arts_reference_set(wildcards):
    row = SAMPLE_BY_NAME[wildcards.sample]
    return row["taxon"] or config["tools"]["arts"]["reference_set"]


rule run_arts:
    input:
        gbff=lambda wc: bakta_file(wc, "gbff"),
        antismash_done="results/{sample}/antismash/.done",
        qc="results/{sample}/qc/bakta_input_check.json"
    output:
        done=touch("results/{sample}/arts/.done")
    params:
        outdir="results/{sample}/arts",
        antismash_dir="results/{sample}/antismash",
        reference=arts_reference_set,
        extra=lambda wc: config["tools"]["arts"]["extra_args"],
        mode=lambda wc: config["execution"]["mode"],
        reuse_existing=lambda wc: str(config["execution"].get("reuse_existing_outputs", True)).lower(),
        image=lambda wc: tool_container("arts"),
        work_root=lambda wc: config["execution"]["work_root"],
        db_root=lambda wc: config["execution"]["container_db_root"],
        ref_dir=arts_ref_dir,
        arts_root_dir=arts_root_dir,
        python_bin=lambda wc: config["tools"]["arts"].get("python", ""),
        script_path=lambda wc: config["tools"]["arts"].get("script", "/opt/arts/artspipeline1.py")
    shell:
        """
        mkdir -p {params.outdir}
        if [ "{params.mode}" = "mock" ]; then
          printf 'contig\tstart\tend\tgene\tcategory\tscore\ncontig_1\t1100\t3000\tresistance_gene_A\tduplication\t85\ncontig_2\t8000\t12000\ttransporter_B\tcore_gene\t71\n' > {params.outdir}/{wildcards.sample}.hits.tsv
        elif [ "{params.mode}" = "docker" ]; then
          docker run --rm \
            -v "$(pwd):{params.work_root}" \
            -v "{params.ref_dir}:{params.db_root}/arts_ref" \
            -v "{params.arts_root_dir}:{params.db_root}/arts" \
            {params.image} \
            {params.work_root}/{input.gbff} \
            {params.db_root}/arts_ref \
            -rd {params.work_root}/{params.outdir} \
            -asp {params.work_root}/{params.antismash_dir} \
            -khmms {params.db_root}/arts/knownresistance.hmm \
            -duf {params.db_root}/arts/dufmodels.hmm \
            {params.extra}
        else
          PYTHON_BIN="{params.python_bin}"
          if [ -z "$PYTHON_BIN" ]; then
            PYTHON_BIN="$(command -v python || command -v python3 || true)"
          fi
          if [ "{params.reuse_existing}" = "true" ] && [ -f "{params.outdir}/tables/knownhits.tsv" ] && [ -f "{params.outdir}/tables/coretable.tsv" ]; then
            printf 'Reusing existing ARTS output in %s\n' "{params.outdir}"
          elif [ -n "$PYTHON_BIN" ] && [ -x "$PYTHON_BIN" ] && [ -f "{params.script_path}" ]; then
            PATH=/usr/local/bin:$PATH "$PYTHON_BIN" "{params.script_path}" {input.gbff} {params.ref_dir} -rd {params.outdir} -asp {params.antismash_dir} -khmms {params.arts_root_dir}/knownresistance.hmm -duf {params.arts_root_dir}/dufmodels.hmm {params.extra}
          else
            printf 'ARTS runtime was not found and no reusable output exists in %s\n' "{params.outdir}" >&2
            exit 1
          fi
        fi
        touch {output.done}
        """
