rule run_antismash:
    input:
        gbff=lambda wc: bakta_file(wc, "gbff"),
        qc="results/{sample}/qc/bakta_input_check.json"
    output:
        done=touch("results/{sample}/antismash/.done")
    params:
        outdir="results/{sample}/antismash",
        extra=lambda wc: config["tools"]["antismash"]["extra_args"],
        mode=lambda wc: config["execution"]["mode"],
        reuse_existing=lambda wc: str(config["execution"].get("reuse_existing_outputs", True)).lower(),
        image=lambda wc: tool_container("antismash"),
        work_root=lambda wc: config["execution"]["work_root"],
        db_root=lambda wc: config["execution"]["container_db_root"],
        db_dir=lambda wc: db_mount_dir("antismash"),
        executable=lambda wc: config["tools"]["antismash"].get("executable", "")
    shell:
        """
        mkdir -p {params.outdir}
        if [ "{params.mode}" = "mock" ]; then
          printf 'contig\tstart\tend\tproduct\tregion_number\ncontig_1\t1000\t15000\tNRPS\t1\ncontig_1\t22000\t34000\tRiPP\t2\n' > {params.outdir}/{wildcards.sample}.regions.tsv
        elif [ "{params.mode}" = "docker" ]; then
          docker run --rm \
            --user "$(id -u):$(id -g)" \
            --entrypoint antismash \
            -v "$(pwd):{params.work_root}" \
            -v "{params.db_dir}:{params.db_root}/antismash" \
            {params.image} \
            {params.work_root}/{input.gbff} \
            --output-dir {params.work_root}/{params.outdir} \
            --databases {params.db_root}/antismash \
            {params.extra}
        else
          TOOL_BIN="{params.executable}"
          if [ -z "$TOOL_BIN" ]; then
            TOOL_BIN="$(command -v antismash || true)"
          fi
          if [ "{params.reuse_existing}" = "true" ] && [ -f "{params.outdir}/index.html" ] && [ -f "{params.outdir}/regions.js" ]; then
            printf 'Reusing existing antiSMASH output in %s\n' "{params.outdir}"
          elif [ -n "$TOOL_BIN" ] && [ -x "$TOOL_BIN" ]; then
            "$TOOL_BIN" {input.gbff} --output-dir {params.outdir} {params.extra}
          else
            printf 'antiSMASH executable was not found and no reusable output exists in %s\n' "{params.outdir}" >&2
            exit 1
          fi
        fi
        touch {output.done}
        """
