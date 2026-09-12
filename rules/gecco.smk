rule run_gecco:
    input:
        fna=lambda wc: bakta_file(wc, "fna"),
        qc="results/{sample}/qc/bakta_input_check.json"
    output:
        done=touch("results/{sample}/gecco/.done")
    params:
        outdir="results/{sample}/gecco",
        extra=lambda wc: config["tools"]["gecco"]["extra_args"],
        mode=lambda wc: config["execution"]["mode"],
        reuse_existing=lambda wc: str(config["execution"].get("reuse_existing_outputs", True)).lower(),
        image=lambda wc: tool_container("gecco"),
        work_root=lambda wc: config["execution"]["work_root"],
        executable=lambda wc: config["tools"]["gecco"].get("executable", "")
    shell:
        """
        mkdir -p {params.outdir}
        if [ "{params.mode}" = "mock" ]; then
          printf 'sequence_id\tstart\tend\tcluster_id\ttype\taverage_p\tmax_p\ncontig_1\t1200\t14800\tgecco_1\tNRPS\t0.81\t0.92\ncontig_1\t50000\t61000\tgecco_2\tTerpene\t0.73\t0.78\n' > {params.outdir}/{wildcards.sample}.clusters.tsv
        elif [ "{params.mode}" = "docker" ]; then
          docker run --rm \
            --user "$(id -u):$(id -g)" \
            -v "$(pwd):{params.work_root}" \
            {params.image} \
            gecco -v run \
            --genome {params.work_root}/{input.fna} \
            --output-dir {params.work_root}/{params.outdir} \
            {params.extra}
        else
          TOOL_BIN="{params.executable}"
          if [ -z "$TOOL_BIN" ]; then
            TOOL_BIN="$(command -v gecco || true)"
          fi
          if [ "{params.reuse_existing}" = "true" ] && ls "{params.outdir}"/*.clusters.tsv >/dev/null 2>&1; then
            printf 'Reusing existing GECCO output in %s\n' "{params.outdir}"
          elif [ -n "$TOOL_BIN" ] && [ -x "$TOOL_BIN" ]; then
            "$TOOL_BIN" -v run --genome {input.fna} --output-dir {params.outdir} {params.extra}
          else
            printf 'GECCO executable was not found and no reusable output exists in %s\n' "{params.outdir}" >&2
            exit 1
          fi
        fi
        touch {output.done}
        """
