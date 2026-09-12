rule run_eggnog:
    input:
        faa=lambda wc: bakta_file(wc, "faa"),
        qc="results/{sample}/qc/bakta_input_check.json"
    output:
        done=touch("results/{sample}/eggnog/.done")
    threads: 4
    params:
        outdir="results/{sample}/eggnog",
        extra=lambda wc: config["tools"]["eggnog"]["extra_args"],
        search_mode=lambda wc: config["tools"]["eggnog"].get("search_mode", "diamond"),
        mode=lambda wc: config["execution"]["mode"],
        reuse_existing=lambda wc: str(config["execution"].get("reuse_existing_outputs", True)).lower(),
        image=lambda wc: tool_container("eggnog"),
        work_root=lambda wc: config["execution"]["work_root"],
        db_root=lambda wc: config["execution"]["container_db_root"],
        db_dir=lambda wc: db_mount_dir("eggnog"),
        executable=lambda wc: config["tools"]["eggnog"].get("executable", "")
    shell:
        """
        mkdir -p {params.outdir}
        if [ "{params.mode}" = "mock" ]; then
          printf '## emapper mock output\n#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tPFAMs\nBFJFPI_00001\tseed1\t1e-50\t250\tCOG0001@1|root\t2|Bacteria\tL\tATP-dependent helicase\tHelicase\tGO:0005524\t-\tko:K03628\tPF00271\nBFJFPI_00002\tseed2\t1e-40\t220\tCOG0002@1|root\t2|Bacteria\tK\tTelomere-associated protein\tTap\t-\t-\t-\tPF00000\n' > {params.outdir}/{wildcards.sample}.emapper.annotations
        elif [ "{params.mode}" = "docker" ]; then
          docker run --rm \
            --user "$(id -u):$(id -g)" \
            -v "$(pwd):{params.work_root}" \
            -v "{params.db_dir}:{params.db_root}/eggnog" \
            {params.image} \
            emapper.py \
            -m {params.search_mode} \
            -i {params.work_root}/{input.faa} \
            --itype proteins \
            --data_dir {params.db_root}/eggnog \
            --output_dir {params.work_root}/{params.outdir} \
            -o {wildcards.sample} \
            --cpu {threads} \
            {params.extra}
        else
          TOOL_BIN="{params.executable}"
          if [ -z "$TOOL_BIN" ]; then
            TOOL_BIN="$(command -v emapper.py || true)"
          fi
          if [ "{params.reuse_existing}" = "true" ] && ls "{params.outdir}"/*.emapper.annotations >/dev/null 2>&1; then
            printf 'Reusing existing eggNOG output in %s\n' "{params.outdir}"
          elif [ -n "$TOOL_BIN" ] && [ -x "$TOOL_BIN" ]; then
            "$TOOL_BIN" \
              -m {params.search_mode} \
              -i {input.faa} \
              --itype proteins \
              --data_dir {params.db_dir} \
              --output_dir {params.outdir} \
              -o {wildcards.sample} \
              --cpu {threads} \
              {params.extra}
          else
            printf 'eggNOG-mapper executable was not found and no reusable output exists in %s\n' "{params.outdir}" >&2
            exit 1
          fi
        fi
        touch {output.done}
        """
