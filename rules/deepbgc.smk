rule run_deepbgc:
    input:
        fna=lambda wc: bakta_file(wc, "fna"),
        qc="results/{sample}/qc/bakta_input_check.json"
    output:
        done=touch("results/{sample}/deepbgc/.done")
    params:
        outdir="results/{sample}/deepbgc",
        extra=lambda wc: config["tools"]["deepbgc"]["extra_args"],
        mode=lambda wc: config["execution"]["mode"],
        reuse_existing=lambda wc: str(config["execution"].get("reuse_existing_outputs", True)).lower(),
        image=lambda wc: tool_container("deepbgc"),
        work_root=lambda wc: config["execution"]["work_root"],
        db_root=lambda wc: config["execution"]["container_db_root"],
        db_dir=lambda wc: db_mount_dir("deepbgc"),
        executable=lambda wc: config["tools"]["deepbgc"].get("executable", "")
    shell:
        """
        mkdir -p {params.outdir}
        if [ "{params.mode}" = "mock" ]; then
          printf 'sequence_id\tnucl_start\tnucl_end\tbgc_candidate_id\tproduct_class\tdetector_score\tclassifier_score\ncontig_1\t900\t15200\tdeepbgc_1\tNRPS\t0.88\t0.84\ncontig_2\t7000\t18000\tdeepbgc_2\tPKS\t0.76\t0.71\n' > {params.outdir}/{wildcards.sample}.bgc.tsv
        elif [ "{params.mode}" = "docker" ]; then
          docker run --rm \
            --user "$(id -u):$(id -g)" \
            -e DEEPBGC_DOWNLOADS_DIR={params.db_root}/deepbgc \
            -v "$(pwd):{params.work_root}" \
            -v "{params.db_dir}:{params.db_root}/deepbgc" \
            {params.image} \
            deepbgc pipeline \
            {params.work_root}/{input.fna} \
            --output {params.work_root}/{params.outdir} \
            {params.extra}
        else
          TOOL_BIN="{params.executable}"
          if [ -z "$TOOL_BIN" ]; then
            TOOL_BIN="$(command -v deepbgc || true)"
          fi
          if [ "{params.reuse_existing}" = "true" ] && [ -f "{params.outdir}/deepbgc.bgc.tsv" ]; then
            printf 'Reusing existing DeepBGC output in %s\n' "{params.outdir}"
          elif [ -n "$TOOL_BIN" ] && [ -x "$TOOL_BIN" ]; then
            DEEPBGC_DOWNLOADS_DIR={params.db_dir} "$TOOL_BIN" pipeline {input.fna} --output {params.outdir} {params.extra} || true
          else
            printf 'DeepBGC executable was not found and no reusable output exists in %s\n' "{params.outdir}" >&2
            exit 1
          fi
        fi
        # DeepBGC occasionally crashes during final rename when no BGCs are found.
        # Ensure the expected outputs exist so downstream parsing can proceed.
        touch {params.outdir}/deepbgc.full.gbk {params.outdir}/deepbgc.bgc.gbk
        touch {output.done}
        """
