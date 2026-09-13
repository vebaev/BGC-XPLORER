/bin/bash: warning: setlocale: LC_ALL: cannot change locale (C.UTF-8)
rule run_dbcan:
    input:
        faa=lambda wc: bakta_file(wc, "faa"),
        gff3=lambda wc: bakta_file(wc, "gff3"),
        qc="results/{sample}/qc/bakta_input_check.json"
    output:
        done=touch("results/{sample}/dbcan/.done")
    threads: BGC_THREADS
    params:
        outdir="results/{sample}/dbcan",
        extra=lambda wc: config["tools"]["dbcan"]["extra_args"],
        mode=lambda wc: config["execution"]["mode"],
        reuse_existing=lambda wc: str(config["execution"].get("reuse_existing_outputs", True)).lower(),
        image=lambda wc: tool_container("dbcan"),
        work_root=lambda wc: config["execution"]["work_root"],
        db_root=lambda wc: config["execution"]["container_db_root"],
        db_dir=lambda wc: db_mount_dir("dbcan"),
        binary=lambda wc: config["tools"]["dbcan"].get("binary", "run_dbcan"),
        hmmer_bin_dir=lambda wc: config["tools"]["dbcan"].get("hmmer_bin_dir", "")
    shell:
        """
        mkdir -p {params.outdir}
        if [ "{params.mode}" = "mock" ]; then
          printf '#cgcid\tPULID\tdbCAN-PUL substrate\tbitscore\tsignature pairs\tdbCAN-sub substrate\tdbCAN-sub substrate score\ncontig_1|CGC1\tPUL0001\tchitin\t185.2\tCAZyme-CAZyme;CAZyme-TC\tchitin\t3.6\n' > {params.outdir}/substrate_prediction.tsv
          printf 'CGC#\tContig ID\tCluster Start\tCluster End\tGenes\tCAZymes\tTC\tTF\tSTP\tSulfatase\tPeptidase\tSignatures\tLength (bp)\nCGC1\tcontig_1\t90000\t118000\t14\t5\t2\t1\t0\t0\t0\t8\t28001\n' > {params.outdir}/cgc_standard_out_summary.tsv
        elif [ "{params.mode}" = "docker" ]; then
          docker run --rm \
            --user "$(id -u):$(id -g)" \
            -v "$(pwd):{params.work_root}" \
            -v "{params.db_dir}:{params.db_root}/dbcan" \
            {params.image} \
            run_dbcan easy_substrate \
            --mode protein \
            --input_raw_data {params.work_root}/{input.faa} \
            --input_gff {params.work_root}/{input.gff3} \
            --output_dir {params.work_root}/{params.outdir} \
            --db_dir {params.db_root}/dbcan \
            --threads {threads} {params.extra}
        else
          if [ -n "{params.hmmer_bin_dir}" ]; then
            export PATH="{params.hmmer_bin_dir}:$PATH"
          fi
          TOOL_BIN="{params.binary}"
          if [ -z "$TOOL_BIN" ]; then
            TOOL_BIN="$(command -v run_dbcan || true)"
          fi
          if [ "{params.reuse_existing}" = "true" ] && ( [ -f "{params.outdir}/substrate_prediction.tsv" ] || [ -f "{params.outdir}/cgc_standard_out_summary.tsv" ] || [ -f "{params.outdir}/cgc_standard_out.tsv" ] ); then
            printf 'Reusing existing dbCAN output in %s\n' "{params.outdir}"
          elif [ -n "$TOOL_BIN" ] && [ -x "$TOOL_BIN" ]; then
            "$TOOL_BIN" easy_substrate \
              --mode protein \
              --input_raw_data {input.faa} \
              --input_gff {input.gff3} \
              --output_dir {params.outdir} \
              --db_dir {params.db_dir} \
              --threads {threads} {params.extra}
          else
            printf 'dbCAN executable was not found and no reusable output exists in %s\n' "{params.outdir}" >&2
            exit 1
          fi
        fi
        touch {output.done}
        """
