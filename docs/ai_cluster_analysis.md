# On-demand AI cluster analysis

The HTML report can call a local AI analysis service through the `AI analysis`
button next to each cluster. The browser never receives the API key.

## Start the local service

Save the NVIDIA API key once in a local ignored file:

```bash
scripts/save_ai_key.sh
```

The key is written to `config/local_ai.env`, which is ignored by git and should
have permissions `600`.

After the report has been generated, start both the AI service and the local
report web server:

```bash
scripts/serve_report_with_ai.sh
```

Then open:

```text
http://127.0.0.1:8000/SOIL_CONTIGS.html
```

To start only the AI service:

```bash
scripts/start_ai_cluster_server.sh --start
scripts/start_ai_cluster_server.sh --start --model deepseek-ai/deepseek-v4-flash
```

If `config/local_ai.env` does not exist, the script asks for the API key without
echoing it and does not write it to a file. The service is started in the
background, so you can close the terminal afterward.

To stop or restart it later:

```bash
scripts/start_ai_cluster_server.sh --stop
scripts/start_ai_cluster_server.sh --restart
scripts/start_ai_cluster_server.sh --restart --model deepseek-ai/deepseek-v4-flash
```

Alternatively, set the key yourself:

```bash
export NVIDIA_API_KEY="YOUR_KEY_HERE"
scripts/start_ai_cluster_server.sh --start
```

Do not paste the API key into the HTML report or any repository file.

You can override the model per run with `--model ...`, or set `NVIDIA_MODEL`
in the environment / `config/local_ai.env`. The report HTML does not need to be
edited.

The report expects:

```text
http://127.0.0.1:8787/analyze_cluster
```

## Behavior

- Requests are on-demand, one cluster at a time.
- Results are cached under `results/{sample}/ai_cluster_reports/{consensus_id}.json`.
- Cached results are returned without another NVIDIA API call.
- The service throttles calls to at most 40 requests per minute.
- If `NVIDIA_API_KEY` is missing, the report shows a clear local-service error.
- If the model returns plain text or markdown instead of JSON, the service
  converts it to a valid compact summary instead of failing.

## If the report shows `Load failed`

This means the browser could not reach the local service at
`http://127.0.0.1:8787/analyze_cluster`.

First, confirm the service is running:

```bash
scripts/start_ai_cluster_server.sh --start
```

If the service is running but the browser still shows `Load failed`, open the
report through a local HTTP server instead of directly from the filesystem:

```bash
scripts/serve_soil_contigs_report.sh
```

Then open:

```text
http://127.0.0.1:8000/SOIL_CONTIGS.html
```

## Input sent to the model

For the selected `consensus_id`, the service sends structured evidence:

- region coordinates, caller count, class/product signals, ARTS hit counts and representative MIBiG comparison metrics;
- antiSMASH, GECCO and DeepBGC overlapping predictions;
- ARTS overlapping evidence;
- dbCAN CGC overlap when available;
- all genes in the cluster with Bakta, eggNOG, dbCAN and ARTS annotations.

The prompt asks the model to return strict JSON with summary, likely function,
biosynthetic logic, key genes, resistance/transport/regulation and recommended
follow-up. The analysis is an interpretation of supplied evidence, not a calibrated
confidence or novelty estimate.

### Scientific interpretation prompt v2.0

The model returns six sections: summary, likely product/function, biosynthetic
logic, key genes, resistance/transport/regulation, and recommended follow-up.
Uncertainty belongs next to the relevant claim, rather than in separate
confidence, caveats, or novelty sections. Interpretations must reference supplied
locus tags or tool evidence; unsupported compounds, activities and citations are
prohibited. No literature search is performed by this endpoint.

Gene coordinates, strand, EC, KEGG KO, and complete supplied annotations are
retained. Empty evidence tables are not treated as verified negative findings:
tool completion is explicitly unknown in this payload. Workflow interpretations
are labeled as preliminary hypotheses. The cache fingerprint includes the input,
prompt/version, model, endpoint and generation parameters. Responses record the
prompt version and fingerprint. A model response remains a computational
interpretation requiring expert review, not experimental validation.
