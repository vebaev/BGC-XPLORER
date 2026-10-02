# On-demand AI cluster summary

The HTML report can call a local AI service through the `AI summary`
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

## What the AI summary does (default mode `summary`)

The model is never asked to interpret the raw evidence. For the selected
`consensus_id`, `scripts/ai_summary.py` first computes numbered facts from the
run's own outputs: location and callers; gene roles, NRPS/PKS modules,
substrate predictions and KnownClusterBlast gene pairs from antiSMASH's JSON;
and the closest cluster's MIBiG 4.0 entry (compound, organism, class, recorded
activity, evidence, reference; read from `db/mibig/mibig_json_4.0`, or
`MIBIG_JSON_DIR`). The model writes three short parts (what the cluster is,
what it contains, what it resembles) and cites the facts of every sentence.

Each sentence is then checked on its own: its numbers, locus tags,
accessions, callers, substrates and compound names must come from the facts it
cites, a number that counts a category must count the same one there, a gene's
attributes and callers must come from a fact about that gene, and interpretive
wording is not allowed. Sentences that fail are removed and counted in the
report; a part left empty is filled with the facts as written by code.
Definitions of technical terms come from a fixed glossary added by code.

A separate hypothesis about what the cluster may make or do is shown in its
own box, labelled "Hypothesis (AI-generated, not verified)". It may use all
gene products of the locus as context and is checked only for invented
identifiers and numbers. Every model-written result ends with a notice that AI
models can make mistakes, naming the model.

No literature or web search is performed. Each stored result records the
model, prompt version, parameters, date and a SHA-256 fingerprint of the
request (prompt, facts, model, endpoint and parameters). The prompt, an
example input and a ten-locus evaluation are in `showcase/ai_summaries/`.

Other modes can be selected with `AI_MODE`: `facts` (one paragraph worded from
a fixed fact list), `evidence` (statement-level citations) and
`interpretation` (the earlier free interpretation, prompt 2.0, unverified).
