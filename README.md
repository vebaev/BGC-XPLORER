# BGC Discovery Workflow

Snakemake workflow that starts from precomputed Bakta outputs and runs:

- antiSMASH
- GECCO
- DeepBGC
- ARTS
- eggNOG-mapper

The workflow then normalizes the tool outputs, builds a simple consensus table, prioritizes candidate regions, and renders an HTML summary report.

## Expected input layout

Each sample must exist under `data/bakta/{sample}/` and include:

- `{sample}.gbff`
- `{sample}.fna`
- `{sample}.faa`
- `{sample}.gff3`
- `{sample}.json`
- `{sample}.tsv`

## Docker quick start

Docker Engine with the Compose plugin is the only host dependency. Copy the
example settings, add your NVIDIA API key, and select any model exposed by the
NVIDIA API catalog:

```bash
cp .env.example .env
# Edit .env, then start BGC-XPLORER.
docker compose up -d
```

The default URL is `http://localhost:8778`. Change `BGC_PORT` in `.env` to
publish the app on another host port. For example:

```dotenv
NVIDIA_API_KEY=replace-with-your-nvidia-api-key
NVIDIA_MODEL=nvidia/nemotron-3.5-lightning-30b-a3b
BGC_PORT=9000
AUTO_PREPARE_DATABASES=true
ARTS_REFERENCE=actinobacteria
BGC_IMAGE=vebaev/bgc-xplorer
BGC_VERSION=latest
```

For local development, build the image from the checked-out source:

```bash
docker compose -f docker-compose.yml -f docker-compose.dev.yml up -d --build
```

The API key is passed only to the backend container and `.env` is ignored by
Git. On startup, the container verifies its bundled executables and the
mounted `db/` directory. Existing valid resources are reused. When required
resources are missing and `AUTO_PREPARE_DATABASES=true`, supported database
downloaders run and validation is repeated before the service starts.

ARTS uses a taxon-specific reference set. If the required ARTS files are not
present under `db/arts/actinobacteria`, startup stops with the missing file
list; populate that reference set and run `docker compose up -d` again.

## Snakemake quick start

1. Edit `config/samples.tsv`
2. Edit `config/config.yaml` if you need custom arguments
3. Run:

```bash
snakemake --use-conda --cores 4
```

## Execution modes

The default `config/config.yaml` uses:

```yaml
execution:
  mode: "mock"
```

That mode creates lightweight placeholder outputs for antiSMASH, GECCO, DeepBGC, and ARTS so the workflow can be tested end-to-end from the provided demo Bakta sample.

Switch to real tool execution with:

```yaml
execution:
  mode: "local"
```

Use Docker-backed execution with:

```yaml
execution:
  mode: "docker"
  work_root: "/work"
  host_db_root: "db"
  container_db_root: "/db"
```

## Databases and reference data

For real functional runs, these external assets are needed:

- `antiSMASH`: database directory under `db/antismash`
- `DeepBGC`: downloaded models and Pfam resources under `db/deepbgc`
- `ARTS`: the Actinobacteria reference under `db/arts/actinobacteria`, prepared automatically
- `eggNOG-mapper`: annotation and DIAMOND databases under `db/eggnog`
- `GECCO`: no separate database directory is assumed by this workflow

Suggested local layout:

```text
db/
  antismash/
  deepbgc/
  arts/
    actinobacteria/
  eggnog/
```

In `docker` mode the workflow mounts:

- project root to `/work`
- `db/antismash` to `/db/antismash`
- `db/deepbgc` to `/db/deepbgc`
- `db/arts/{taxon}` to `/db/arts_ref`
- `db/eggnog` to `/db/eggnog`

You would then run:

```bash
snakemake --cores 4
```

## Database bootstrap

The repository now includes a database scaffold and manifest:

- `db/manifest.yaml`
- `scripts/prepare_databases.sh`
- `scripts/fetch_databases.sh`
- `scripts/build_tool_images.sh`
- `scripts/fetch_arts.sh`
- `scripts/check_databases.py`
- `docs/database_setup.md`

Prepare the directory layout with:

```bash
bash scripts/prepare_databases.sh
```

Attempt automated downloads where supported:

```bash
bash scripts/fetch_databases.sh
```

The default `ARTS_REFERENCE=actinobacteria` is installed from the pinned ARTS
source bundled in the application image. A filesystem lock prevents concurrent
containers from preparing the same database tree. Completed resources are
validated on every start and skipped; an interrupted download is retried on the
next start.

Build local runtime images for GECCO and DeepBGC:

```bash
bash scripts/build_tool_images.sh
```

Fetch eggNOG-mapper databases where supported:

```bash
bash scripts/fetch_eggnog_db.sh
```

Then validate whether the required resources are actually present:

```bash
python3 scripts/check_databases.py
```

Check whether each tool can run locally on the current host, or whether the
workflow will rely on already existing raw outputs:

```bash
python3 scripts/check_tool_runtimes.py
```

Expected scaffold:

```text
db/
  manifest.yaml
  antismash/
  deepbgc/
  arts/
    actinobacteria/
  eggnog/
```

## Notes

- The parser layer is intentionally defensive because these tools emit different schemas across versions.
- `ARTS` installation is often the trickiest part; the current environment file is a scaffold and may need project-specific pinning.
- The workflow is designed so raw tool outputs stay isolated and all cross-tool logic happens in the summary layer.
