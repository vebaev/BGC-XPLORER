/bin/bash: warning: setlocale: LC_ALL: cannot change locale (C.UTF-8)
# BGC Discovery Workflow

Snakemake workflow that starts from a bacterial genome FASTA, annotates it with
Bakta, and runs:

- antiSMASH
- GECCO
- DeepBGC
- ARTS
- eggNOG-mapper

The workflow then normalizes the tool outputs, builds a simple consensus table, prioritizes candidate regions, and renders an HTML summary report.

## Expected input layout

The web interface accepts one uncompressed `.fa`, `.fasta`, or `.fna` file. It
stores the normalized input as:

```text
data/fasta/{sample}.fasta
```

Bakta runs as the first workflow step and writes its derived annotation files
under `data/bakta/{sample}/`.

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
BGC_THREADS=8
AUTO_PREPARE_DATABASES=true
BAKTA_DB_TYPE=light
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

`BAKTA_DB_TYPE` accepts `light` or `full`. The default `light` database needs
substantially less disk space. The selected database is stored under the
external `./db/bakta` volume and is downloaded only when no valid local copy is
present. Changing the setting to `full` downloads and then uses the full set;
existing databases are not updated automatically.

`BGC_THREADS` is a positive integer shared by Snakemake and tools that expose
parallel execution controls. It sets `snakemake --cores` and the worker/CPU
count for Bakta, antiSMASH, GECCO, eggNOG-mapper, and dbCAN. Tools without a
thread option still run under the same Snakemake core budget. The default is 4.

ARTS uses a taxon-specific reference set. If the required ARTS files are not
present under `db/arts/actinobacteria`, startup stops with the missing file
list; populate that reference set and run `docker compose up -d` again.

## Snakemake quick start

1. Put each input at `data/fasta/{sample}.fasta`
2. Edit `config/samples.tsv`
3. Set `BAKTA_DB_TYPE=light` or `BAKTA_DB_TYPE=full` and optionally `BGC_THREADS`
4. Run:

```bash
BGC_THREADS=8 snakemake --use-conda --cores 8
```

## Execution modes

The default `config/config.yaml` uses:

```yaml
execution:
  mode: "mock"
```

That mode creates lightweight placeholder outputs for the downstream BGC tools.
Bakta still requires its selected external database and a FASTA input.

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

- `Bakta`: `db/bakta/db-light` or `db/bakta/db`, selected with `BAKTA_DB_TYPE`
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

## Reproducibility records

Every completed sample produces `results/{sample}/provenance.json`. It records
the BGC-XPLORER version and commit, image reference, tool versions, selected AI
model, ARTS reference, effective configuration, SHA-256 checksums for inputs and
result tables, and SHA-256 checksums for required database files. Database hashes
are cached under the external database directory and reused while file size and
modification time remain unchanged. API keys and other secret configuration
fields are redacted. The HTML report links to the corresponding provenance file.

Successfully validated databases remain fixed. Startup does not query for or
install newer database releases; removing an individual database directory is
the explicit way to request a fresh installation.

## Notes

- The parser layer is intentionally defensive because these tools emit different schemas across versions.
- `ARTS` installation is often the trickiest part; the current environment file is a scaffold and may need project-specific pinning.
- The workflow is designed so raw tool outputs stay isolated and all cross-tool logic happens in the summary layer.
