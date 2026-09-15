<p align="center">
  <img src="logo.jpg" alt="BGC-XPLORER logo" width="420">
</p>

<h1 align="center">BGC-XPLORER</h1>

<p align="center">
  A reproducible web application for discovering, comparing and interpreting biosynthetic gene clusters in bacterial genomes.
</p>

BGC-XPLORER accepts a bacterial genome in FASTA format and runs an integrated natural-product discovery workflow. It annotates the genome, combines predictions from several BGC callers, adds functional and resistance evidence, prioritizes promising regions and produces an interactive HTML report. Optional NVIDIA-hosted AI analysis explains individual clusters from the evidence collected by the workflow.

<p align="center">
  <img src="docs/images/bgc-xplorer-report.png" alt="BGC-XPLORER report showing consensus BGC metrics, a summary, and product and class signal charts" width="1000">
</p>

## What it does

- Accepts `.fa`, `.fasta` and `.fna` genome files through a web interface.
- Annotates the genome with **Bakta**.
- Detects BGC candidates with **antiSMASH**, **GECCO** and **DeepBGC**.
- Adds **ARTS**, **eggNOG-mapper** and **dbCAN** evidence.
- Merges compatible predictions into consensus regions and ranks the candidates.
- Generates interactive gene maps, summary tables and reproducibility metadata.
- Uses a configurable NVIDIA AI model for evidence-grounded cluster interpretation.

## Quick start with Docker

You need [Docker Engine](https://docs.docker.com/engine/install/) with the Docker Compose plugin. The analysis tools are already included in the container. Reference databases are stored outside the image and downloaded automatically during the first start.

```bash
git clone https://github.com/vebaev/BGC-XPLORER.git
cd BGC-XPLORER
cp .env.example .env
```

Open `.env` and replace the placeholder with your [NVIDIA API key](https://build.nvidia.com/):

```dotenv
NVIDIA_API_KEY=your-nvidia-api-key
NVIDIA_MODEL=nvidia/nemotron-3-ultra-550b-a55b
BGC_PORT=8778
BGC_THREADS=8
BAKTA_DB_TYPE=light
```

Start the application:

```bash
docker compose up -d
```

Open **http://localhost:8778**, upload a bacterial genome FASTA file and start the analysis. Follow startup or analysis logs with:

```bash
docker compose logs -f
```

Stop the application with `docker compose down`. Your inputs, databases and results remain in the mounted `data/`, `db/` and `results/` directories.

> The first start can take considerable time and disk space because the biological reference databases must be downloaded. Completed databases are reused and are not automatically upgraded on later starts.

## Configuration

The main settings live in `.env`:

| Variable | Default | Purpose |
| --- | --- | --- |
| `NVIDIA_API_KEY` | required | API key used by the optional AI interpretation service. |
| `NVIDIA_MODEL` | Set in `.env` | NVIDIA-hosted model used in reports; the example above uses Nemotron 3 Ultra. |
| `BGC_PORT` | `8778` | Port exposed on the host. |
| `BGC_THREADS` | `8` in `.env.example` | CPU limit shared by Snakemake and supported tools. |
| `BAKTA_DB_TYPE` | `light` | Bakta database: `light` for a smaller download or `full` for maximum coverage. |
| `AUTO_PREPARE_DATABASES` | `true` | Downloads missing databases before starting the app. |
| `ARTS_REFERENCE` | `actinobacteria` | Taxon-specific ARTS reference set. |
| `BGC_IMAGE` | `vebaev/bgc-xplorer` | GHCR image owner and name. |
| `BGC_VERSION` | `latest` | Container image tag to run. |

To use another port, model or CPU limit, edit `.env` and recreate the service:

```bash
docker compose up -d --force-recreate
```

## Reference databases

The container keeps reference data in the mounted `./db/` directory, outside the image. With `AUTO_PREPARE_DATABASES=true`, startup checks the selected **Bakta** database and the required **antiSMASH**, **DeepBGC**, **ARTS** (Actinobacteria) and **eggNOG-mapper** resources. Missing or incomplete databases are downloaded and validated before the web app starts. You can follow progress with `docker compose logs -f`.

Choose `BAKTA_DB_TYPE=light` for the smaller Bakta database or `BAKTA_DB_TYPE=full` for the full database. Both can coexist under `./db/bakta/`. The first download may take time and require substantial disk space. Later starts reuse valid databases without checking for new versions; an interrupted download is retried on the next start. Removing a specific database directory explicitly requests a fresh installation.

For the directory layout and manual preparation commands, see [Database setup](docs/database_setup.md).

## Results and reproducibility

Each analysis is stored under `results/<sample>/`. The main report is:

```text
results/<sample>/report/<sample>.html
```

Every completed sample also includes `provenance.json`, which records the application version, commit, tool versions, selected model and ARTS reference, effective configuration, and checksums for inputs, result tables and databases. API keys are excluded.

The workflow produces computational hypotheses. BGC classes, biological activities and AI interpretations require expert review and experimental validation.

## Local development

Build and run directly from the checked-out source:

```bash
cp .env.example .env
# Add NVIDIA_API_KEY to .env first.
docker compose -f docker-compose.yml -f docker-compose.dev.yml up -d --build
```

Run the unit tests with:

```bash
python3 -m unittest discover -s tests
```

Detailed setup and implementation notes are available in [Database setup](docs/database_setup.md) and [AI cluster analysis](docs/ai_cluster_analysis.md).

## Citation

If you use BGC-XPLORER in research, cite the associated publication when its citation becomes available and report the release tag or commit used for the analysis.

## License

The project license and the licenses of bundled third-party tools and downloaded databases govern redistribution and use. Review them before publishing or redistributing a derived container image.
