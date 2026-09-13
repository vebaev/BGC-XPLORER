# Database Setup

This project expects external resources for Bakta, `antiSMASH`, `DeepBGC`, `ARTS`, and `eggNOG-mapper`.

## Current layout

```text
db/
  manifest.yaml
  bakta/
    db-light/  # when BAKTA_DB_TYPE=light
    db/        # when BAKTA_DB_TYPE=full
  antismash/
  deepbgc/
  arts/
    actinobacteria/
  eggnog/
```

## Bakta

Set `BAKTA_DB_TYPE=light` or `BAKTA_DB_TYPE=full` in `.env` before starting
Docker Compose. Startup validates the selected directory and runs Bakta's
official downloader only if it is absent or incomplete:

```bash
BAKTA_DB_TYPE=light BGC_DB_ROOT="$PWD/db" bash scripts/fetch_bakta_db.sh
```

The official downloader selects the newest database compatible with the pinned
Bakta runtime, verifies the archive, extracts it, and initializes AMRFinderPlus.
Once validation succeeds, later starts do not contact the release service or
check for updates. Light and full databases may coexist in the external volume.

## antiSMASH

Official documentation says local Bioconda installs should run `download-antismash-databases`, and Docker users can use the official helper script `download_antismash_databases_docker`. The full standalone Docker image includes required databases already. Source: [antiSMASH install docs](https://docs.antismash.secondarymetabolites.org/install/)

In this repository:

```bash
bash scripts/fetch_antismash_db.sh
```

That script downloads the official Docker database helper and populates `db/antismash`. It is safe to re-run when a long download was interrupted.

## DeepBGC

The project README says to download trained models and the Pfam database before use:

```bash
deepbgc download
deepbgc info
```

Source: [DeepBGC README](https://github.com/Merck/deepbgc)

In this repository:

```bash
bash scripts/fetch_deepbgc_db.sh
```

That script runs `deepbgc download` only if the `deepbgc` command is already installed locally.

For this repository, the preferred runtime is the local Docker image:

```bash
bash scripts/build_tool_images.sh
docker run --rm aktino/deepbgc:latest deepbgc --help
```

## ARTS

ARTS expects a precomputed reference directory, and this workflow maps that to:

```text
db/arts/actinobacteria/
```

Source: [ARTS README](https://github.com/ZiemertLab/ARTS)

The default reference is `actinobacteria`. The Docker image contains the reference
archives from the pinned ARTS source revision, and startup extracts them into the
external database directory automatically. For a local source checkout, run:

```bash
bash scripts/fetch_arts.sh
```

Set `ARTS_REFERENCE=actinobacteria` explicitly to override the environment. This
release rejects other values because their reference archives are not bundled or
individually versioned yet.

## eggNOG-mapper

The official repository documents installation with `pip install eggnog-mapper` and database setup with:

```bash
download_eggnog_data.py --data_dir /path/to/eggnog-data
```

The same documentation shows protein annotation with:

```bash
emapper.py -m diamond -i proteins.fa --itype proteins --data_dir /path/to/eggnog-data -o my_annotation
```

Source: [eggNOG-mapper README](https://github.com/eggnogdb/eggnog-mapper)

In this repository:

```bash
bash scripts/fetch_eggnog_db.sh
```

That script uses either a local `download_eggnog_data.py` installation or the local Docker image `aktino/eggnog:latest` if available.

## Combined bootstrap

To prepare all required databases:

```bash
bash scripts/fetch_databases.sh
```

The combined bootstrap holds a filesystem lock, validates resources before each
download, and skips valid databases. Failed or interrupted downloads remain
incomplete and are retried at the next invocation. antiSMASH intentionally uses
the official `latest` database helper.

## Validation

After resources are in place:

```bash
python3 scripts/check_databases.py
```
