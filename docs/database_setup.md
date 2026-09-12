# Database Setup

This project expects external resources for `antiSMASH`, `DeepBGC`, `ARTS`, and `eggNOG-mapper`.

## Current layout

```text
db/
  manifest.yaml
  antismash/
  deepbgc/
  arts/
    actinobacteria/
  eggnog/
```

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

This repository does not auto-download ARTS references because they are taxon-specific and need to be placed deliberately.

To fetch the official ARTS repository and unpack the bundled Actinobacteria reference zips:

```bash
bash scripts/fetch_arts.sh
```

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

If you want a best-effort combined setup:

```bash
bash scripts/fetch_databases.sh
```

## Validation

After resources are in place:

```bash
python3 scripts/check_databases.py
```
