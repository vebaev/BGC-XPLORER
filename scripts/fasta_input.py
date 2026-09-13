"""Pure helpers for the FASTA-only upload flow."""

import re
from pathlib import Path


FASTA_EXTENSIONS = (".fa", ".fasta", ".fna")


def fasta_suffix(filename):
    suffix = Path(filename).suffix.lower()
    return suffix if suffix in FASTA_EXTENSIONS else None


def normalize_sample_name(value):
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())
    name = name.strip("._-")
    if not name:
        raise ValueError("Sample name must contain a letter or number")
    return name

