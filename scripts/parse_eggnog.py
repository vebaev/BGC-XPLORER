from pathlib import Path

import pandas as pd

from common import write_tsv


DEFAULT_COLUMNS = [
    "query",
    "seed_ortholog",
    "evalue",
    "score",
    "eggNOG_OGs",
    "max_annot_lvl",
    "COG_category",
    "Description",
    "Preferred_name",
    "GOs",
    "EC",
    "KEGG_ko",
    "PFAMs",
]


def find_annotations_file(outdir, sample):
    candidates = [
        outdir / "{sample}.emapper.annotations".format(sample=sample),
        outdir / "emapper.annotations",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    matches = sorted(outdir.glob("*.emapper.annotations"))
    return matches[0] if matches else None


def load_emapper_annotations(path):
    header = None
    rows = []
    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#query"):
                header = line.lstrip("#").split("\t")
                continue
            if line.startswith("#"):
                continue
            rows.append(line.split("\t"))
    if header is None:
        header = list(DEFAULT_COLUMNS)
    if not rows:
        return pd.DataFrame(columns=header)
    width = len(header)
    normalized = []
    for row in rows:
        if len(row) < width:
            row = row + [""] * (width - len(row))
        elif len(row) > width:
            row = row[:width]
        normalized.append(row)
    return pd.DataFrame(normalized, columns=header)


sample = snakemake.wildcards.sample
outdir = Path(snakemake.input.done).parent
annotations_path = find_annotations_file(outdir, sample)
bakta = pd.read_csv(snakemake.input.bakta, sep="\t")

if annotations_path is None:
    merged = pd.DataFrame(columns=[
        "sample", "locus_tag", "contig", "start", "end", "strand", "gene",
        "bakta_product", "seed_ortholog", "evalue", "score", "eggNOG_OGs",
        "max_annot_lvl", "COG_category", "eggnog_description", "preferred_name",
        "gos", "ec", "kegg_ko", "pfams", "source_file",
    ])
else:
    annotations = load_emapper_annotations(annotations_path)
    if "query" not in annotations.columns:
        annotations["query"] = ""
    annotations = annotations.rename(columns={
        "query": "locus_tag",
        "Description": "eggnog_description",
        "Preferred_name": "preferred_name",
        "GOs": "gos",
        "EC": "ec",
        "PFAMs": "pfams",
    })
    if "product" in bakta.columns:
        bakta = bakta.rename(columns={"product": "bakta_product"})
    merged = bakta.merge(annotations, on="locus_tag", how="left")
    merged["sample"] = sample
    merged["source_file"] = str(annotations_path)
    keep = [
        "sample", "locus_tag", "contig", "start", "end", "strand", "gene",
        "bakta_product", "seed_ortholog", "evalue", "score", "eggNOG_OGs",
        "max_annot_lvl", "COG_category", "eggnog_description", "preferred_name",
        "gos", "ec", "KEGG_ko", "pfams", "source_file",
    ]
    for column in keep:
        if column not in merged.columns:
            merged[column] = ""
    merged = merged[keep]

write_tsv(merged, snakemake.output[0])
