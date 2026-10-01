# Showcase genome: supplementary tables and figures

Scripts that build the supplementary material for the showcase genome
(*Streptomyces albidoflavus* STR-S001, GenBank CP158056, sample `Soil_1`) from a
finished BGC-XPLORER run. They read only the workflow's outputs, so every table
and figure can be rebuilt after a re-run.

| Script | Output |
|---|---|
| `loci_table.py` | one row per locus: coordinates, agreed length and gene-map window, callers and category (multi-caller with antiSMASH; GECCO + DeepBGC only; single caller), genes covered by at least two callers, biosynthetic core gene, MIBiG representative, ARTS hits, MIBiG classes and DeepBGC activity |
| `caller_upset.py` | UpSet plot of the caller combinations behind the loci; loci supported only by GECCO and DeepBGC, both Pfam-based, are set apart |

```bash
python showcase/loci_table.py --results results --sample Soil_1 --out Soil_1_loci.xlsx
python showcase/caller_upset.py --results results --sample Soil_1 --out Soil_1_upset.pdf
```

Both need pandas; the table needs openpyxl for `.xlsx`, the plot needs
matplotlib (the environment in `validation/environment.yaml` has all three).
