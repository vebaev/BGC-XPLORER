import pandas as pd

from common import BGC_COLUMNS, add_consensus_label, load_table_if_exists, write_tsv
from core_gene_evidence import bakta_core_match, parse_records


def to_num(series):
    return pd.to_numeric(series, errors="coerce")


def interval_overlap(a_start, a_end, b_start, b_end):
    left = max(a_start, b_start)
    right = min(a_end, b_end)
    return max(0, right - left + 1)


def interval_length(start, end):
    if pd.isna(start) or pd.isna(end):
        return 0
    return int(end) - int(start) + 1


def reciprocal_overlap(a_start, a_end, b_start, b_end):
    overlap = interval_overlap(a_start, a_end, b_start, b_end)
    if overlap <= 0:
        return 0.0, 0.0, 0
    a_len = max(interval_length(a_start, a_end), 1)
    b_len = max(interval_length(b_start, b_end), 1)
    return overlap / float(a_len), overlap / float(b_len), overlap


def has_valid_interval(row):
    return not pd.isna(row["start_num"]) and not pd.isna(row["end_num"])


def choose_first_nonempty(series):
    values = series.fillna("").astype(str)
    values = values[values.str.strip() != ""]
    if values.empty:
        return ""
    return values.iloc[0]


def unique_preserve_order(values):
    seen = []
    for value in values:
        text = str(value).strip()
        if not text or text.lower() == "nan":
            continue
        if text not in seen:
            seen.append(text)
    return seen


def unique_join(values, sep=","):
    return sep.join(unique_preserve_order(values))


def parse_fasta_lengths(path):
    lengths = {}
    current = None
    size = 0
    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    lengths[current] = size
                current = line[1:].split()[0]
                size = 0
            else:
                size += len(line)
    if current is not None:
        lengths[current] = size
    return lengths


def containment(left_row, right_row):
    """Fraction of the shorter region that the longer one covers.

    Callers disagree on boundaries mostly because they disagree on how far to
    extend past the core genes, so the length ratio - not the biology - drives
    reciprocal overlap. Containment is stable under that asymmetry: a short
    prediction sitting inside a long one scores 1.0 either way round.
    """
    left_frac, right_frac, overlap_bp = reciprocal_overlap(
        left_row["start_num"], left_row["end_num"], right_row["start_num"], right_row["end_num"]
    )
    if overlap_bp <= 0:
        return 0.0, 0
    return max(left_frac, right_frac), overlap_bp


def overlapping_genes(row, genes_df):
    if genes_df.empty:
        return genes_df.iloc[0:0]
    local = genes_df[
        (genes_df["contig"].astype(str) == str(row["contig"]))
        & (genes_df["start_num"] <= float(row["end_num"]))
        & (genes_df["end_num"] >= float(row["start_num"]))
    ]
    if local.empty:
        return genes_df.iloc[0:0]
    return local


def resolve_core_gene_records(row, bakta):
    """Resolve predictor evidence to Bakta locus tags, falling back per region."""
    resolved = {}
    records = parse_records(row.get("core_gene_records", ""))
    for record in records:
        locus = str(record.get("identifier", "")).strip()
        if locus and locus in set(bakta["locus_tag"].fillna("").astype(str)):
            resolved[locus] = record
            continue
        if record.get("start") is None or record.get("end") is None:
            continue
        candidates = bakta[
            (bakta["contig"].astype(str) == str(row["contig"]))
            & (bakta["start_num"] <= float(record["end"]))
            & (bakta["end_num"] >= float(record["start"]))
        ].copy()
        if candidates.empty:
            continue
        candidates["mapping_overlap"] = candidates.apply(
            lambda gene: interval_overlap(
                gene["start_num"], gene["end_num"], record["start"], record["end"]
            ), axis=1,
        )
        best = candidates.sort_values(
            ["mapping_overlap", "locus_tag"], ascending=[False, True]
        ).iloc[0]
        locus = str(best.get("locus_tag", "")).strip()
        if locus:
            resolved[locus] = record
    if resolved:
        return resolved

    # The predictor supplied no resolvable core role for this region.
    for _, gene in overlapping_genes(row, bakta).iterrows():
        evidence = bakta_core_match(gene)
        locus = str(gene.get("locus_tag", "")).strip()
        if evidence and locus:
            resolved[locus] = {
                "identifier": locus,
                "start": int(gene["start_num"]),
                "end": int(gene["end_num"]),
                "source": "bakta_fallback",
                "evidence": evidence,
            }
    return resolved


# Rule-based callers report full regions and explicit core-gene roles, so they
# make the most reliable anchors; the ML callers follow.
TOOL_PRIORITY = {"antismash": 0, "gecco": 1, "deepbgc": 2}


def anchor_sort_key(row):
    """Order candidates so the same input always produces the same groups."""
    start = row["start_num"]
    return (
        TOOL_PRIORITY.get(str(row["tool"]).strip().lower(), len(TOOL_PRIORITY)),
        -int(row.get("length_bp", 0) or 0),
        str(row["contig"]),
        float(start) if not pd.isna(start) else 0.0,
        str(row["bgc_id"]),
    )


def boundary_sort_key(row):
    """Pick the member that defines the reported boundaries: the longest one.

    The anchor is chosen for caller reliability, not for extent, so a short
    antiSMASH or GECCO region can anchor a group that also holds a much longer
    prediction. Taking the anchor's coordinates would discard that extent; taking
    the longest member keeps boundaries a caller actually reported without the
    invented span that min/max across the group produces.
    """
    return (
        -int(row.get("length_bp", 0) or 0),
        TOOL_PRIORITY.get(str(row["tool"]).strip().lower(), len(TOOL_PRIORITY)),
        str(row["bgc_id"]),
    )


def should_merge(anchor_row, other_row, min_containment=0.80):
    """Decide whether a candidate describes the same locus as the group anchor.

    Every member is judged against the anchor alone. Shared core genes are
    reported as evidence but never merge on their own: they used to join
    regions that do not overlap at all.
    """
    if anchor_row["contig"] != other_row["contig"]:
        return False, "", 0.0, 0
    if not has_valid_interval(anchor_row) or not has_valid_interval(other_row):
        return False, "", 0.0, 0
    fraction, overlap_bp = containment(anchor_row, other_row)
    if fraction < min_containment:
        return False, "", fraction, overlap_bp
    relationship = "nested" if fraction >= 0.999 else "containment"
    return True, relationship, fraction, overlap_bp


def summarize_shared_core_genes(subset):
    if "core_gene_details" not in subset.columns:
        return ""
    per_tool = []
    for tool_name, local in subset.groupby("tool"):
        genes = sorted({
            locus
            for details in local["core_gene_details"]
            for locus in (details or {}).keys()
        })
        if genes:
            per_tool.append("{tool}:{genes}".format(tool=tool_name, genes=",".join(genes[:6])))
    return "; ".join(per_tool)


def summarize_core_evidence(subset):
    summaries = []
    for _, row in subset.iterrows():
        for locus, record in (row.get("core_gene_details", {}) or {}).items():
            summaries.append("{locus} [{source}: {evidence}]".format(
                locus=locus,
                source=record.get("source", row.get("tool", "")),
                evidence=record.get("evidence", ""),
            ))
    return "; ".join(unique_preserve_order(summaries))


def biological_interpretation(values):
    text = ",".join(str(value) for value in values if str(value).strip()).lower()
    findings = []
    if any(term in text for term in ["nrps", "nrp", "nrps-like"]):
        findings.append("NRPS/non-ribosomal peptide candidate; often linked to siderophores, antibiotics, lipopeptides, or other peptide natural products.")
    if any(term in text for term in ["t1pks", "pks", "polyketide"]):
        findings.append("PKS/polyketide candidate; commonly associated with antibiotics, pigments, toxins, or signaling metabolites.")
    if "terpene" in text:
        findings.append("Terpene candidate; may encode volatile, pigment, membrane-active, or signaling terpenoid metabolites.")
    if any(term in text for term in ["ripp", "azole-containing-ripp", "rre-containing"]):
        findings.append("RiPP candidate; may encode ribosomally synthesized and post-translationally modified peptides such as antimicrobial peptides.")
    if "siderophore" in text:
        findings.append("Siderophore candidate; likely involved in iron acquisition and metal competition.")
    if "saccharide" in text:
        findings.append("Saccharide/glycan-related candidate; may encode exopolysaccharide, glycosylated metabolite, or cell-surface carbohydrate biosynthesis.")
    if "betalactone" in text or "beta-lactone" in text:
        findings.append("Beta-lactone candidate; this class can include enzyme inhibitors and other reactive natural products.")
    if "ectoine" in text:
        findings.append("Ectoine/osmoprotectant candidate; likely involved in osmotic or stress adaptation.")
    if "lanthipeptide" in text:
        findings.append("Lanthipeptide candidate; a RiPP subtype often associated with antimicrobial activity.")
    if "antibacterial" in text:
        findings.append("DeepBGC predicts antibacterial activity for at least one overlapping candidate.")
    if "cytotoxic" in text:
        findings.append("DeepBGC predicts cytotoxic activity for at least one overlapping candidate.")
    if not findings:
        findings.append("Unclassified BGC-like candidate; biological function needs manual inspection of genes, domains, and nearest known clusters.")
    return " ".join(dict.fromkeys(findings))


def add_contig_edge_fields(row, contig_lengths):
    contig_length = int(contig_lengths.get(row["contig_id"], 0) or 0)
    start = int(row["start"]) if row["start"] != "" and not pd.isna(row["start"]) else 0
    end = int(row["end"]) if row["end"] != "" and not pd.isna(row["end"]) else 0
    left_distance = max(start - 1, 0) if start else ""
    right_distance = max(contig_length - end, 0) if contig_length and end else ""
    row["contig_length"] = contig_length
    row["distance_to_left_edge"] = left_distance
    row["distance_to_right_edge"] = right_distance
    return row


antismash = load_table_if_exists(snakemake.input.antismash, BGC_COLUMNS)
gecco = load_table_if_exists(snakemake.input.gecco, antismash.columns.tolist())
deepbgc = load_table_if_exists(snakemake.input.deepbgc, antismash.columns.tolist())
arts = load_table_if_exists(snakemake.input.arts, [
    "sample", "tool", "contig", "start", "end", "feature", "score", "evidence", "source_file"
])
bakta = pd.read_csv(snakemake.input.bakta, sep="\t") if snakemake.input.bakta else pd.DataFrame()

contig_lengths = parse_fasta_lengths(snakemake.input.fna)
merge_threshold = float(snakemake.config["consensus"].get("min_containment", 0.80))

if not bakta.empty:
    for column in ["contig", "type", "start", "end", "strand", "locus_tag", "gene", "product", "dbxrefs"]:
        if column not in bakta.columns:
            bakta[column] = ""
    bakta["start_num"] = to_num(bakta["start"])
    bakta["end_num"] = to_num(bakta["end"])
    bakta = bakta[bakta["type"].astype(str) == "cds"].copy()
else:
    bakta = pd.DataFrame(columns=["contig", "start_num", "end_num", "locus_tag", "gene", "product", "dbxrefs"])

bgcs = pd.concat([antismash, gecco, deepbgc], ignore_index=True)
if not bgcs.empty:
    bgcs["start_num"] = to_num(bgcs["start"])
    bgcs["end_num"] = to_num(bgcs["end"])
    bgcs["length_bp"] = (bgcs["end_num"] - bgcs["start_num"] + 1).fillna(0).astype(int)
    bgcs["contig_id"] = bgcs["contig"].astype(str)
    bgcs["product_annotation"] = bgcs["product"].astype(str)
    bgcs["core_gene_details"] = bgcs.apply(lambda row: resolve_core_gene_records(row, bakta), axis=1)
    bgcs = bgcs.sort_values(["contig", "start_num", "end_num", "tool"], kind="mergesort").reset_index(drop=True)
else:
    bgcs["start_num"] = pd.Series(dtype=float)
    bgcs["end_num"] = pd.Series(dtype=float)
    bgcs["core_gene_details"] = pd.Series(dtype=object)

groups = []
claimed = set()
contig_index = (
    {str(contig): list(frame.index) for contig, frame in bgcs.groupby("contig_id", sort=False)}
    if not bgcs.empty else {}
)
anchor_order = sorted(bgcs.index, key=lambda index: anchor_sort_key(bgcs.loc[index])) if not bgcs.empty else []

# Anchor linkage: the strongest unclaimed prediction defines the locus and every
# other member is judged against it. Single linkage used to chain A-B-C together
# through a bridging region, stretching boundaries past anything a caller called.
for anchor_idx in anchor_order:
    if anchor_idx in claimed:
        continue
    anchor = bgcs.loc[anchor_idx]
    claimed.add(anchor_idx)
    members = [anchor_idx]
    edges = []
    for other_idx in contig_index.get(str(anchor["contig_id"]), []):
        if other_idx in claimed:
            continue
        merged, relationship, fraction, overlap_bp = should_merge(
            anchor, bgcs.loc[other_idx], merge_threshold
        )
        if not merged:
            continue
        claimed.add(other_idx)
        members.append(other_idx)
        edges.append((other_idx, relationship, fraction, overlap_bp))
    groups.append((anchor_idx, members, edges))


def group_position(group):
    members = bgcs.loc[group[1]]
    starts = members["start_num"].dropna()
    return (
        str(members.iloc[0]["contig_id"]),
        float(starts.min()) if not starts.empty else 0.0,
    )


groups.sort(key=group_position)

consensus_rows = []
overlap_rows = []

for group_id, (anchor_idx, members, edges) in enumerate(groups, start=1):
    anchor = bgcs.loc[anchor_idx]
    for other_idx, relationship, fraction, overlap_bp in edges:
        other = bgcs.loc[other_idx]
        shared_core = sorted(
            set((anchor.get("core_gene_details", {}) or {}).keys())
            & set((other.get("core_gene_details", {}) or {}).keys())
        )
        shared_evidence = []
        for locus in shared_core:
            for source_row in (anchor, other):
                record = (source_row.get("core_gene_details", {}) or {}).get(locus, {})
                shared_evidence.append("{0}:{1}".format(
                    record.get("source", source_row.get("tool", "")),
                    record.get("evidence", ""),
                ))
        overlap_rows.append({
            "sample": anchor["sample"],
            "group_id": group_id,
            "tool_a": anchor["tool"],
            "bgc_id_a": anchor["bgc_id"],
            "tool_b": other["tool"],
            "bgc_id_b": other["bgc_id"],
            "contig": anchor["contig"],
            "overlap_bp": overlap_bp,
            "containment": round(fraction, 4),
            "overlap_relationship": relationship,
            "shared_core_genes": ",".join(shared_core),
            "shared_core_evidence": "; ".join(unique_preserve_order(shared_evidence)),
        })

    relationships = [relationship for _, relationship, _, _ in edges]
    subset = bgcs.loc[members].copy()
    tools = sorted(set(subset["tool"]))
    contig = choose_first_nonempty(subset["contig"])
    valid = subset.dropna(subset=["start_num", "end_num"])
    union_start = int(valid["start_num"].min()) if not valid.empty else ""
    union_end = int(valid["end_num"].max()) if not valid.empty else ""
    if valid.empty:
        start_val, end_val, boundary_id = "", "", ""
    else:
        boundary = valid.loc[sorted(valid.index, key=lambda index: boundary_sort_key(valid.loc[index]))[0]]
        start_val = int(boundary["start_num"])
        end_val = int(boundary["end_num"])
        boundary_id = str(boundary["bgc_id"])
    region = {
        "sample": anchor["sample"],
        "consensus_id": "{sample}_consensus_{group_id}".format(sample=anchor["sample"], group_id=group_id),
        "contig": contig,
        "contig_id": contig,
        "start": start_val,
        "end": end_val,
        "length_bp": interval_length(start_val, end_val) if start_val != "" and end_val != "" else "",
        "support_tools": ",".join(tools),
        "supporting_tools": ",".join(tools),
        "support_count": len(tools),
        "num_supporting_tools": len(tools),
        "candidate_ids": ",".join(subset["bgc_id"].astype(str)),
        "anchor_tool": anchor["tool"],
        "anchor_bgc_id": anchor["bgc_id"],
        "boundary_bgc_id": boundary_id,
        "union_start": union_start,
        "union_end": union_end,
        "bgc_types": unique_join(subset["bgc_type"].tolist()),
        "products": unique_join(subset["product"].tolist()),
        "product_annotations": unique_join(subset["product"].tolist()),
        "biological_interpretation": biological_interpretation(
            list(subset["bgc_type"].fillna("").astype(str)) + list(subset["product"].fillna("").astype(str))
        ),
        "overlap_relationship": unique_join(relationships if relationships else ["single_tool"]),
        "core_gene_support": summarize_shared_core_genes(subset),
        "core_gene_evidence": summarize_core_evidence(subset),
    }
    region = add_contig_edge_fields(region, contig_lengths)
    consensus_rows.append(region)

consensus = pd.DataFrame(consensus_rows)
overlap = pd.DataFrame(overlap_rows, columns=[
    "sample", "group_id", "tool_a", "bgc_id_a", "tool_b", "bgc_id_b",
    "contig", "overlap_bp", "containment", "overlap_relationship",
    "shared_core_genes", "shared_core_evidence",
])
consensus = add_consensus_label(consensus)

write_tsv(consensus, snakemake.output.consensus)
write_tsv(overlap, snakemake.output.overlap)
