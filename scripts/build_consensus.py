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


def classify_overlap_relationship(left_row, right_row, overlap_fraction=0.30):
    left_rec, right_rec, overlap_bp = reciprocal_overlap(
        left_row["start_num"], left_row["end_num"], right_row["start_num"], right_row["end_num"]
    )
    if overlap_bp <= 0:
        return ""
    if left_rec >= overlap_fraction and right_rec >= overlap_fraction:
        return "reciprocal_overlap"
    if (
        (left_row["start_num"] <= right_row["start_num"] and left_row["end_num"] >= right_row["end_num"])
        or
        (right_row["start_num"] <= left_row["start_num"] and right_row["end_num"] >= left_row["end_num"])
    ):
        return "nested"
    return ""


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


def should_merge(left_row, right_row, overlap_fraction=0.30):
    if left_row["tool"] == right_row["tool"]:
        return False, ""
    if left_row["contig"] != right_row["contig"]:
        return False, ""
    if not has_valid_interval(left_row) or not has_valid_interval(right_row):
        return False, ""
    left_core = set((left_row.get("core_gene_details", {}) or {}).keys())
    right_core = set((right_row.get("core_gene_details", {}) or {}).keys())
    if left_core and right_core and (left_core & right_core):
        return True, "shared_core_gene"
    relationship = classify_overlap_relationship(left_row, right_row, overlap_fraction)
    return bool(relationship), relationship


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
merge_threshold = float(snakemake.config["consensus"].get("overlap_fraction", 0.30))

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

consensus_rows = []
overlap_rows = []
group_id = 0
visited = set()

for idx, row in bgcs.iterrows():
    if idx in visited:
        continue
    group_id += 1
    members = [idx]
    visited.add(idx)
    relationships = []

    expanded = True
    while expanded:
        expanded = False
        for jdx, other in bgcs.iterrows():
            if jdx in visited:
                continue
            merge_this = False
            local_relationships = []
            for member_idx in members:
                member = bgcs.loc[member_idx]
                merged, relationship = should_merge(member, other, merge_threshold)
                if merged:
                    merge_this = True
                    local_relationships.append((member, other, relationship))
            if merge_this:
                members.append(jdx)
                visited.add(jdx)
                expanded = True
                for member, other, relationship in local_relationships:
                    relationships.append(relationship)
                    _, _, overlap_bp = reciprocal_overlap(
                        member["start_num"], member["end_num"], other["start_num"], other["end_num"]
                    )
                    shared_core = sorted(
                        set((member.get("core_gene_details", {}) or {}).keys())
                        & set((other.get("core_gene_details", {}) or {}).keys())
                    )
                    shared_evidence = []
                    for locus in shared_core:
                        for source_row in (member, other):
                            record = (source_row.get("core_gene_details", {}) or {}).get(locus, {})
                            shared_evidence.append("{0}:{1}".format(
                                record.get("source", source_row.get("tool", "")),
                                record.get("evidence", ""),
                            ))
                    overlap_rows.append({
                        "sample": member["sample"],
                        "group_id": group_id,
                        "tool_a": member["tool"],
                        "bgc_id_a": member["bgc_id"],
                        "tool_b": other["tool"],
                        "bgc_id_b": other["bgc_id"],
                        "contig": member["contig"],
                        "overlap_bp": overlap_bp,
                        "overlap_relationship": relationship,
                        "shared_core_genes": ",".join(shared_core),
                        "shared_core_evidence": "; ".join(unique_preserve_order(shared_evidence)),
                    })

    subset = bgcs.loc[members].copy()
    tools = sorted(set(subset["tool"]))
    contig = choose_first_nonempty(subset["contig"])
    start_val = int(subset["start_num"].dropna().min()) if not subset["start_num"].dropna().empty else ""
    end_val = int(subset["end_num"].dropna().max()) if not subset["end_num"].dropna().empty else ""
    region = {
        "sample": row["sample"],
        "consensus_id": "{sample}_consensus_{group_id}".format(sample=row["sample"], group_id=group_id),
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
    "contig", "overlap_bp", "overlap_relationship", "shared_core_genes",
    "shared_core_evidence",
])
consensus = add_consensus_label(consensus)

write_tsv(consensus, snakemake.output.consensus)
write_tsv(overlap, snakemake.output.overlap)
