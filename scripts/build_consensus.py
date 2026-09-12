import pandas as pd

from common import add_consensus_label, load_table_if_exists, write_tsv


CORE_GENE_PATTERNS = [
    "nrps", "non-ribosomal", "polyketide", "pks", "terpene", "lanthipeptide",
    "lassopeptide", "thiopeptide", "ripp", "bacteriocin", "siderophore",
    "synthetase", "cyclase", "transferase", "methyltransferase", "oxidase",
    "dehydrogenase", "monooxygenase", "dioxygenase", "halogenase", "thioesterase",
    "acyl carrier", "peptide synthetase", "prenyl", "glycosyltransferase",
]

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


def classify_overlap_relationship(left_row, right_row):
    left_rec, right_rec, overlap_bp = reciprocal_overlap(
        left_row["start_num"], left_row["end_num"], right_row["start_num"], right_row["end_num"]
    )
    if overlap_bp <= 0:
        return ""
    if left_rec >= 0.30 and right_rec >= 0.30:
        return "reciprocal_overlap"
    if (
        (left_row["start_num"] <= right_row["start_num"] and left_row["end_num"] >= right_row["end_num"])
        or
        (right_row["start_num"] <= left_row["start_num"] and right_row["end_num"] >= left_row["end_num"])
    ):
        return "nested"
    return ""


def is_biosynthetic_core_gene(row):
    text = " ".join(
        [
            str(row.get("gene", "")).strip(),
            str(row.get("product", "")).strip(),
            str(row.get("dbxrefs", "")).strip(),
        ]
    ).lower()
    if not text:
        return False
    return any(pattern in text for pattern in CORE_GENE_PATTERNS)


def overlapping_gene_set(row, genes_df):
    if genes_df.empty:
        return set()
    local = genes_df[
        (genes_df["contig"].astype(str) == str(row["contig"]))
        & (genes_df["start_num"] <= float(row["end_num"]))
        & (genes_df["end_num"] >= float(row["start_num"]))
    ]
    if local.empty:
        return set()
    return set(local["locus_tag"].fillna("").astype(str).tolist())


def should_merge(left_row, right_row):
    if left_row["tool"] == right_row["tool"]:
        return False, ""
    if left_row["contig"] != right_row["contig"]:
        return False, ""
    if not has_valid_interval(left_row) or not has_valid_interval(right_row):
        return False, ""
    left_core = set(left_row.get("core_genes", set()) or set())
    right_core = set(right_row.get("core_genes", set()) or set())
    if left_core and right_core and (left_core & right_core):
        return True, "shared_core_gene"
    relationship = classify_overlap_relationship(left_row, right_row)
    return bool(relationship), relationship


def build_tool_core_support(subset):
    support = []
    for tool_name in ["antismash", "gecco", "deepbgc"]:
        local = subset[subset["tool"] == tool_name]
        if local.empty:
            continue
        support.append("{tool}:{count}".format(tool=tool_name, count=len(local)))
    return "; ".join(support)


def summarize_shared_core_genes(subset):
    if "core_genes" not in subset.columns:
        return ""
    tool_sets = []
    for _, local in subset.groupby("tool"):
        if local.empty:
            continue
        genes = set(local.iloc[0].get("core_genes", set()) or set())
        if genes:
            tool_sets.append(genes)
    if not tool_sets:
        return ""
    shared = set.intersection(*tool_sets) if len(tool_sets) > 1 else set(tool_sets[0])
    if shared:
        return ",".join(sorted(shared))
    per_tool = []
    for tool_name, local in subset.groupby("tool"):
        genes = sorted(set(local.iloc[0].get("core_genes", set()) or set()))
        if genes:
            per_tool.append("{tool}:{genes}".format(tool=tool_name, genes=",".join(genes[:6])))
    return "; ".join(per_tool)


def score_boundary_confidence(num_supporting_tools, overlap_relationships, edge_truncated):
    relationships = set(overlap_relationships)
    if num_supporting_tools >= 3 and "reciprocal_overlap" in relationships and not edge_truncated:
        return "high"
    if num_supporting_tools >= 2 and not edge_truncated:
        return "medium"
    if edge_truncated and num_supporting_tools >= 2:
        return "medium_edge"
    return "low"


def summarize_arts_hits(hit_rows):
    if hit_rows.empty:
        return ""
    evidence_labels = {
        "known_hit": "known resistance hit",
        "duf_hit": "DUF-associated hit",
        "core_gene": "core-gene signal",
        "duplication_signal": "duplication signal",
    }
    evidence_counts = []
    for evidence, group in hit_rows.groupby("evidence", dropna=False):
        label = evidence_labels.get(str(evidence).strip(), str(evidence).strip() or "ARTS signal")
        features = unique_preserve_order(group["feature"].tolist())[:3]
        feature_text = ", ".join(features) if features else "unspecified features"
        evidence_counts.append(
            "{count} {label} ({features})".format(
                count=len(group),
                label=label,
                features=feature_text,
            )
        )
    return "ARTS overlap: " + "; ".join(evidence_counts)


def build_why_prioritized(row):
    reasons = []
    support_count = int(row.get("num_supporting_tools", 0) or 0)
    arts_hits = int(row.get("arts_hits", 0) or 0)
    support_tools = unique_preserve_order(str(row.get("supporting_tools", "")).split(","))
    dereplication_status = str(row.get("dereplication_status", "")).strip()
    novelty_score = row.get("novelty_score", "")
    edge_truncated = bool(row.get("edge_truncated", False))
    products_text = str(row.get("products", "")).lower()

    if support_count >= 3:
        reasons.append("{0} tools agree".format(support_count))
    elif support_count == 2:
        reasons.append("2 tools agree")
    else:
        reasons.append("single-tool only")

    if arts_hits > 0:
        reasons.append("ARTS support ({0} hits)".format(arts_hits))
    else:
        reasons.append("no ARTS support")

    if "shared_core_gene" in str(row.get("overlap_relationship", "")):
        reasons.append("shared core-gene support")

    if dereplication_status == "novel_candidate":
        reasons.append("no close MIBiG match")
    elif dereplication_status == "divergent":
        reasons.append("divergent from known MIBiG cluster")
    elif dereplication_status == "related":
        reasons.append("related to known MIBiG cluster")
    elif dereplication_status == "known-like":
        reasons.append("known-like MIBiG match")
    elif novelty_score != "" and not pd.isna(novelty_score):
        reasons.append("MIBiG status unknown")

    if "antibacterial" in products_text:
        reasons.append("antibacterial signal")
    if "cytotoxic" in products_text:
        reasons.append("cytotoxic signal")
    if edge_truncated:
        reasons.append("contig-edge region")
    if support_tools:
        reasons.append("callers: {0}".format(", ".join(support_tools)))

    return " + ".join(reasons)


def build_why_not_prioritized(row):
    reasons = []
    support_count = int(row.get("num_supporting_tools", 0) or 0)
    arts_hits = int(row.get("arts_hits", 0) or 0)
    edge_truncated = bool(row.get("edge_truncated", False))
    dereplication_status = str(row.get("dereplication_status", "")).strip()
    boundary_confidence = str(row.get("boundary_confidence", "")).strip()
    priority_class = str(row.get("priority_class", "")).strip()

    if support_count <= 1:
        reasons.append("single-tool only")
    elif support_count == 2:
        reasons.append("only 2 tools support this locus")

    if arts_hits == 0:
        reasons.append("no ARTS support")
    elif support_count <= 1 and arts_hits > 0:
        reasons.append("ARTS support without cross-tool confirmation")

    if boundary_confidence in ["low", "medium_edge"]:
        reasons.append("uncertain region boundaries")

    if edge_truncated:
        reasons.append("near contig edge; possible partial BGC")

    if dereplication_status == "known-like":
        reasons.append("strong similarity to known MIBiG cluster")
    elif dereplication_status == "related":
        reasons.append("close to a known MIBiG family")

    if priority_class == "low" and not reasons:
        reasons.append("limited supporting evidence")

    return "; ".join(unique_preserve_order(reasons))


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
    edge_threshold_bp = int(snakemake.config["consensus"].get("edge_threshold_bp", 5000))
    contig_length = int(contig_lengths.get(row["contig_id"], 0) or 0)
    start = int(row["start"]) if row["start"] != "" and not pd.isna(row["start"]) else 0
    end = int(row["end"]) if row["end"] != "" and not pd.isna(row["end"]) else 0
    left_distance = max(start - 1, 0) if start else ""
    right_distance = max(contig_length - end, 0) if contig_length and end else ""
    edge_truncated = bool(
        contig_length
        and (
            (start and start <= edge_threshold_bp)
            or (end and (contig_length - end) <= edge_threshold_bp)
        )
    )
    row["contig_length"] = contig_length
    row["distance_to_left_edge"] = left_distance
    row["distance_to_right_edge"] = right_distance
    row["edge_truncated"] = edge_truncated
    row["possible_partial_BGC"] = edge_truncated
    return row


antismash = load_table_if_exists(snakemake.input.antismash, [
    "sample", "tool", "contig", "start", "end", "strand", "bgc_id", "bgc_type", "product", "score", "confidence", "source_file"
])
gecco = load_table_if_exists(snakemake.input.gecco, antismash.columns.tolist())
deepbgc = load_table_if_exists(snakemake.input.deepbgc, antismash.columns.tolist())
arts = load_table_if_exists(snakemake.input.arts, [
    "sample", "tool", "contig", "start", "end", "feature", "score", "evidence", "source_file"
])
bakta = pd.read_csv(snakemake.input.bakta, sep="\t") if snakemake.input.bakta else pd.DataFrame()

contig_lengths = parse_fasta_lengths(snakemake.input.fna)

if not bakta.empty:
    for column in ["contig", "type", "start", "end", "strand", "locus_tag", "gene", "product", "dbxrefs"]:
        if column not in bakta.columns:
            bakta[column] = ""
    bakta["start_num"] = to_num(bakta["start"])
    bakta["end_num"] = to_num(bakta["end"])
    bakta = bakta[bakta["type"].astype(str) == "cds"].copy()
    bakta["is_core_gene"] = bakta.apply(is_biosynthetic_core_gene, axis=1)
    core_genes = bakta[bakta["is_core_gene"]].copy()
else:
    core_genes = pd.DataFrame(columns=["contig", "start_num", "end_num", "locus_tag"])

bgcs = pd.concat([antismash, gecco, deepbgc], ignore_index=True)
if not bgcs.empty:
    bgcs["start_num"] = to_num(bgcs["start"])
    bgcs["end_num"] = to_num(bgcs["end"])
    bgcs["length_bp"] = (bgcs["end_num"] - bgcs["start_num"] + 1).fillna(0).astype(int)
    bgcs["contig_id"] = bgcs["contig"].astype(str)
    bgcs["product_annotation"] = bgcs["product"].astype(str)
    bgcs["core_genes"] = bgcs.apply(lambda row: overlapping_gene_set(row, core_genes), axis=1)
    bgcs = bgcs.sort_values(["contig", "start_num", "end_num", "tool"], kind="mergesort").reset_index(drop=True)
else:
    bgcs["start_num"] = pd.Series(dtype=float)
    bgcs["end_num"] = pd.Series(dtype=float)
    bgcs["core_genes"] = pd.Series(dtype=object)

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
                merged, relationship = should_merge(member, other)
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
        "core_gene_support": summarize_shared_core_genes(subset) or build_tool_core_support(subset),
    }
    region = add_contig_edge_fields(region, contig_lengths)
    region["boundary_confidence"] = score_boundary_confidence(
        region["num_supporting_tools"],
        unique_preserve_order(relationships if relationships else ["single_tool"]),
        region["edge_truncated"],
    )
    consensus_rows.append(region)

consensus = pd.DataFrame(consensus_rows)
overlap = pd.DataFrame(overlap_rows)
consensus = add_consensus_label(consensus)

write_tsv(consensus, snakemake.output.consensus)
write_tsv(overlap, snakemake.output.overlap)
