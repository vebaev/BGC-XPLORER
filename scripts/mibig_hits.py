"""Read method-specific antiSMASH 8 MIBiG hits from its structured JSON."""

import math
import re


MIBIG_ID_RE = re.compile(r"^(BGC\d{7}\.\d+)")


def novelty_score_for_priority(value):
    try:
        score = float(value)
    except (TypeError, ValueError):
        return 0.0
    return score if math.isfinite(score) else 0.0


def mibig_method_priority(source, peptide_similarity, strong_peptide_threshold=80.0):
    if source == "comparippson_html":
        try:
            similarity = float(peptide_similarity)
        except (TypeError, ValueError):
            similarity = 0.0
        return -1 if similarity >= strong_peptide_threshold else 1
    return {"knownclusterblast": 0, "clustercompare_mibig": 2}.get(source, 3)


def extract_structured_mibig_hits(payload):
    hits = []
    for record in payload.get("records", []):
        contig_id = str(record.get("id", ""))
        modules = record.get("modules", {}) or {}
        clusterblast = modules.get("antismash.modules.clusterblast", {}) or {}
        known = clusterblast.get("knowncluster", {}) or {}
        for result in known.get("results", []):
            region_id = str(result.get("region_number", ""))
            for ranked in result.get("ranking", []):
                if len(ranked) != 2:
                    continue
                reference, score = ranked
                mibig_id = str(reference.get("accession", ""))
                if not MIBIG_ID_RE.match(mibig_id):
                    continue
                core_hits = int(score.get("core_gene_hits", 0) or 0)
                if core_hits < 1:
                    continue
                hits.append({
                    "contig_id": contig_id,
                    "region_id": region_id,
                    "query_label": "region_{0}".format(region_id),
                    "best_mibig_id": mibig_id,
                    "best_mibig_product": str(reference.get("description", "")),
                    "best_mibig_class": str(reference.get("cluster_type", "")),
                    "mibig_similarity": None,
                    "match_score": int(score.get("similarity", 0) or 0),
                    "score_metric": "KnownClusterBlast empirical score",
                    "matched_genes": int(score.get("hits", 0) or 0),
                    "core_gene_hits": core_hits,
                    "evidence_source": "knownclusterblast",
                })
                break

        compare = modules.get("antismash.modules.cluster_compare", {}) or {}
        mibig = (compare.get("db_results", {}) or {}).get("MIBiG", {}) or {}
        for region_id, modes in (mibig.get("by_region", {}) or {}).items():
            regional = (modes.get("RegionToRegion_RiQ", {}) or {}).get("scores_by_region", {}) or {}
            choices = []
            for reference, score in regional.items():
                match = MIBIG_ID_RE.match(str(reference))
                if match:
                    choices.append((float(score), match.group(1)))
            if not choices:
                continue
            best_score, best_id = max(choices, key=lambda item: (item[0], item[1]))
            hits.append({
                "contig_id": contig_id,
                "region_id": str(region_id),
                "query_label": "region_{0}".format(region_id),
                "best_mibig_id": best_id,
                "best_mibig_product": "",
                "best_mibig_class": "",
                "mibig_similarity": None,
                "match_score": best_score,
                "score_metric": "ClusterCompare score (0-1)",
                "matched_genes": None,
                "core_gene_hits": None,
                "evidence_source": "clustercompare_mibig",
            })
    return hits
