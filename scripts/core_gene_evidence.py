"""Structured biosynthetic core-gene evidence shared by predictor parsers."""

import re


CORE_RECORD_SEPARATOR = ";;"
CORE_FIELD_SEPARATOR = "|"

# Domain families with a direct role in scaffold formation. Broad tailoring
# domains are intentionally excluded.
CORE_BIOSYNTHETIC_PFAMS = {
    "PF00109": "beta-ketoacyl synthase",
    "PF00698": "acyltransferase",
    "PF02801": "beta-ketoacyl synthase C-terminal",
    "PF00501": "AMP-binding enzyme",
    "PF00550": "phosphopantetheine carrier protein",
    "PF00975": "thioesterase",
    "PF01397": "terpene synthase",
    "PF03936": "terpene synthase metal-binding",
    "PF05147": "lanthionine synthetase C-terminal",
    "PF14028": "lanthipeptide dehydratase",
}

# Bakta is a fallback only. These phrases identify scaffold-forming enzymes;
# generic terms such as synthetase, oxidase, transferase and cyclase are absent.
BAKTA_CORE_PATTERNS = (
    ("non-ribosomal peptide synthetase", "non-ribosomal peptide synthetase"),
    ("nonribosomal peptide synthetase", "nonribosomal peptide synthetase"),
    ("polyketide synthase", "polyketide synthase"),
    ("terpene synthase", "terpene synthase"),
    ("lanthipeptide synthetase", "lanthipeptide synthetase"),
    ("lanthionine synthetase", "lanthionine synthetase"),
    ("lasso peptide synthetase", "lasso peptide synthetase"),
    ("thiopeptide synthetase", "thiopeptide synthetase"),
    ("bacteriocin synthetase", "bacteriocin synthetase"),
    ("siderophore synthetase", "siderophore synthetase"),
    ("peptide synthetase", "peptide synthetase"),
)


def clean(value):
    text = str(value or "").strip()
    return "" if text.lower() in {"nan", "none"} else text


def first_value(value):
    if isinstance(value, list):
        return clean(value[0]) if value else ""
    return clean(value)


def parse_location(location):
    match = re.search(r"\[<?(\d+):>?(\d+)\]", clean(location))
    if not match:
        return "", ""
    # antiSMASH JSON/Biopython locations are zero-based, end-exclusive.
    return int(match.group(1)) + 1, int(match.group(2))


def make_record(identifier, start, end, source, evidence):
    fields = [clean(identifier), clean(start), clean(end), clean(source), clean(evidence)]
    return CORE_FIELD_SEPARATOR.join(field.replace(CORE_FIELD_SEPARATOR, "/") for field in fields)


def join_records(records):
    return CORE_RECORD_SEPARATOR.join(dict.fromkeys(record for record in records if clean(record)))


def parse_records(value):
    records = []
    for raw in clean(value).split(CORE_RECORD_SEPARATOR):
        if not raw:
            continue
        fields = raw.split(CORE_FIELD_SEPARATOR, 4)
        if len(fields) != 5:
            continue
        identifier, start, end, source, evidence = fields
        records.append({
            "identifier": identifier,
            "start": int(float(start)) if start else None,
            "end": int(float(end)) if end else None,
            "source": source,
            "evidence": evidence,
        })
    return records


def antismash_core_records(features):
    records = []
    for feature in features or []:
        if feature.get("type") != "CDS":
            continue
        qualifiers = feature.get("qualifiers", {}) or {}
        kinds = [clean(item).lower() for item in qualifiers.get("gene_kind", [])]
        functions = [clean(item) for item in qualifiers.get("gene_functions", [])]
        biosynthetic_functions = [item for item in functions if item.lower().startswith("biosynthetic ")]
        if "biosynthetic" not in kinds and not biosynthetic_functions:
            continue
        start, end = parse_location(feature.get("location", ""))
        locus = first_value(qualifiers.get("locus_tag", ""))
        evidence = "; ".join(biosynthetic_functions) or "gene_kind=biosynthetic"
        records.append(make_record(locus, start, end, "antismash", evidence))
    return records


def domain_core_record(identifier, start, end, source, domains):
    supported = []
    for domain in domains:
        domain = clean(domain).split(".")[0]
        if domain in CORE_BIOSYNTHETIC_PFAMS:
            supported.append("{0} ({1})".format(domain, CORE_BIOSYNTHETIC_PFAMS[domain]))
    if not supported:
        return ""
    return make_record(identifier, start, end, source, "; ".join(dict.fromkeys(supported)))


def bakta_core_match(row):
    text = " ".join(clean(row.get(field)) for field in ("gene", "product", "dbxrefs")).lower()
    for phrase, label in BAKTA_CORE_PATTERNS:
        if phrase in text:
            return label
    return ""
