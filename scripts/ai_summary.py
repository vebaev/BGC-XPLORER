"""Readable cluster summaries built on facts computed by code, with an AI hypothesis kept apart.

Code reads the run's own outputs (the locus table, antiSMASH's JSON with gene
roles, NRPS/PKS domains, modules, substrate predictions and KnownClusterBlast
gene pairs, and the MIBiG 4.0 entry of the closest known cluster) and states
each finding as a numbered fact. The model writes a short description in
three parts (what the cluster is, what it contains, what it resembles), citing
for every sentence the facts it uses. Each sentence is checked on its own:
every number, locus tag, accession, caller, substrate, domain and compound
name must come from the facts it cites, a locus tag's attributes must come
from a fact about that gene, and interpretive wording is not allowed. A
sentence that fails is removed; a part left empty is filled from the facts as
written by code. Short definitions of technical terms come from a fixed
glossary and are added by code, not by the model.

The model may also offer a hypothesis about what the cluster makes or does.
It is shown separately and labelled as an unverified AI hypothesis; it is
checked only for identifiers, numbers and names that are not in the facts, and
for being worded as a hypothesis.
"""
import json
import os
import re

import ai_evidence
from ai_facts import CALLER_NAMES, _caller_list, _clean, _core_support, _digits, _int, caller_groups

PROMPT_VERSION = "5.8-summary"
SECTIONS = (("overview", "What it is"), ("genes", "What it contains"), ("similar", "What it resembles"))
MAX_CORE_GENES = 6
MAX_ADDITIONAL_GENES = 5
MAX_CONTEXT_GENES = 80
EDGE_DISTANCE_BP = 10000

GLOSSARY = {
    "module": "the unit of an assembly-line enzyme that adds one building block to the growing chain",
    "adenylation domain": "selects and activates the amino acid a module adds",
    "epimerization domain": "converts an amino acid to its D-form",
    "condensation domain": "joins the amino acid of its module to the growing chain",
    "thioesterase": "an enzyme type that releases the finished chain from NRPS and PKS assembly lines",
    "MbtH-like protein": "a small helper protein of NRPS adenylation domains",
    "MbtH family protein": "a small helper protein of NRPS adenylation domains",
    "cytochrome P450": "a family of oxidising enzymes",
    "ketosynthase domain": "extends a polyketide chain by one unit",
    "acyltransferase domain": "selects the extender unit a polyketide module adds",
    "glycosyltransferase": "an enzyme type that attaches sugars",
    "methyltransferase": "an enzyme type that adds methyl groups",
    "terpene synthase": "an enzyme type that folds a linear precursor into a terpene skeleton",
}
TERM_PATTERNS = {term: re.compile(r"\b" + re.escape(term) + r"s?\b", re.IGNORECASE) for term in GLOSSARY}

DOMAIN_NAMES = {"AMP-binding": "A", "PCP": "PCP", "ACP": "ACP", "Epimerization": "E", "Thioesterase": "TE",
                "PKS_KS": "KS", "PKS_AT": "AT", "PKS_KR": "KR", "PKS_DH": "DH", "PKS_ER": "ER",
                "Heterocyclization": "Cy", "cMT": "cMT", "nMT": "nMT", "oMT": "oMT", "TD": "TD"}
AMINO_ACIDS = {"Ala": "alanine", "Arg": "arginine", "Asn": "asparagine", "Asp": "aspartate", "Cys": "cysteine",
               "Gln": "glutamine", "Glu": "glutamate", "Gly": "glycine", "His": "histidine", "Ile": "isoleucine",
               "Leu": "leucine", "Lys": "lysine", "Met": "methionine", "Phe": "phenylalanine", "Pro": "proline",
               "Ser": "serine", "Thr": "threonine", "Trp": "tryptophan", "Tyr": "tyrosine", "Val": "valine",
               "Orn": "ornithine", "Pip": "pipecolate", "Hpg": "4-hydroxyphenylglycine", "Dhpg": "3,5-dihydroxyphenylglycine",
               "Dab": "2,4-diaminobutyrate", "Bht": "beta-hydroxytyrosine", "Aad": "2-aminoadipate"}
EXTENDERS = {"mal": "malonyl-CoA", "mmal": "methylmalonyl-CoA", "emal": "ethylmalonyl-CoA",
             "mxmal": "methoxymalonyl-ACP", "ohmal": "hydroxymalonyl-ACP"}
CLASS_PHRASES = {"NRPS": "nonribosomal peptide (NRPS)", "NRPS-like": "NRPS-like", "T1PKS": "type I polyketide (T1PKS)",
                 "T2PKS": "type II polyketide (T2PKS)", "T3PKS": "type III polyketide (T3PKS)", "terpene": "terpene",
                 "lanthipeptide-class-i": "class I lanthipeptide", "lanthipeptide-class-ii": "class II lanthipeptide",
                 "lanthipeptide-class-iii": "class III lanthipeptide", "RiPP-like": "RiPP-like",
                 "siderophore": "siderophore", "NI-siderophore": "NRPS-independent siderophore",
                 "ectoine": "ectoine", "butyrolactone": "butyrolactone", "melanin": "melanin",
                 "Polyketide": "polyketide", "NRP": "nonribosomal peptide", "Terpene": "terpene", "RiPP": "RiPP",
                 "Saccharide": "saccharide", "Alkaloid": "alkaloid"}
ROLE_WORDS = (("biosynthetic", "core biosynthetic"), ("biosynthetic-additional", "additional biosynthetic"),
              ("regulatory", "regulatory"), ("transport", "transport"), ("resistance", "resistance"))
INTERPRETIVE = re.compile(
    r"\b(likely|suggests?|suggesting|may|might|probably|possibly|could|appears?|presumably|"
    r"putative|produces?|producing|synthesi[sz]es?|responsible|novel|antibiotics?|role|function(?:s|al)?)\b",
    re.IGNORECASE)
HEDGE = re.compile(r"\b(may|might|could|possibly|perhaps|suggests?|hypothes[ie]s|consistent with|would)\b", re.IGNORECASE)
COMPOUND_LIKE = re.compile(r"\b[a-z][a-z-]*(?:mycin|micin|cidin|statin|lactam|peptin|bactin|chelin)s?\b", re.IGNORECASE)
AMINO_CODE = re.compile(r"\b(" + "|".join(AMINO_ACIDS) + r")\b")
CALLER_PATTERN = re.compile(r"\b(antismash|gecco|deepbgc)\b", re.IGNORECASE)
CLOSENESS = re.compile(r"\b(close(?:ly)?|moderate(?:ly)?|distant(?:ly)?|weak(?:ly)?|strong(?:ly)?|partial(?:ly)?|"
                       r"high|low|good|poor|remote(?:ly)?|near|loose(?:ly)?)\b", re.IGNORECASE)
STOP_WORDS = set("a an the of and or as at in on to for with by from is are be its it this that these those each all "
                 "both total about only also which who whose than then there their them while whereas including "
                 "comprising plus gene genes protein proteins cluster locus".split())
WORD = re.compile(r"[a-z][a-z-]*")
EXCLUSION = re.compile(r"\b(remaining|the rest|other(?:s)?|further|additional)\b", re.IGNORECASE)
PUBMED = re.compile(r"\b(?:PubMed|PMID)[:\s]*(\d{6,9})\b", re.IGNORECASE)


def _num(value):
    n = ai_evidence._as_number(_clean(value))
    return None if n is None else (int(n) if n == int(n) else n)


def _fmt(n):
    return "{0:,}".format(n) if isinstance(n, int) else str(n)


def _plural(n, word, plural=None):
    return "{0} {1}".format(n, word if n == 1 else (plural or word + "s"))


def _join(items):
    items = [i for i in items if i]
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + " and " + items[-1] if items else ""


def _location(text):
    numbers = [int(n) for n in re.findall(r"\d+", text or "")]
    return (numbers[0], numbers[1]) if len(numbers) >= 2 else (None, None)


def fact(fid, section, sentence, values=(), genes=(), terms=(), names=(), detail="", callers=(), total=False):
    """A numbered fact: its sentence as code writes it, the checkable values it holds, the genes it is about,
    the glossary terms it uses, the names (substrates, compounds) it contains, and detail shown to the reader
    on hover but not given to the model."""
    return {"id": fid, "section": section, "sentence": sentence, "values": [str(v) for v in values if v not in ("", None)],
            "genes": list(genes), "terms": list(terms), "names": [n for n in names if n], "detail": detail,
            "callers": list(callers), "total": total}


# ---------------------------------------------------------------- inputs

def antismash_json_path(results_dir, sample):
    return os.path.join(str(results_dir), sample, "antismash", sample + ".json")


def load_antismash(path):
    if not path or not os.path.exists(path):
        return None
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def mibig_entry(mibig_dir, accession):
    if not mibig_dir or not accession:
        return None
    path = os.path.join(mibig_dir, accession.split(".")[0] + ".json")
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def default_mibig_dir():
    for path in (os.environ.get("MIBIG_JSON_DIR", ""), "/db/mibig/mibig_json_4.0"):
        if path and os.path.isdir(path):
            return path
    return ""


def _record_for(antismash, contig):
    for record in (antismash or {}).get("records", []):
        if record.get("id") == contig or record.get("name") == contig:
            return record
    records = (antismash or {}).get("records", [])
    return records[0] if len(records) == 1 else None


def antismash_locus(antismash, contig, start, end):
    """Gene roles, domains, modules and predictions for the genes of one locus, and its KnownClusterBlast result."""
    record = _record_for(antismash, contig)
    if not record:
        return None
    genes, domains_by_tag = {}, {}
    for feature in record.get("features", []):
        q = feature.get("qualifiers", {})
        if feature.get("type") == "CDS":
            a, b = _location(feature.get("location"))
            if a is None or not start <= (a + b) / 2 <= end:
                continue
            tag = q.get("locus_tag", [""])[0]
            core_types = sorted({g.split(") ", 1)[1].split(":")[0] for g in q.get("gene_functions", [])
                                 if g.startswith("biosynthetic (") and ") " in g})
            genes[tag] = {"kind": q.get("gene_kind", [""])[0], "product": q.get("product", [""])[0], "core_types": core_types,
                          "domains": [d.split("Domain: ", 1)[1].split(" (", 1)[0] for d in q.get("NRPS_PKS", [])
                                      if d.startswith("Domain: ")], "modules": []}
    regions = []
    for feature in record.get("features", []):
        q = feature.get("qualifiers", {})
        if feature.get("type") == "aSModule":
            tags = q.get("locus_tags", [])
            if tags and tags[0] in genes:
                genes[tags[0]]["modules"].append({"type": q.get("type", [""])[0], "domains": q.get("domains", [])})
        elif feature.get("type") == "region":
            a, b = _location(feature.get("location"))
            if a is not None and a < end and b > start:
                regions.append(int(q.get("region_number", ["0"])[0]))
    consensus = record.get("modules", {}).get("antismash.modules.nrps_pks", {}).get("consensus", {})
    for gene in genes.values():
        for module in gene["modules"]:
            module["substrates"] = [("A" if "AMP-binding" in d else "AT", consensus.get(d))
                                    for d in module["domains"] if "AMP-binding" in d or "PKS_AT" in d]
            module["epimerization"] = any("Epimerization" in d for d in module["domains"])
    known = record.get("modules", {}).get("antismash.modules.clusterblast", {}).get("knowncluster", {})
    kcb = [r for r in known.get("results", []) if r.get("region_number") in regions]
    return {"genes": genes, "kcb": kcb}


def kcb_match(kcb_results, accession):
    """The KnownClusterBlast ranking entry for one MIBiG accession, with each query gene's best pair."""
    base = accession.split(".")[0]
    for result in kcb_results:
        for cluster, score in result.get("ranking", []):
            if cluster.get("accession", "").split(".")[0] == base:
                best = {}
                for query, _, subject in score.get("pairings", []):
                    tag = query.split("|")[4] if query.count("|") >= 4 else query
                    if tag not in best or subject.get("perc_ident", 0) > best[tag]["perc_ident"]:
                        best[tag] = subject
                return {"similarity": score.get("similarity"), "hits": score.get("hits"),
                        "core_gene_hits": score.get("core_gene_hits"), "pairs": best}
    return None


# ---------------------------------------------------------------- facts

def _substrate(code):
    if not code or code in ("X", "pk"):
        return ""
    return AMINO_ACIDS.get(code, EXTENDERS.get(code, code))


def _short_domain(name):
    return "C" if name.startswith("Condensation") else DOMAIN_NAMES.get(name, "")


def _arch(gene):
    return "–".join(s for s in (_short_domain(d) for d in gene["domains"]) if s)


def _counted(names):
    """'phenylalanine, no substrate and valine'; with repeats 'malonyl-CoA for 5 and no substrate for 1'."""
    names = [n or "no substrate" for n in names]
    order = sorted(set(names), key=names.index)
    if len(order) == len(names):
        return _join(names)
    if len(order) == 1:
        return "{0} for all {1}".format(order[0], len(names))
    return _join(["{0} for {1}".format(n, names.count(n)) for n in order])


def _substrate_phrase(substrates):
    """What antiSMASH predicts for a set of (domain kind, code) pairs, by kind of domain."""
    parts = []
    for kind, label in (("A", "adenylation"), ("AT", "acyltransferase")):
        names = [_substrate(code) for k, code in substrates if k == kind]
        if names:
            parts.append("{0} for its {1} domain{2}".format(_counted(names), label, "s" if len(names) > 1 else ""))
    return _join(parts)


def _match_word(similarity, identities):
    """How close a KnownClusterBlast match is, from the share of the MIBiG cluster's genes matched and the
    median identity of the matched genes."""
    idents = sorted(identities)
    median = idents[len(idents) // 2] if idents else 0
    if similarity >= 75 and median >= 80:
        return "close"
    if similarity >= 25 and median >= 50:
        return "moderate"
    return "distant"


def _gene_terms(product):
    return [t for t, p in TERM_PATTERNS.items() if p.search(product or "") and t not in ("module",)]


def build_facts(payload, antismash=None, mibig_dir=""):
    region = payload.get("region", {}) or {}
    genes = payload.get("genes", []) or []
    label = _clean(payload.get("consensus_label")) or _clean(payload.get("consensus_id")).replace("_consensus_", "_")
    contig = _clean(region.get("contig"))
    start, end, length = _int(region.get("start")), _int(region.get("end")), _int(region.get("length_bp"))
    callers = [c.strip().lower() for c in _clean(region.get("support_tools")).split(",") if c.strip()]
    facts, n = [], {"overview": 0, "genes": 0, "similar": 0}

    def add(section, sentence, **kw):
        n[section] += 1
        facts.append(fact("{0}{1}".format(section[0].upper(), n[section]), section, sentence, **kw))

    # ---- what it is
    types = [t.strip() for t in re.split(r"[,;/]", _clean(region.get("bgc_types")))
             if t.strip() and t.strip().lower() not in ("unknown", "other", "none")]
    phrases = []
    for t in types:
        phrase = CLASS_PHRASES.get(t, t)
        if not any(phrase in other or other in phrase for other in phrases):
            phrases.append(phrase)
    who = _caller_list(callers) + (" (all three detection tools)" if len(callers) == 3 else " only" if len(callers) == 1 else "")
    text = ("{0} is a {1} cluster, found by {2}.".format(label, _join(phrases), who) if phrases else
            "{0} is a biosynthetic gene cluster, found by {1}; no compound class is assigned to it.".format(label, who))
    add("overview", text,
        values=[CALLER_NAMES.get(c, c) for c in callers] + ([3] if len(callers) == 3 else []), names=types)
    in_span = [g for g in genes if start is not None and _int(g.get("gene_start")) is not None and
               start <= (_int(g.get("gene_start")) + _int(g.get("gene_end"))) / 2 <= end]
    kb = int(round(length / 1000.0)) if length else None
    add("overview", "It spans {0} kb on {1} ({2}–{3}) and contains {4}.".format(
        kb, contig, _fmt(start), _fmt(end), _plural(len(in_span), "gene")), values=[kb, start, end, len(in_span), contig])
    edge = _int(region.get("nearest_contig_edge_bp"))
    if edge is not None and edge < EDGE_DISTANCE_BP:
        add("overview", "It lies {0} bp from the end of its contig, so it may be incomplete.".format(_fmt(edge)),
            values=[edge])

    # ---- what it contains
    local = antismash_locus(antismash, contig, start, end) if antismash and "antismash" in callers else None
    as_genes = local["genes"] if local else {}
    products = {_clean(g.get("locus_tag")): _clean(g.get("bakta_product")) for g in genes}
    roles = {}
    for tag, gene in as_genes.items():
        if gene["kind"]:
            roles.setdefault(gene["kind"], []).append(tag)
    if roles:
        core, extra_n = len(roles.get("biosynthetic", [])), len(roles.get("biosynthetic-additional", []))
        parts = []
        if core or extra_n:
            parts.append("{0} as biosynthetic ({1} core and {2} additional)".format(_plural(core + extra_n, "gene"), core, extra_n))
        parts += ["{0} as {1}".format(len(roles[k]), word) for k, word in ROLE_WORDS[2:] if roles.get(k)]
        add("genes", "antiSMASH marks {0}.".format(_join(parts)),
            values=[core + extra_n, core, extra_n] + [len(roles[k]) for k, _ in ROLE_WORDS[2:] if roles.get(k)])
    else:
        counts = {}
        for g in in_span:
            counts[_clean(g.get("gene_category"))] = counts.get(_clean(g.get("gene_category")), 0) + 1
        parts = ["{0} {1}".format(counts[k], w) for k, w in (("regulator", "regulatory"), ("transporter", "transport"),
                 ("tailoring_enzyme", "tailoring-enzyme"), ("resistance", "resistance")) if counts.get(k)]
        if parts:
            add("genes", "By the annotation categories of BGC-XPLORER, it holds {0} genes.".format(_join(parts)),
                values=[counts[k] for k in ("regulator", "transporter", "tailoring_enzyme", "resistance") if counts.get(k)])

    # Only core genes inside the agreed span count as core genes of the locus; callers sometimes name core genes
    # in the wider region they predicted, and those are reported apart.
    span_tags = {_clean(g.get("locus_tag")) for g in in_span}
    support_all = _core_support(region.get("core_gene_support"))
    support = {t: c for t, c in support_all.items() if t in span_tags}
    outside = {t: c for t, c in support_all.items() if t not in span_tags}
    core_tags = list(support)[:MAX_CORE_GENES]
    if support:
        by = sorted({c for tags in support.values() for c in tags})
        if len(by) == 1:
            text = "{0} names {1} in the locus.".format(_caller_list(by), _plural(len(support), "core gene"))
        elif len(support) == 1:
            text = "1 core gene is named in the locus."
        else:
            text = "{0} core genes are named in the locus, each by at least one of {1}{2}.".format(
                len(support), "the three detection tools, " if len(by) == 3 else "", _caller_list(by))
        add("genes", text, values=[len(support)] + [CALLER_NAMES.get(c, c) for c in by])
    else:
        add("genes", "No detection tool names a core biosynthetic gene inside the agreed span of this locus.")
    if outside:
        by = sorted({c for tags in outside.values() for c in tags})
        add("genes", "{0} named by {1} lie{2} outside the agreed span of the locus.".format(
            _plural(len(outside), "further core gene"), _caller_list(by), "s" if len(outside) == 1 else ""),
            values=[len(outside)])
    for tag in core_tags:
        product = products.get(tag) or as_genes.get(tag, {}).get("product", "")
        gene = as_genes.get(tag)
        names = [CALLER_NAMES.get(c, c) for c in support[tag]]
        text = "{0}{1} is named as a core gene by {2}".format(
            tag, ", annotated by Bakta as {0},".format(product) if product else "", _caller_list(support[tag]))
        values, terms, gnames, detail = [tag] + names, _gene_terms(product), [], ""
        if gene and gene.get("core_types") and gene["kind"] == "biosynthetic":
            text += "; antiSMASH classes it as a core {0} gene".format(_join(gene["core_types"]))
            gnames += gene["core_types"]
        if gene and gene["modules"]:
            m = len(gene["modules"])
            text += "; antiSMASH finds {0} in it".format(_plural(m, "module"))
            detail = "antiSMASH domain order: " + _arch(gene)
            pairs = [p for mod in gene["modules"] for p in mod["substrates"]]
            known = [_substrate(c) for _, c in pairs if _substrate(c)]
            if pairs:
                text += ", and predicts {0}".format(_substrate_phrase(pairs))
            gene_epi = sum(mod["epimerization"] for mod in gene["modules"])
            if gene_epi:
                text += "; {0} of its modules carr{1} an epimerization domain".format(gene_epi, "ies" if gene_epi == 1 else "y")
            values += [m]
            terms += ["module"] + (["epimerization domain"] if any(mod["epimerization"] for mod in gene["modules"]) else [])
            gnames += known
        add("genes", text + ".", values=values, genes=[tag], terms=terms, names=gnames, detail=detail,
            callers=support[tag])

    nrps = [g for t, g in as_genes.items() if g["modules"]]
    if len(nrps) > 1:
        modules = [mod for g in nrps for mod in g["modules"]]
        subs = [_substrate(c) for mod in modules for _, c in mod["substrates"]]
        known = [x for x in subs if x]
        epi = sum(mod["epimerization"] for mod in modules)
        tags = [t for t, g in as_genes.items() if g["modules"]]
        kinds = sorted({k.upper() for g in nrps for mod in g["modules"] for k in [mod["type"]] if k})
        text = "Together these {0} {5}genes with modules ({1}) hold {2} in all; antiSMASH predicts a substrate for {3} ({4})".format(
            len(tags), _join(tags), _plural(len(modules), "module"), _plural(len(known), "module"),
            _counted(known) if known else "none", (_join(kinds) + " ") if kinds else "")
        if len(subs) > len(known):
            text += " and none for {0}".format(_plural(len(subs) - len(known), "module"))
        if len(modules) > len(subs):
            text += "; {0} no adenylation or acyltransferase domain".format(
                _plural(len(modules) - len(subs), "module has", "modules have"))
        if epi:
            text += "; {0} an epimerization domain".format(_plural(epi, "module carries", "modules carry"))
        add("genes", text + ".", values=tags + [len(tags), len(modules), len(known), len(subs) - len(known), epi], genes=tags,
            total=True, terms=["module"] + (["epimerization domain"] if epi else []), names=known + kinds)

    extra = [t for t in roles.get("biosynthetic-additional", []) if t not in core_tags]
    extra = [t for t in extra if products.get(t) and "hypothetical" not in products[t].lower()]
    # Genes whose product the glossary explains come first; then the rest in genome order.
    extra = sorted(extra, key=lambda t: (not _gene_terms(products[t]), extra.index(t)))[:MAX_ADDITIONAL_GENES]
    for tag in extra:
        add("genes", "{0} encodes {1} and is marked by antiSMASH as an additional biosynthetic gene.".format(
            tag, products[tag]), values=[tag], genes=[tag], terms=_gene_terms(products[tag]))

    listed = ["{0} {1}".format(_clean(g.get("locus_tag")), _clean(g.get("bakta_product")) or "no product")
              for g in in_span][:MAX_CONTEXT_GENES]
    if listed:
        facts.append(fact("H1", "context", "Gene products in the locus, in genome order: " + "; ".join(listed) + "."))

    known_hits, duf_hits = _int(region.get("arts_known_hits")) or 0, _int(region.get("arts_duf_hits")) or 0
    if known_hits:
        add("genes", "ARTS finds {0} to known resistance models.".format(_plural(known_hits, "hit")), values=[known_hits])

    # ---- what it resembles
    mibig = payload.get("mibig_dereplication", {}) or {}
    accession = _clean(mibig.get("best_mibig_id") or region.get("best_mibig_id"))
    if not accession:
        if "antismash" in callers:
            add("similar", "antiSMASH found no similar cluster in MIBiG.")
        else:
            add("similar", "This locus has not been compared with MIBiG: the comparison is made by antiSMASH, "
                           "which did not detect it.")
        return facts
    entry = mibig_entry(mibig_dir, accession) or {}
    acc = accession.split(".")[0]
    compounds = [c.get("name") for c in entry.get("compounds", []) if c.get("name")]
    compound = _join(compounds[:3]) or _clean(mibig.get("best_mibig_product")).split("/")[0]
    organism = (entry.get("taxonomy") or {}).get("name", "")
    add("similar", "The closest MIBiG entry is {0} ({1}){2}.".format(
        accession.split(".")[0], compound, ", from {0}".format(organism) if organism else ""),
        values=[accession.split(".")[0]], names=compounds[:3] + [compound, organism])
    match = kcb_match(local["kcb"], accession) if local else None
    if match:
        match["pairs"] = {t: p for t, p in match["pairs"].items() if t in span_tags}
    if match and match["pairs"]:
        idents = [p.get("perc_ident") for p in match["pairs"].values() if p.get("perc_ident") is not None]
        sim = match["similarity"]
        word = _match_word(sim, idents)
        span = "{0} %".format(min(idents)) if min(idents) == max(idents) else "{0}–{1} %".format(min(idents), max(idents))
        add("similar", "{0} of the locus {1} a counterpart in {2}, at {3} identity; antiSMASH finds counterparts for "
                       "{4} % of that cluster's genes. Overall this is a {5} match.".format(
                           _plural(len(match["pairs"]), "gene"), "has" if len(match["pairs"]) == 1 else "have",
                           accession.split(".")[0], span, sim, word),
            values=[len(match["pairs"]), min(idents), max(idents), sim, accession.split(".")[0]], names=[word])
        for tag in core_tags:
            pair = match["pairs"].get(tag)
            if pair:
                annotation = pair.get("annotation", "").replace("_", " ")
                add("similar", "{0} matches {1}{2} of {4} at {3} % identity.".format(
                    tag, pair.get("name"), " ({0})".format(annotation) if annotation else "", pair.get("perc_ident"),
                    accession.split(".")[0]),
                    values=[tag, pair.get("name"), pair.get("perc_ident")], genes=[tag], names=[annotation])
    else:
        metric, score = _clean(mibig.get("score_metric")), _num(mibig.get("match_score"))
        if metric and score is not None and "ClusterCompare" in metric:
            score = round(float(score), 2)
            word = "close" if score >= 0.8 else "moderate" if score >= 0.5 else "distant"
            add("similar", "antiSMASH's ClusterCompare gives this entry a score of {0} on a scale of 0 to 1. "
                           "Overall this is a {1} match.".format(score, word), values=[score], names=[word])
        elif metric and score is not None:
            add("similar", "antiSMASH reports a {0} of {1} for this entry.".format(metric, score), values=[score])
    classes = [c.get("class") for c in (entry.get("biosynthesis") or {}).get("classes", []) if c.get("class")]
    if classes:
        add("similar", "MIBiG lists the compound class of {0} as {1}.".format(acc, _join(classes)), values=[acc], names=classes)
    activities = sorted({_activity(b) for c in entry.get("compounds", [])
                         for b in c.get("bioactivities", []) or [] if isinstance(b, dict) and b.get("observed")} - {""})
    if activities:
        add("similar", "MIBiG records {0} activity for {1} ({2}).".format(_join(activities), compound, acc), values=[acc],
            names=activities + [compound])
    methods = sorted({e.get("method", "") for locus in entry.get("loci", []) for e in locus.get("evidence", [])} - {""})
    refs = [r.split(":", 1)[1] for r in entry.get("legacy_references", []) if r.startswith("pubmed:")]
    if methods:
        add("similar", "MIBiG gives {0} as evidence for {1}{2}.".format(
            _join([m.lower() for m in methods]), acc, " (PubMed {0})".format(refs[0]) if refs else ""), values=refs[:1] + [acc])
    quality = entry.get("quality")
    if quality and quality != "high":
        add("similar", "MIBiG rates the annotation quality of {0} as {1}.".format(acc, quality), values=[acc])
    return facts


def _activity(bioactivity):
    """MIBiG 4.0 stores an activity name either as text or as {"activity": text}."""
    name = bioactivity.get("name", "")
    return (name.get("activity", "") if isinstance(name, dict) else str(name or "")).strip()


def gene_roles(payload, antismash=None):
    """Rows for the "Genes by role" table in the report."""
    region = payload.get("region", {}) or {}
    start, end = _int(region.get("start")), _int(region.get("end"))
    local = antismash_locus(antismash, _clean(region.get("contig")), start, end) if antismash else None
    support = _core_support(region.get("core_gene_support"))
    rows = []
    for g in payload.get("genes", []) or []:
        tag = _clean(g.get("locus_tag"))
        a, b = _int(g.get("gene_start")), _int(g.get("gene_end"))
        if start is None or a is None or not start <= (a + b) / 2 <= end:
            continue
        kind = (local or {}).get("genes", {}).get(tag, {}).get("kind", "")
        rows.append({"locus_tag": tag, "product": _clean(g.get("bakta_product")),
                     "antismash_role": kind, "category": _clean(g.get("gene_category")),
                     "core_by": _caller_list(support[tag]) if tag in support else ""})
    return rows


# ---------------------------------------------------------------- the model

SYSTEM_PROMPT = """You describe one biosynthetic gene cluster for a biologist, from numbered facts computed by software.
Write three short parts:
- "overview": 1-2 sentences, what the cluster is (class, which tools found it, size).
- "genes": 2-5 sentences, what it contains: how many biosynthetic genes there are and of which kinds, the core genes and
  their make-up (modules, domains, predicted substrates), and the most informative additional genes. Do not list every
  gene; choose. Prefer words a biologist reads easily (e.g. "phenylalanine" over "Phe").
- "similar": 1-3 sentences, the closest known cluster in MIBiG and how close it is, or that there is none.
Rules for these three parts:
- Write for a reader, in natural, varied prose, not a list. Sentences of 15-30 words; one idea per sentence.
- Start "genes" with the overall picture (how many biosynthetic genes, of which kinds), then the core genes (describe at
  most three one by one and use the "Together" fact for the rest), then at most four additional genes, grouped by kind (e.g. "two oxidising enzymes, a cytochrome P450 and an L-asparagine oxygenase").
- Say modules and substrates in words ("three modules, predicted to load phenylalanine, an unidentified residue and
  valine"); give a domain order only if it adds something. Write "an unidentified residue" for "no substrate".
- In "similar", say how close the match is in words first, then the numbers.
- After every sentence give the ids of the facts it uses. Use only values from those facts: keep every number, locus
  tag, accession, organism, compound, substrate and caller name exactly as in the facts (numbers up to twenty may be
  written as words). Do not add numbers, genes, names or comparisons of your own.
- Do not interpret: do not say what the cluster makes, what a gene does in this cluster, or whether it is new.
- Use technical terms such as module, adenylation domain, epimerization domain or thioesterase without defining them;
  the software adds the definitions.
Then a separate "hypothesis": 1-3 sentences on what the cluster might make or do, worded as a hypothesis (may, could,
possibly) and based on the facts and on the gene products listed as context; consider the whole set of genes, not
only the core genes, and say so if the genes point to primary metabolism rather than a specialised metabolite. Cite
the fact ids it rests on. It will be shown labelled as an unverified AI
hypothesis. Do not invent identifiers or numbers in it.
Return ONLY JSON matching the schema."""

_SENTENCES = {"type": "array", "minItems": 1, "items": {
    "type": "object", "additionalProperties": False, "required": ["text", "facts"],
    "properties": {"text": {"type": "string", "minLength": 10},
                   "facts": {"type": "array", "minItems": 1, "items": {"type": "string"}}}}}
JSON_SCHEMA = {
    "type": "object", "additionalProperties": False, "required": ["overview", "genes", "similar", "hypothesis"],
    "properties": {"overview": _SENTENCES, "genes": _SENTENCES, "similar": _SENTENCES,
                   "hypothesis": {"type": "object", "additionalProperties": False, "required": ["text", "facts"],
                                  "properties": {"text": {"type": "string"},
                                                 "facts": {"type": "array", "items": {"type": "string"}}}}}}


def build_user_prompt(facts):
    lines = []
    for key, title in SECTIONS:
        lines.append("{0} facts:".format(key))
        lines.extend("  [{0}] {1}".format(f["id"], f["sentence"]) for f in facts if f["section"] == key)
    context = [f for f in facts if f["section"] == "context"]
    if context:
        lines.append("context for the hypothesis only (do not use it in the three parts):")
        lines.extend("  [{0}] {1}".format(f["id"], f["sentence"]) for f in context)
    return "\n".join(lines)


def normalize(raw):
    raw = raw if isinstance(raw, dict) else {}
    out = {}
    for key, _ in SECTIONS:
        items = raw.get(key) if isinstance(raw.get(key), list) else []
        out[key] = [{"text": _clean(i.get("text")), "facts": [str(x) for x in i.get("facts", []) if x]}
                    for i in items if isinstance(i, dict) and _clean(i.get("text"))]
    hyp = raw.get("hypothesis") if isinstance(raw.get("hypothesis"), dict) else {}
    out["hypothesis"] = {"text": _clean(hyp.get("text")), "facts": [str(x) for x in hyp.get("facts", []) if x]}
    return out


# ---------------------------------------------------------------- checks

def _numbers(text):
    return set(ai_evidence._numbers(_digits(text)))


def _fact_numbers(facts):
    return set().union(*[_numbers(f["sentence"]) for f in facts]) if facts else set()


_MIBIG_NAMES = {}


def mibig_compound_names(mibig_dir=None):
    """Every one-word compound name in MIBiG (lower case), so a compound the facts do not hold is caught."""
    mibig_dir = default_mibig_dir() if mibig_dir is None else mibig_dir
    if mibig_dir not in _MIBIG_NAMES:
        names = set()
        if mibig_dir and os.path.isdir(mibig_dir):
            for name in os.listdir(mibig_dir):
                if name.endswith(".json"):
                    with open(os.path.join(mibig_dir, name), encoding="utf-8") as handle:
                        for compound in json.load(handle).get("compounds", []):
                            for word in [compound.get("name", "")] + compound.get("synonyms", []):
                                word = word.strip().lower()
                                if len(word) >= 5 and re.fullmatch(r"[a-z][a-z0-9-]*", word):
                                    names.add(word)
        _MIBIG_NAMES[mibig_dir] = names
    return _MIBIG_NAMES[mibig_dir]


def _vocabulary(facts, known_compounds=None):
    """Names that must not appear in a sentence unless a cited fact holds them: substrates, compounds,
    organisms and classes (the facts' names), every residue and extender unit the software knows, and every
    one-word compound name in MIBiG."""
    words = set(mibig_compound_names() if known_compounds is None else known_compounds)
    for f in facts:
        for name in f["names"]:
            words.update(w.lower() for w in re.findall(r"[A-Za-z][\w-]{3,}", name))
    words.update(a.lower() for a in AMINO_ACIDS.values())
    words.update(e.lower() for e in EXTENDERS.values())
    return words


def _after(text, end, count=4):
    """The first content words after a position: what a number counts."""
    window = re.split(r"[,;.:()\[\]]|\d", text[end:end + 60])[0]
    words = [w.rstrip("s") for w in WORD.findall(window.lower())]
    return [w for w in words if w not in STOP_WORDS][:count]


TAG_RANGE = re.compile(r"\b([A-Z][A-Z0-9]{2,11}_)(\d{4,6})\s*[–-]\s*(\d{4,6})\b")


def _expand_tag_ranges(text):
    """'CJLEIP_03606–03609' is two locus tags, not a tag and a number."""
    return TAG_RANGE.sub(lambda m: "{0}{1} to {0}{2}".format(m.group(1), m.group(2), m.group(3)), text)


def _number_uses(text):
    """(number, words that follow it) for every number in a text; number words count as numbers."""
    text = _digits(text)
    stripped = ai_evidence.NAME_WITH_DIGITS.sub(lambda m: " " * len(m.group(0)),
                                                ai_evidence.IDENTIFIER.sub(lambda m: " " * len(m.group(0)), text))
    uses = []
    for match in ai_evidence.NUMBER.finditer(stripped):
        value = ai_evidence._as_number(match.group(0))
        if value is not None:
            uses.append((value, _after(stripped, match.end())))
    return uses


NAMING = re.compile(r"\b(?:named|identified|called|detected|found|marked|flagged|recogni[sz]ed|predicted)\b"
                    r"(?:\s+as\s+(?:an?\s+)?(?:core|biosynthetic)[\w\s-]*?)?\s+(?:by|in)\s+(?:both\s+|all\s+three\s+of\s+)?"
                    r"((?:antiSMASH|DeepBGC|GECCO)(?:\s*(?:,|and|,\s*and)\s*(?:antiSMASH|DeepBGC|GECCO))*)", re.IGNORECASE)


def _naming_groups(text):
    """Locus tag -> sets of callers that a naming phrase ('named by antiSMASH and GECCO') gives it: the phrase
    applies to the locus tags written before it in the same sentence, since the previous naming phrase."""
    assigned = {}
    for sentence in re.split(r"(?<=[.;])\s+", text):
        last = 0
        for match in NAMING.finditer(sentence):
            callers = {c.lower() for c in CALLER_PATTERN.findall(match.group(1))}
            for tag in re.findall(ai_evidence.LOCUS_TAG, sentence[last:match.start()]):
                assigned.setdefault(tag, []).append(callers)
            last = match.end()
    return assigned


CATEGORY_STEMS = ("regulat", "transport", "tailor", "resist", "biosynth", "core", "addition", "modul", "substrat",
                  "domain", "hit", "epimeri", "cazyme", "mobile")


def _category(words):
    if any(w.startswith("epimeri") for w in words):
        return "epimeri"
    for word in words[:1]:
        for stem in CATEGORY_STEMS:
            if word.startswith(stem):
                return stem
    return None


def _check_number_binding(text, cited):
    """A number that counts a category of genes or domains must count the same category in a cited fact:
    '3 resistance genes' needs a fact with '3 resistance', not just any 3. In a sentence about particular genes the
    fact must be about those genes (or about no gene in particular), so a total is not given to one gene."""
    tags = set(re.findall(ai_evidence.LOCUS_TAG, text))
    usable = [f for f in cited if not tags or not f["genes"] or set(f["genes"]) <= tags]
    fact_uses = [(v, _category(w)) for f in usable for v, w in _number_uses(f["sentence"])]
    reasons = []
    for value, words in _number_uses(text):
        category = _category(words)
        if category is None or value != int(value):
            continue
        if (value, category) not in fact_uses:
            reasons.append("the number {0:g} counts '{1}' here, but not in the cited facts".format(value, words[0]))
    return reasons


def _listed_counts(text):
    """Numbers that count the locus tags listed right after them in the same phrase ("the two genes X and Y")."""
    digits = _digits(text)
    stripped = ai_evidence.NAME_WITH_DIGITS.sub(lambda m: " " * len(m.group(0)),
                                                ai_evidence.IDENTIFIER.sub(lambda m: " " * len(m.group(0)), digits))
    allowed = set()
    for match in ai_evidence.NUMBER.finditer(stripped):
        value = ai_evidence._as_number(match.group(0))
        phrase = re.split(r";|, and |, while |, whereas |\.\s", digits[match.end():] + " ")[0]
        by_prefix = {}
        for tag in re.findall(ai_evidence.LOCUS_TAG, phrase):
            by_prefix.setdefault(tag.rsplit("_", 1)[0], set()).add(tag)
        if value is not None and value in {float(len(tags)) for tags in by_prefix.values()}:
            allowed.add(value)
    return allowed


def check_sentence(item, facts, vocabulary, section=""):
    by_id = {f["id"]: f for f in facts}
    text, reasons = _expand_tag_ranges(item["text"]), []
    cited = [by_id[i] for i in item["facts"] if i in by_id and by_id[i]["section"] != "context"]
    if not cited:
        return ["cites no known fact"]
    pool = " ".join(f["sentence"] for f in cited)
    pool_lower = pool.lower()
    for number in _numbers(text) - _fact_numbers(cited) - _listed_counts(text):
        reasons.append("the number {0:g} is not in the cited facts".format(number))
    for identifier in set(ai_evidence.IDENTIFIER.findall(text)):
        if identifier not in pool:
            reasons.append("{0} is not in the cited facts".format(identifier))
    for caller in {c.lower() for c in CALLER_PATTERN.findall(text)}:
        if caller not in pool_lower:
            reasons.append("{0} is not in the cited facts".format(caller))
    for word in {w.lower() for w in re.findall(r"[A-Za-z][\w-]{3,}", text)} & vocabulary:
        if word not in pool_lower:
            reasons.append("'{0}' is not in the cited facts".format(word))
    for word in {w.lower() for w in COMPOUND_LIKE.findall(text)}:
        if word not in pool_lower and not any(TERM_PATTERNS[t].search(word) for t in GLOSSARY):
            reasons.append("'{0}' is not in the cited facts".format(word))
    match = INTERPRETIVE.search(text)
    if match and match.group(0).lower() not in pool_lower:
        reasons.append("interpretive wording: '{0}'".format(match.group(0)))
    if section == "similar":
        for word in {w.lower() for w in CLOSENESS.findall(text)}:
            if not re.search(r"\b" + re.escape(word[:-2] if word.endswith("ly") else word), pool_lower):
                reasons.append("'{0}' judges the match, but the cited facts do not".format(word))
    reasons.extend(_check_number_binding(text, cited))
    if any(f.get("total") for f in cited) and EXCLUSION.search(text):
        reasons.append("'{0}' gives a total for all these genes to only some of them".format(EXCLUSION.search(text).group(0)))
    # Where a sentence says which callers name a core gene, they must be exactly the callers that name it.
    known_callers = {f["genes"][0]: set(f["callers"]) for f in facts if f.get("callers") and f["genes"]}
    for tag, groups in _naming_groups(text).items():
        for named in groups:
            if tag in known_callers and named != known_callers[tag]:
                reasons.append("{0} is said to be named by {1}, but is named by {2}".format(
                    tag, ", ".join(sorted(named)), ", ".join(sorted(known_callers[tag]))))
    for code in set(AMINO_CODE.findall(text)):
        if code not in pool:
            reasons.append("'{0}' is not in the cited facts".format(code))
    # A gene's attributes must come from a fact about that gene (or from a fact about no gene).
    tags = set(re.findall(ai_evidence.LOCUS_TAG, text))
    if tags:
        allowed = " ".join(f["sentence"] for f in cited if not f["genes"] or set(f["genes"]) & tags).lower()
        for name in {n.lower() for f in cited for n in f["names"]}:
            if name in text.lower() and name not in allowed:
                reasons.append("'{0}' is given to a gene whose fact does not hold it".format(name))
    return reasons


def check_hypothesis(hyp, facts):
    text, reasons = _expand_tag_ranges(hyp.get("text", "")), []
    if not text:
        return ["no hypothesis"]
    pool = " ".join(f["sentence"] for f in facts)
    if not HEDGE.search(text):
        reasons.append("not worded as a hypothesis")
    for number in _numbers(text) - _fact_numbers(facts):
        reasons.append("the number {0:g} is not in the facts".format(number))
    for identifier in set(ai_evidence.IDENTIFIER.findall(text)):
        if identifier not in pool:
            reasons.append("{0} is not in the facts".format(identifier))
    if len(re.findall(r"[.!?](?:\s|$)", text)) > 3:
        reasons.append("longer than three sentences")
    return reasons


def template_section(facts, section):
    return [{"text": f["sentence"], "facts": [f["id"]], "source": "template"} for f in facts if f["section"] == section]


def _stands_alone(text, match):
    """A term is defined only where it stands as a noun of its own: not inside brackets, not inside a longer
    product name, and not already followed by a bracket."""
    before, after = text[:match.start()], text[match.end():]
    if before.count("(") > before.count(")"):
        return False
    return bool(re.match(r"\s*(?:$|[,.;:]|and\b|or\b|in\b|with\b|is\b|are\b|that\b|which\b|of\b|genes?\b|domains?\b)", after)) \
        and not after.lstrip().startswith("(")


def add_definitions(sentences):
    """The first use of each glossary term gets its definition in brackets, written by code."""
    seen = set()
    for item in sentences:
        text = item["text"]
        for term, pattern in TERM_PATTERNS.items():
            if term in seen:
                continue
            for match in pattern.finditer(text):
                if _stands_alone(text, match):
                    text = text[:match.end()] + " (" + GLOSSARY[term] + ")" + text[match.end():]
                    seen.add(term)
                    break
        item["shown_text"] = text
    return sentences


def result(analysis, facts):
    """What is stored and shown: verified sentences per part, the template where a part is empty, the hypothesis
    if it passed its check."""
    vocabulary = _vocabulary(facts)
    shown, checks, removed = {}, [], 0
    all_sentences = []
    for key, title in SECTIONS:
        kept = []
        for index, item in enumerate(analysis.get(key, [])):
            reasons = check_sentence(item, facts, vocabulary, key)
            checks.append({"section": key, "index": index, "text": item["text"], "verified": not reasons,
                           "reasons": reasons})
            if reasons:
                removed += 1
            else:
                kept.append(dict(item, source="model"))
        if not kept:
            kept = template_section(facts, key)
        shown[key] = kept
        all_sentences.extend(kept)
    add_definitions(all_sentences)
    hyp = analysis.get("hypothesis", {}) or {}
    hyp_reasons = check_hypothesis(hyp, facts)
    return {"facts": facts, "summary_sections": [{"key": k, "title": t, "sentences": shown[k]} for k, t in SECTIONS],
            "sentence_checks": checks, "sentences_removed": removed,
            "hypothesis": {"text": hyp.get("text", ""), "facts": hyp.get("facts", []), "shown": not hyp_reasons,
                           "reasons": hyp_reasons}}
