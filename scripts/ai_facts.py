"""Fact-card summaries: the facts are chosen by code, the model only words them.

The model is never given the whole locus record. Code selects the handful of
facts that matter for a biosynthetic gene cluster locus (where it is, which
callers support it, its class, its core genes, its closest MIBiG entry, ARTS
hits, flags and what is missing) and states each as a short sentence with the
values it must keep. The model rewrites those sentences into one readable
paragraph. The paragraph is accepted only when every fact's values are in it,
no number, identifier or caller name appears that the facts do not contain,
each core gene appears in one sentence with its callers, and it uses no
interpretive wording. Otherwise the report shows the facts as written by the
template, so text that failed the check is never shown as the summary.
"""
import csv
import os
import re

import ai_evidence
from bgc_classes import harmonised_classes

PROMPT_VERSION = "4.1-facts"
MAX_CORE_GENES = 5
EDGE_DISTANCE_BP = 10000
CALLER_NAMES = {"antismash": "antiSMASH", "gecco": "GECCO", "deepbgc": "DeepBGC"}

SYSTEM_PROMPT = """You turn a short list of facts about one biosynthetic gene cluster locus into one readable paragraph.
Rules:
- Use every fact, and nothing else. Do not add any information, interpretation, function, product, activity or
  comparison that is not in the facts.
- Combine the facts into a few well-formed sentences instead of repeating them one per line: for example put the
  location, callers and class in one sentence, and core genes named by the same callers in one sentence
  ("CJLEIP_01147 and CJLEIP_01149 (both D-alanine--poly(phosphoribitol) ligase subunit DltA) are core genes named
  by antiSMASH, DeepBGC and GECCO").
- Keep every number, locus tag, accession, compound name, class and caller name exactly as written in the facts.
  Do not round or reformat numbers, and do not add numbers of your own.
- Name the callers right after the core gene or group of core genes they name, in the same sentence; never move a
  caller to another gene, and keep each gene's product next to the gene.
- Do not use words such as likely, suggests, may, might, probably, possibly, could, indicates, appears,
  putative, produces, synthesises, responsible, novel, antibiotic, activity.
- 3 to 6 sentences, plain scientific English.
Return ONLY JSON: {"summary": "..."}"""

JSON_SCHEMA = {
    "type": "object",
    "additionalProperties": False,
    "properties": {"summary": {"type": "string", "minLength": 40}},
    "required": ["summary"],
}

FORBIDDEN = re.compile(
    r"\b(likely|suggests?|suggesting|may|might|probably|possibly|could|indicates?|indicating|appears?|"
    r"presumably|putative|produces?|producing|synthesi[sz]es?|responsible|novel|antibiotics?|activity)\b",
    re.IGNORECASE)
CALLER_PATTERN = re.compile(r"\b(antismash|gecco|deepbgc)\b", re.IGNORECASE)


def _clean(value):
    text = "" if value is None else str(value).strip()
    return "" if text.lower() in {"", "nan", "none", "-"} else text


def _int(value):
    try:
        return int(float(str(value).replace(",", "")))
    except ValueError:
        return None


def _number(value):
    """A number as the facts will print it: integers without a trailing .0."""
    n = ai_evidence._as_number(_clean(value))
    if n is None:
        return _clean(value)
    return str(int(n)) if n == int(n) else str(n)


def _plural(n, word, plural=None):
    return "{0} {1}".format(n, word if str(n) in ("1", "1.0") else (plural or word + "s"))


def _callers(text):
    return [c.strip().lower() for c in _clean(text).split(",") if c.strip()]


def _caller_list(callers):
    names = [CALLER_NAMES.get(c, c) for c in callers]
    return names[0] if len(names) == 1 else ", ".join(names[:-1]) + " and " + names[-1]


def _core_support(text):
    support = {}
    for part in _clean(text).split(";"):
        if ":" in part:
            caller, tags = part.split(":", 1)
            for tag in tags.split(","):
                if tag.strip():
                    support.setdefault(tag.strip(), []).append(caller.strip().lower())
    return support


def fact(key, sentence, tokens, together=(), callers=()):
    """One fact: its sentence, the values the paragraph must keep, the values that must share a sentence,
    and, for a core gene, the callers that name it."""
    return {"key": key, "sentence": sentence, "tokens": [t for t in tokens if t],
            "together": [t for t in together if t], "callers": list(callers)}


def load_conflicts(results_dir, sample, consensus_id):
    """Locus tags flagged with an annotation conflict in cluster_genes.tsv for one locus."""
    path = os.path.join(str(results_dir), sample, "summary", "cluster_genes.tsv")
    if not os.path.exists(path):
        return []
    with open(path, newline="") as handle:
        return sorted({row.get("locus_tag", "") for row in csv.DictReader(handle, delimiter="\t")
                       if row.get("consensus_id") == consensus_id and _clean(row.get("annotation_conflict"))})


def _genes_in_span(genes, start, end):
    """Genes whose midpoint lies inside the reported (agreed) span."""
    lo, hi = _int(start), _int(end)
    count = 0
    for gene in genes:
        a, b = _int(gene.get("gene_start") or 0), _int(gene.get("gene_end") or 0)
        if lo is not None and hi is not None and a is not None and b is not None and lo <= (a + b) / 2 <= hi:
            count += 1
    return count


def build_fact_card(payload, conflicts=()):
    """The facts for one locus, from the full payload of build_cluster_payload."""
    region = payload.get("region", {}) or {}
    genes = payload.get("genes", []) or []
    facts = []

    contig, start, end = _clean(region.get("contig")), _number(region.get("start")), _number(region.get("end"))
    length = _number(region.get("length_bp"))
    in_span = _genes_in_span(genes, start, end)
    facts.append(fact("location", "The locus spans {0}:{1}-{2} ({3} bp, {4} genes).".format(
        contig, start, end, length, in_span), [contig, start, end, length, str(in_span)]))

    callers = _callers(region.get("support_tools"))
    if callers:
        facts.append(fact("callers", "It is supported by {0} ({1} of 3 callers).".format(
            _caller_list(callers), len(callers)), [CALLER_NAMES.get(c, c) for c in callers] + [str(len(callers))]))

    classes = [c for c in harmonised_classes(region.get("bgc_types")) if c != "Unknown"]
    if classes:
        facts.append(fact("class", "Its MIBiG compound class{0} {1} {2}.".format(
            "es" if len(classes) > 1 else "", "are" if len(classes) > 1 else "is", _caller_list(classes)), classes))
    else:
        facts.append(fact("class", "No caller assigned a compound class.", ["class"]))

    support = _core_support(region.get("core_gene_support"))
    products = {_clean(g.get("locus_tag")): _clean(g.get("bakta_product")) for g in genes}
    midpoints = {_clean(g.get("locus_tag")): ((_int(g.get("gene_start") or 0) or 0) + (_int(g.get("gene_end") or 0) or 0)) / 2
                 for g in genes}
    lo, hi = _int(start), _int(end)
    if support:
        tags = list(support)
        for tag in tags[:MAX_CORE_GENES]:
            names = [CALLER_NAMES.get(c, c) for c in support[tag]]
            product = products.get(tag, "")
            mid = midpoints.get(tag)
            outside = mid is not None and lo is not None and hi is not None and not lo <= mid <= hi
            facts.append(fact("core:" + tag, "{0}{1} is a core gene named by {2}{3}.".format(
                tag, " ({0})".format(product) if product else "", _caller_list(support[tag]),
                ", outside the agreed span" if outside else ""),
                [tag, product] + names, together=[tag] + names, callers=support[tag]))
        if len(tags) > MAX_CORE_GENES:
            rest = len(tags) - MAX_CORE_GENES
            facts.append(fact("core:more", "{0} further core gene{1} listed in the evidence table.".format(
                rest, " is" if rest == 1 else "s are"), [str(rest)]))
    else:
        facts.append(fact("core:none", "No caller marks a biosynthetic core gene in this locus.", ["core"]))

    mibig = payload.get("mibig_dereplication", {}) or {}
    accession = _clean(mibig.get("best_mibig_id") or region.get("best_mibig_id"))
    if accession:
        compound = _clean(mibig.get("best_mibig_product") or region.get("best_mibig_product"))
        first_compound = compound.split("/")[0].strip()
        metric = _clean(mibig.get("score_metric") or region.get("score_metric"))
        score = _number(mibig.get("match_score") or region.get("match_score"))
        matched = _number(mibig.get("matched_genes") or region.get("matched_genes"))
        core_hits = _number(mibig.get("core_gene_hits") or region.get("core_gene_hits"))
        if score and ai_evidence._as_number(score) is not None and float(score) < 1:
            score = "{0:.2f}".format(float(score))  # similarity fractions, shown to two decimals
        parts = ["The closest MIBiG entry is {0}{1}".format(accession, " ({0})".format(compound) if compound else "")]
        detail = []
        if metric and score:
            detail.append("{0} {1}".format(metric, score))
        if matched:
            detail.append(_plural(matched, "matched gene") + (", {0} of them core".format(core_hits) if core_hits else ""))
        sentence = parts[0] + (": " + ", ".join(detail) if detail else "") + "."
        facts.append(fact("mibig", sentence, [accession, first_compound, score, matched, core_hits],
                          together=[accession, first_compound, score]))
    else:
        facts.append(fact("mibig", "No MIBiG comparison is reported.", ["MIBiG"]))

    known, duf = _int(region.get("arts_known_hits")) or 0, _int(region.get("arts_duf_hits")) or 0
    if known or duf:
        facts.append(fact("arts", "ARTS reports {0} and {1}.".format(
            _plural(known, "known resistance hit"), _plural(duf, "DUF hit")),
                          ["ARTS", str(known), str(duf)], together=["ARTS", str(known), str(duf)]))
    else:
        facts.append(fact("arts", "ARTS reports no hits.", ["ARTS"]))

    if conflicts:
        facts.append(fact("conflict", "{0} gene{1} carr{2} an annotation conflict (a dbCAN family on a biosynthetic core gene): {3}.".format(
            len(conflicts), "" if len(conflicts) == 1 else "s", "ies" if len(conflicts) == 1 else "y", ", ".join(conflicts)),
            ["conflict", str(len(conflicts))] + list(conflicts), together=["conflict", str(len(conflicts))]))

    edge = _int(region.get("nearest_contig_edge_bp"))
    if edge is not None and edge < EDGE_DISTANCE_BP:
        facts.append(fact("edge", "The locus lies {0} bp from a contig edge.".format(edge), ["edge", str(edge)],
                          together=["edge", str(edge)]))
    return facts


def template_summary(facts):
    return " ".join(f["sentence"] for f in facts)


def build_user_prompt(facts):
    return "Facts:\n" + "\n".join("- " + f["sentence"] for f in facts)


NUMBER_WORDS = {w: str(i) for i, w in enumerate(
    "zero one two three four five six seven eight nine ten eleven twelve thirteen fourteen fifteen sixteen "
    "seventeen eighteen nineteen twenty".split())}


def _digits(text):
    """Numbers written as words become digits, so 'Six genes' counts as 6 genes."""
    return re.sub(r"\b(" + "|".join(NUMBER_WORDS) + r")\b", lambda m: NUMBER_WORDS[m.group(1).lower()],
                  text, flags=re.IGNORECASE)


def _contains(text, token):
    """A number must occur as a number (never inside a locus tag); words are matched case-insensitively."""
    n = ai_evidence._as_number(token)
    if n is not None:
        return n in set(ai_evidence._numbers(_digits(text)))
    return token.lower() in text.lower()


def verify(summary, facts):
    """Accept the paragraph only if it keeps every fact and adds nothing checkable."""
    reasons = []
    text = _clean(summary)
    if not text:
        return {"verified": False, "reasons": ["empty summary"], "facts": len(facts), "facts_kept": 0}
    if FORBIDDEN.search(text):
        reasons.append("interpretive wording: '{0}'".format(FORBIDDEN.search(text).group(0)))
    kept = 0
    for f in facts:
        missing = [t for t in f["tokens"] if not _contains(text, t)]
        if missing:
            reasons.append("fact '{0}' lost: {1}".format(f["key"], ", ".join(missing)))
        else:
            kept += 1
    tokens = [t for f in facts for t in f["tokens"]] + [f["sentence"] for f in facts]
    pool = " ".join(tokens)
    allowed_numbers = set(ai_evidence._numbers(pool)) | {3.0}
    for number in ai_evidence._numbers(_digits(text)):
        if number not in allowed_numbers:
            reasons.append("the number {0:g} is not in the facts".format(number))
    for identifier in set(ai_evidence.IDENTIFIER.findall(text)):
        if identifier not in pool:
            reasons.append("{0} is not in the facts".format(identifier))
    allowed_callers = {c.lower() for c in CALLER_PATTERN.findall(pool)}
    for caller in {c.lower() for c in CALLER_PATTERN.findall(text)}:
        if caller not in allowed_callers:
            reasons.append("{0} is not among the callers in the facts".format(caller))
    sentences = re.split(r"(?<=[.;])\s+", text)
    groups = caller_groups(text)
    for f in facts:
        if f["callers"]:
            reasons.extend(_check_gene_callers(text, f, groups))
        elif f["together"] and not any(all(_contains(s, t) for t in f["together"]) for s in sentences):
            reasons.append("fact '{0}': {1} are not in one sentence".format(f["key"], ", ".join(f["together"])))
    return {"verified": not reasons, "reasons": reasons, "facts": len(facts), "facts_kept": kept}


def caller_groups(text):
    """Map each locus tag to the callers named after it in its sentence.

    A run of locus tags followed by caller names gives those callers to every
    tag in the run: "A (x) and B (y) are core genes named by GECCO" names both.
    """
    assigned = {}
    for sentence in re.split(r"(?<=[.;])\s+", text):
        pending, callers = [], None
        for match in re.finditer(ai_evidence.LOCUS_TAG + r"|\b(?:antismash|gecco|deepbgc)\b", sentence, re.IGNORECASE):
            token = match.group(0)
            if CALLER_PATTERN.fullmatch(token):
                if pending:
                    callers = set() if callers is None else callers
                    callers.add(token.lower())
            else:
                if callers:
                    for tag in pending:
                        assigned.setdefault(tag, []).append(callers)
                    pending, callers = [], None
                pending.append(token)
        if pending and callers:
            for tag in pending:
                assigned.setdefault(tag, []).append(callers)
    return assigned


def _check_gene_callers(text, f, groups=None):
    """The callers named after a core gene (or its group of genes) must be exactly the callers that name it."""
    tag = f["together"][0]
    groups = caller_groups(text) if groups is None else groups
    if tag not in groups:
        return ["{0} is not followed by the callers that name it".format(tag)]
    return ["{0} is said to be named by {1}, but is named by {2}".format(
        tag, ", ".join(sorted(named)), ", ".join(sorted(f["callers"])))
        for named in groups[tag] if named != set(f["callers"])]


def normalize(analysis):
    if isinstance(analysis, dict):
        return {"summary": _clean(analysis.get("summary"))}
    return {"summary": _clean(analysis)}


def result(analysis, facts):
    """What is stored and shown: the model's paragraph if verified, otherwise the template."""
    check = verify(analysis.get("summary", ""), facts)
    shown = analysis.get("summary") if check["verified"] else template_summary(facts)
    return {"facts": facts, "fact_check": check, "shown_summary": shown,
            "shown_source": "model" if check["verified"] else "template"}
