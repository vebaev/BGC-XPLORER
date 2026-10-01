"""One vocabulary for the class and activity labels the three callers report.

antiSMASH, GECCO and DeepBGC name the same chemistry differently (NRPS, NRP,
NRPS-like; T1PKS, Polyketide), and DeepBGC adds an activity prediction
(antibacterial, cytotoxic) to the same product field. Counting the raw labels
splits one class across several names and mixes activity with class, so the
report maps classes onto the MIBiG compound classes and keeps activity apart.
"""

ACTIVITY_TERMS = {
    "antibacterial", "antifungal", "cytotoxic", "inhibitor",
    "antibacterial-cytotoxic", "antibacterial-inhibitor",
}
ACTIVITIES = ("antibacterial", "antifungal", "cytotoxic", "inhibitor")

MIBIG_CLASSES = ("NRP", "Polyketide", "RiPP", "Terpene", "Saccharide", "Alkaloid", "Other", "Unknown")
_BASE = {name.lower(): name for name in MIBIG_CLASSES}

RIPP_TERMS = {
    "lap", "bottromycin", "cyanobactin", "glycocin", "linaridin", "proteusin", "rre-containing",
    "microviridin", "guanidinotides", "spliceotide", "crocagin", "darobactin", "methanobactin",
    "triceptide", "lipolanthine", "thioamitides", "redox-cofactor", "head_to_tail",
    "cyclic-lactone-autoinducer",
}
SACCHARIDE_TERMS = {"amglyccycl", "aminoglycoside", "oligosaccharide", "saccharide"}
POLYKETIDE_TERMS = {"polyketide", "hgle-ks", "arylpolyene", "ladderane"}


def split_labels(value):
    """Comma- or semicolon-separated labels, without blanks and NaN."""
    text = "" if value is None else str(value)
    parts = []
    for chunk in text.replace(";", ",").split(","):
        label = chunk.strip()
        if label and label.lower() != "nan":
            parts.append(label)
    return parts


def mibig_class(label):
    """The MIBiG class for one caller label, or None for an activity label."""
    term = str(label).strip().lower()
    if not term or term == "unknown":
        return "Unknown"
    if term in ACTIVITY_TERMS:
        return None
    if term in _BASE:
        return _BASE[term]
    if term.startswith("nrp") or term.endswith("-nrp"):
        return "NRP"
    if "pks" in term or term in POLYKETIDE_TERMS or term.startswith("enediyne"):
        return "Polyketide"
    if "terpene" in term:
        return "Terpene"
    if "ripp" in term or "lanthipeptide" in term or term.endswith("peptide") or term in RIPP_TERMS:
        return "RiPP"
    if term in SACCHARIDE_TERMS:
        return "Saccharide"
    if "alkaloid" in term:
        return "Alkaloid"
    return "Other"


def harmonised_classes(value):
    """MIBiG classes behind a label list, in first-seen order, each once.

    A hybrid written as "NRP-Polyketide" counts for both classes.
    """
    classes = []
    for label in split_labels(value):
        pieces = label.split("-")
        if len(pieces) > 1 and all(piece.strip().lower() in _BASE for piece in pieces):
            names = [_BASE[piece.strip().lower()] for piece in pieces]
        else:
            names = [mibig_class(label)]
        for name in names:
            if name and name not in classes:
                classes.append(name)
    return classes


def activity_labels(value):
    """Predicted activities in a label list; "antibacterial-cytotoxic" gives both."""
    found = []
    for label in split_labels(value):
        term = label.lower()
        if term not in ACTIVITY_TERMS:
            continue
        for part in term.split("-"):
            if part in ACTIVITIES and part not in found:
                found.append(part)
    return found


def count_per_locus(series, labels_of, limit=None):
    """(label, number of loci) pairs, most frequent first; a locus counts once per label."""
    counts = {}
    for value in series.fillna(""):
        for label in labels_of(value):
            counts[label] = counts.get(label, 0) + 1
    ordered = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    return ordered[:limit] if limit else ordered
