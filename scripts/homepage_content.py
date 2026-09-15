from html import escape


HERO_TITLE = "Explore biosynthetic gene clusters"
HERO_BODY = (
    "BGC-XPLORER turns an assembled bacterial genome into an integrated, "
    "reproducible analysis of biosynthetic potential. It combines multi-tool "
    "BGC discovery with functional and resistance context, searchable evidence, "
    "gene maps, and optional AI-supported interpretation."
)

WORKFLOW_STEPS = (
    ("Genome input", "Upload an assembled bacterial genome, MAG, or plasmid as nucleotide FASTA.", "1"),
    ("Genome annotation", "Bakta prepares the standardized genes and proteins used downstream.", "2"),
    ("BGC discovery", "antiSMASH, GECCO, and DeepBGC independently identify candidate regions.", "3"),
    ("Biological context", "eggNOG, dbCAN, ARTS, and MIBiG add function, substrates, resistance, and known-cluster evidence.", "4"),
    ("Integrated report", "Compare, filter, and sort candidate loci with gene maps, provenance, and optional AI analysis.", "5"),
)

RESULT_FEATURES = (
    ("Grouped candidate loci", "Overlapping calls are grouped while each original caller signal remains visible.", "◎"),
    ("Searchable evidence", "Filter and sort loci by caller, ARTS signal, MIBiG comparison, location, and class.", "✦"),
    ("Biological context", "Resistance signals, functional annotations, substrates, and MIBiG evidence are connected to each region.", "⌘"),
    ("Interactive gene maps", "Explore gene order, predicted functions, and cluster-level biological interpretation.", "↔"),
    ("Reproducible outputs", "Download tables and provenance with tool, database, model, and execution metadata.", "✓"),
    ("AI-supported interpretation", "Request focused AI analysis for selected clusters directly from the report.", "◇"),
)


def section_header(title, description=""):
    markup = (
        "<div class='section-head'><h2>{title}</h2>"
        "<span class='section-accent'></span></div>"
    ).format(title=escape(title))
    if description:
        markup += "<p class='section-description muted'>{0}</p>".format(escape(description))
    return markup
