from html import escape


PRIMARY_GLANCE_LABELS = (
    "Total loci",
    "Multi-caller loci",
    "MIBiG comparisons",
    "ARTS known hits",
)

DONUT_COLORS = (
    "#b8a9e8",
    "#a8d8c7",
    "#a9cce8",
    "#f3c6a8",
    "#e8b4c4",
    "#9fd6d2",
    "#e3d7a6",
    "#cfd4dc",
)


def bakta_database_label(bakta_metadata):
    """'full 6.0' or 'light 6.0': the Bakta database type and version the run used, from Bakta's own output."""
    db = ((bakta_metadata or {}).get("bakta", {}).get("version", {}) or {}).get("db", {}) or {}
    label = " ".join(str(db[k]) for k in ("type", "version") if db.get(k))
    return label or "unknown"


def reproducibility_panel(provenance, bakta_metadata=None):
    application = provenance.get("application", {})
    generation_model = str(provenance.get("ai", {}).get("model", "not-configured"))
    return (
        "<section class='panel'><div class='section-head'><h2>Reproducibility</h2>"
        "<span class='section-accent'></span></div>"
        "<p>Version: <strong>{version}</strong> · Commit: <code>{commit}</code> · "
        "Bakta database: <strong>{bakta_db}</strong> · "
        "ARTS reference: <strong>{arts}</strong> · "
        "<span class='model-status'><span id='ai-model-label'>Active AI model:</span> "
        "<code id='active-ai-model' data-generation-model='{model}'>Checking active model...</code></span></p>"
        "</section>"
    ).format(
        version=escape(str(application.get("version", "unknown"))),
        commit=escape(str(application.get("git_commit", "unknown"))),
        bakta_db=escape(bakta_database_label(bakta_metadata)),
        arts=escape(str(provenance.get("arts_reference", "unknown"))),
        model=escape(generation_model, quote=True),
    )
