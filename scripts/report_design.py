from html import escape


PRIMARY_GLANCE_LABELS = (
    "Grouped loci",
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
)


def reproducibility_panel(provenance):
    application = provenance.get("application", {})
    generation_model = str(provenance.get("ai", {}).get("model", "not-configured"))
    return (
        "<section class='panel'><div class='section-head'><h2>Reproducibility</h2>"
        "<span class='section-accent'></span></div>"
        "<p>Version: <strong>{version}</strong> · Commit: <code>{commit}</code> · "
        "ARTS reference: <strong>{arts}</strong> · "
        "<span class='model-status'><span id='ai-model-label'>Active AI model:</span> "
        "<code id='active-ai-model' data-generation-model='{model}'>Checking active model...</code></span></p>"
        "</section>"
    ).format(
        version=escape(str(application.get("version", "unknown"))),
        commit=escape(str(application.get("git_commit", "unknown"))),
        arts=escape(str(provenance.get("arts_reference", "unknown"))),
        model=escape(generation_model, quote=True),
    )
