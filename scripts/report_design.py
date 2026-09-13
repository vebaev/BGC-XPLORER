from html import escape


PRIMARY_GLANCE_LABELS = (
    "Consensus",
    "Multi-tool",
    "High-confidence",
    "High-interest",
)


def reproducibility_panel(provenance):
    application = provenance.get("application", {})
    generation_model = str(provenance.get("ai", {}).get("model", "not-configured"))
    return (
        "<section class='panel'><div class='section-head'><h2>Reproducibility</h2>"
        "<span class='section-accent'></span></div>"
        "<p>Version: <strong>{version}</strong> · Commit: <code>{commit}</code> · "
        "ARTS reference: <strong>{arts}</strong></p>"
        "<p class='model-status'><span id='ai-model-label'>Active AI model:</span> "
        "<code id='active-ai-model' data-generation-model='{model}'>Checking active model...</code></p>"
        "</section>"
    ).format(
        version=escape(str(application.get("version", "unknown"))),
        commit=escape(str(application.get("git_commit", "unknown"))),
        arts=escape(str(provenance.get("arts_reference", "unknown"))),
        model=escape(generation_model, quote=True),
    )
