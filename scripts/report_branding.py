import base64
import mimetypes
from html import escape
from pathlib import Path


def image_data_uri(candidates):
    for candidate in candidates:
        path = Path(candidate)
        if not path.is_file():
            continue
        media_type = mimetypes.guess_type(str(path))[0] or "application/octet-stream"
        encoded = base64.b64encode(path.read_bytes()).decode("ascii")
        return "data:{0};base64,{1}".format(media_type, encoded)
    return ""


def report_home_link(home_url="/"):
    return (
        "<nav class='report-nav' aria-label='Report navigation'>"
        "<a class='report-home-link' href='{url}' target='_top'>"
        "<span aria-hidden='true'>&larr;</span> Back to BGC-XPLORER"
        "</a></nav>"
    ).format(url=escape(home_url, quote=True))
