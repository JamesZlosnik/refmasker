"""Report writers: figures, HTML and PDF."""

from __future__ import annotations

import logging
from pathlib import Path

from refmasker.detectors import categories
from refmasker.report import html, pdf, plots

log = logging.getLogger(__name__)

FORMATS = {"html", "pdf", "png"}


def parse_formats(spec: str) -> set[str]:
    chosen = {f.strip().lower() for f in spec.split(",") if f.strip()}
    if "none" in chosen:
        return set()
    bad = chosen - FORMATS
    if bad:
        raise ValueError(
            f"Unknown report format(s): {', '.join(sorted(bad))}. Choose from html, pdf, png"
        )
    return chosen


def write_reports(data, outdir: Path, formats: set[str]) -> None:
    if not formats:
        log.info("No report formats requested")
        return
    labels = {c.name: c.label for c in categories().values()}
    figdir = outdir / "figures"
    figdir.mkdir(exist_ok=True)

    lin = plots.linear_figure(data, labels)
    circ = plots.circular_figure(data, labels)
    try:
        lin_paths = plots.save(lin, figdir / "linear", formats)
        circ_paths = plots.save(circ, figdir / "circular", formats)
        if "pdf" in formats:
            pdf.write_pdf(data, outdir / "report.pdf", [lin, circ], labels)
            log.info("Wrote %s", outdir / "report.pdf")
    finally:
        plots.close(lin)
        plots.close(circ)
    if "html" in formats:
        html.write_html(data, outdir / "report.html", [lin_paths["svg"], circ_paths["svg"]], labels)
        log.info("Wrote %s", outdir / "report.html")
