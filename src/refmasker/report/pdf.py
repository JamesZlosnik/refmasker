"""PDF report built with matplotlib (no extra dependencies)."""

from __future__ import annotations

import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from refmasker.report.text import SNIPPY_USAGE, WHY
from refmasker.safeio import open_replacing

A4 = (8.27, 11.69)  # page size in inches (matplotlib figure sizes are in inches)
ROWS_PER_PAGE = 48  # region table rows that fit on one page


def _text_page(pdf, lines: list[tuple[str, dict]]):
    """Draw lines of text down one A4 page.

    Each line has a style dict with optional keys: size (pt), weight, family, wrap
    (characters per line) and after (extra space below, as a fraction of the page height).
    """
    fig = plt.figure(figsize=A4)
    y = 0.95  # vertical position as a fraction of the page height, from the bottom
    for text, style in lines:
        size = style.get("size", 9)
        wrap = style.get("wrap", 100)
        for chunk in textwrap.wrap(text, wrap) or [""]:
            fig.text(
                0.07,
                y,
                chunk,
                fontsize=size,
                weight=style.get("weight", "normal"),
                family=style.get("family", "sans-serif"),
                va="top",
            )
            # Move down one line: points / 72 = inches, / page height = fraction of the
            # page, x 1.45 for line spacing.
            y -= size / 72 / A4[1] * 1.45
        y -= style.get("after", 0.004)
    pdf.savefig(fig)
    plt.close(fig)


def _table_page(pdf, title: str, header: list[str], rows: list[list[str]], widths=None, note=""):
    """Draw one table on an A4 page, with a title, an optional note and a bold header row."""
    fig = plt.figure(figsize=A4)
    fig.text(0.07, 0.96, title, fontsize=12, weight="bold", va="top")
    if note:
        fig.text(0.07, 0.935, note, fontsize=8, va="top", color="#555555")
    ax = fig.add_axes((0.07, 0.05, 0.86, 0.86))
    ax.axis("off")
    if rows:
        tbl = ax.table(
            cellText=rows, colLabels=header, loc="upper left", colWidths=widths, cellLoc="left"
        )
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(7)
        tbl.scale(1, 1.15)
        for (r, _c), cell in tbl.get_celld().items():
            cell.set_linewidth(0.3)
            if r == 0:
                cell.set_facecolor("#eeeeee")
                cell.set_text_props(weight="bold")
    pdf.savefig(fig)
    plt.close(fig)


def write_pdf(data, path: Path, figures, labels: dict[str, str]) -> None:
    """Write report.pdf with matplotlib alone (no PDF library needed).

    Pages: summary text, the figures, the repeat sensitivity table, the region table (split
    over as many pages as needed) and run details.
    """
    total = next(r for r in data.summary if r.record == "ALL" and r.category == "total")
    with open_replacing(path, "wb") as fh, PdfPages(fh) as pdf:
        # Page 1: headline numbers, warnings, per-category amounts, reasons and snippy usage.
        lines: list[tuple[str, dict]] = [
            (
                f"Reference mask for {data.reference.name}",
                {"size": 16, "weight": "bold", "after": 0.01},
            ),
            (
                f"refmasker {data.versions.get('refmasker', '')} · {data.started} · "
                f"built for {data.params['read_length']} bp reads",
                {"size": 9, "after": 0.02},
            ),
            (
                f"Masked {total.bp:,} of {total.length:,} bp ({total.pct:.2f}%) in "
                f"{total.regions:,} regions across {len(data.records)} record(s).",
                {"size": 11, "weight": "bold", "after": 0.015},
            ),
        ]
        for w in data.warnings:
            lines.append((f"Warning: {w}", {"size": 9, "after": 0.01}))
        lines.append(("Masked per category", {"size": 12, "weight": "bold", "after": 0.006}))
        for r in (x for x in data.summary if x.record == "ALL"):
            name = (
                "Total (each base counted once)"
                if r.category == "total"
                else labels.get(r.category, r.category)
            )
            lines.append(
                (f"{name}: {r.bp:,} bp ({r.pct:.2f}%), {r.regions:,} regions", {"size": 9})
            )
        lines.append(("", {"after": 0.01}))
        lines.append(
            ("Why each category is masked", {"size": 12, "weight": "bold", "after": 0.006})
        )
        for cat in data.categories:
            lines.append(
                (f"{labels.get(cat, cat)}: {WHY.get(cat, '')}", {"size": 8.5, "after": 0.006})
            )
        lines.append(("", {"after": 0.01}))
        lines.append(("Using the mask", {"size": 12, "weight": "bold", "after": 0.006}))
        lines.append((SNIPPY_USAGE, {"size": 8, "family": "monospace", "wrap": 110}))
        _text_page(pdf, lines)

        # One page per figure (linear, then circular).
        for fig in figures:
            pdf.savefig(fig, bbox_inches="tight")

        if data.sensitivity:
            _table_page(
                pdf,
                "Repeat threshold sensitivity",
                ["Min length (bp)", "Min identity (%)", "Regions", "bp masked", "% genome", "Used"],
                [
                    [
                        str(s.min_len),
                        f"{s.min_id:g}",
                        f"{s.regions:,}",
                        f"{s.bp:,}",
                        f"{s.pct:.2f}",
                        "yes" if s.chosen else "",
                    ]
                    for s in data.sensitivity
                ],
            )

        # Region table. Long names and gene lists are cut short to fit the page.
        rows = [
            [
                r.chrom,
                f"{r.start + 1:,}",
                f"{r.end:,}",
                f"{r.length:,}",
                labels.get(r.category, r.category),
                r.source,
                r.name[:38],
                (", ".join(r.genes))[:48],
            ]
            for r in data.regions
        ]
        header = ["Record", "Start", "End", "Length", "Category", "Source", "Name", "Genes"]
        widths = [0.13, 0.08, 0.08, 0.07, 0.1, 0.09, 0.2, 0.25]
        pages = max(1, -(-len(rows) // ROWS_PER_PAGE))  # -(-a // b) rounds a / b up
        for p in range(pages):
            chunk = rows[p * ROWS_PER_PAGE : (p + 1) * ROWS_PER_PAGE]
            _table_page(
                pdf,
                f"Regions ({p + 1}/{pages})",
                header,
                chunk,
                widths,
                note="Start is 1-based; names and gene lists are truncated (full list in mask.annotated.bed).",
            )

        # Last page: run details for provenance.
        details = [
            ("Run details", {"size": 12, "weight": "bold", "after": 0.008}),
            (f"Reference: {data.reference}", {"size": 8.5}),
            (f"md5 {data.reference_md5}  sha256 {data.reference_sha256}", {"size": 7.5}),
            (f"Annotation: {data.annotation or 'none'}", {"size": 8.5}),
            (f"Plasmid records: {', '.join(data.plasmids) or 'none'}", {"size": 8.5}),
            (f"Command: {data.command}", {"size": 8.5}),
            (f"Runtime: {data.runtime_s:.0f} s", {"size": 8.5, "after": 0.01}),
        ]
        details += [(f"{k}: {v}", {"size": 8.5}) for k, v in data.params.items() if v is not None]
        details.append(("", {"after": 0.01}))
        details += [(f"{k} {v}", {"size": 8.5}) for k, v in data.versions.items()]
        _text_page(pdf, details)
