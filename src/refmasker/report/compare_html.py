"""Self-contained HTML report for `refmasker compare` (figures inlined, works offline)."""

from __future__ import annotations

import logging
from html import escape
from pathlib import Path

from refmasker import __version__
from refmasker.report.html import CSS, JS, strip_svg_header
from refmasker.report.plots import (
    COMPARE_COLOURS,
    close,
    compare_bar_figure,
    compare_tracks_figure,
    save_png,
    svg_string,
)
from refmasker.report.text import fmt_bp
from refmasker.safeio import open_replacing

log = logging.getLogger(__name__)

# Extra script for this page: the region tables share the generic sorter from the run report
# and add a text filter of their own.
COMPARE_JS = """
function filterTable(id) {
  const q = document.getElementById('q-' + id).value.toLowerCase();
  let shown = 0;
  document.querySelectorAll('#' + id + ' tbody tr').forEach(r => {
    const ok = !q || r.textContent.toLowerCase().includes(q);
    r.style.display = ok ? '' : 'none';
    if (ok) shown++;
  });
  document.getElementById('shown-' + id).textContent = shown;
}
"""

# How the two kinds of difference differ, spelled out because conflating them is the easy
# mistake to make when reading this page.
DIFFERENCE_NOTE = (
    "Two different questions are answered separately below. <em>Bases only in one mask</em> "
    "is a per-base difference: a region the other mask half covers counts half. "
    "<em>Regions the other mask misses entirely</em> counts whole regions that share no base "
    "with the other mask, which is what you want when asking “what would I lose by "
    "switching?”"
)


def _swatch(key: str) -> str:
    """A small square in one of the comparison colours, matching the figures."""
    return f'<span class="swatch" style="background:{COMPARE_COLOURS[key]}"></span>'


def _pct(value: float | None) -> str:
    """Percentages are only available with a reference FASTA; show a dash without one."""
    return "&mdash;" if value is None else f"{value:.2f}%"


def _tiles(c, e) -> str:
    """The headline numbers. Percentages of the genome appear only with a reference."""
    tiles = [
        ("Masked by A", f"{fmt_bp(c.bp_a)} bp", _pct(c.pct(c.bp_a))),
        ("Masked by B", f"{fmt_bp(c.bp_b)} bp", _pct(c.pct(c.bp_b))),
        ("Masked by both", f"{fmt_bp(c.bp_both)} bp", _pct(c.pct(c.bp_both))),
        ("Jaccard (bp)", f"{c.jaccard:.3f}", "1.000 = identical"),
        ("Of B, also in A", f"{100 * c.frac_b_covered:.1f}%", f"{fmt_bp(c.bp_only_b)} bp not in A"),
        ("Of A, also in B", f"{100 * c.frac_a_in_b:.1f}%", f"{fmt_bp(c.bp_only_a)} bp not in B"),
    ]
    out = ["<div class='tiles'>"]
    for key, value, sub in tiles:
        out.append(
            f"<div class='tile'><div class='k'>{e(key)}</div><div class='v'>{value}</div>"
            f"<div class='k' style='text-transform:none'>{sub}</div></div>"
        )
    out.append("</div>")
    return "".join(out)


def _record_table(c, e) -> str:
    """Per-sequence breakdown, longest sequence first."""
    known = bool(c.lengths)
    head = ["Sequence", "Length" if known else "Masked extent"]
    head += ["A (bp)", "B (bp)", "Both (bp)", "Only A (bp)", "Only B (bp)", "Jaccard"]
    if known:
        head += ["A (%)", "B (%)"]
    out = ["<table><thead><tr>"]
    for i, h in enumerate(head):
        cls = " class='n'" if i else ""
        out.append(f"<th{cls} onclick='sortTable(this)'>{e(h)}</th>")
    out.append("</tr></thead><tbody>")
    for r in c.per_record:
        cells = [
            f"<td>{e(r.record)}</td>",
            f"<td class='n' data-v='{r.span}'>{fmt_bp(r.span)}</td>",
            f"<td class='n' data-v='{r.bp_a}'>{fmt_bp(r.bp_a)}</td>",
            f"<td class='n' data-v='{r.bp_b}'>{fmt_bp(r.bp_b)}</td>",
            f"<td class='n' data-v='{r.bp_both}'>{fmt_bp(r.bp_both)}</td>",
            f"<td class='n' data-v='{r.bp_only_a}'>{fmt_bp(r.bp_only_a)}</td>",
            f"<td class='n' data-v='{r.bp_only_b}'>{fmt_bp(r.bp_only_b)}</td>",
            f"<td class='n' data-v='{r.jaccard}'>{r.jaccard:.3f}</td>",
        ]
        if known:
            cells += [
                f"<td class='n'>{_pct(r.pct(r.bp_a))}</td>",
                f"<td class='n'>{_pct(r.pct(r.bp_b))}</td>",
            ]
        out.append("<tr>" + "".join(cells) + "</tr>")
    out.append("</tbody></table>")
    return "".join(out)


def _region_table(regions, table_id: str, e) -> str:
    """A sortable, searchable table of merged regions one mask has and the other does not."""
    if not regions:
        return "<p>None: every region of this mask overlaps the other.</p>"
    out = [
        f"<div class='controls'><input id='q-{table_id}' placeholder='filter...' "
        f"oninput=\"filterTable('{table_id}')\">"
        f"<span class='sub'>showing <span id='shown-{table_id}'>{len(regions)}</span> "
        f"of {len(regions)}</span></div>",
        f"<div class='scroll'><table id='{table_id}'><thead><tr>",
        "<th onclick='sortTable(this)'>Sequence</th>",
        "<th class='n' onclick='sortTable(this)'>Start</th>",
        "<th class='n' onclick='sortTable(this)'>End</th>",
        "<th class='n' onclick='sortTable(this)'>Length (bp)</th>",
        "</tr></thead><tbody>",
    ]
    for r in sorted(regions, key=lambda x: -x.length):
        out.append(
            f"<tr><td>{e(r.chrom)}</td><td class='n' data-v='{r.start}'>{fmt_bp(r.start)}</td>"
            f"<td class='n' data-v='{r.end}'>{fmt_bp(r.end)}</td>"
            f"<td class='n' data-v='{r.length}'>{fmt_bp(r.length)}</td></tr>"
        )
    out.append("</tbody></table></div>")
    return "".join(out)


def write_compare_html(c, path: Path, meta: dict, formats: set[str]) -> None:
    """Write the comparison report as one self-contained page.

    `meta` carries what the page says about its own inputs: the two file names and paths, the
    reference if one was given, the renames applied, and when it was written. Every value
    from the inputs is HTML-escaped, so a sequence name cannot break the page. PNG copies of
    the figures are written beside the page when 'png' is among `formats`.
    """
    e = escape
    name_a, name_b = meta["name_a"], meta["name_b"]
    # Figures are rendered to strings, so a comparison leaves exactly one file behind unless
    # PNGs were asked for. An empty pair of masks has nothing to draw.
    svgs: list[str] = []
    if c.bp_a or c.bp_b:
        figures = {
            "bar": compare_bar_figure(c, name_a, name_b),
            "tracks": compare_tracks_figure(c, name_a, name_b),
        }
        try:
            svgs = [svg_string(fig) for fig in figures.values()]
            if "png" in formats:
                # PNG only: the page inlines the SVG from memory, so no SVG file is needed.
                for kind, fig in figures.items():
                    log.info("Wrote %s", save_png(fig, f"{path.with_suffix('')}.{kind}.png"))
        finally:
            for fig in figures.values():
                close(fig)

    parts: list[str] = []
    add = parts.append
    add("<!doctype html><html lang='en'><head><meta charset='utf-8'>")
    add("<meta name='viewport' content='width=device-width, initial-scale=1'>")
    add(f"<title>refmasker compare: {e(name_a)} vs {e(name_b)}</title>")
    add(f"<style>{CSS}</style></head><body><main>")
    add(f"<h1>Mask comparison: {e(name_a)} vs {e(name_b)}</h1>")
    add(
        f"<p class='sub'>refmasker {e(__version__)} &middot; {e(meta['written'])} &middot; "
        f"<strong>A</strong> = {e(name_a)} &middot; <strong>B</strong> = {e(name_b)}</p>"
    )

    for w in c.warnings:
        add(f"<div class='warn'><strong>Warning.</strong> {e(w)}</div>")
    if not c.lengths:
        add(
            "<div class='warn'>No reference FASTA was given (<code>--ref</code>), so the "
            "sequence lengths are unknown. Percentages of the genome are left out, and the "
            "figures run to the last masked base rather than to the end of each sequence.</div>"
        )

    add(_tiles(c, e))

    if svgs:
        add("<h2>How much the masks share</h2>")
        add(
            f"<p>{_swatch('both')}Bases both masks cover, {_swatch('only_a')}bases only A "
            f"covers, {_swatch('only_b')}bases only B covers. Each bar is one whole mask.</p>"
        )
        add(f"<div class='figure'>{strip_svg_header(svgs[0])}</div>")
        add("<h2>Where they agree and differ</h2>")
        add(
            "<p>One row per mask along each sequence. A stretch is green where both masks "
            "cover it and in the mask's own colour where only that mask does.</p>"
        )
        add(f"<div class='figure'>{strip_svg_header(svgs[1])}</div>")
    else:
        add("<h2>Figures</h2><p>Both masks are empty, so there is nothing to draw.</p>")

    add("<h2>By sequence</h2>")
    add(f"<p>{DIFFERENCE_NOTE}</p>")
    add(_record_table(c, e))

    add(f"<h2>Regions in A that B misses entirely ({len(c.only_a):,})</h2>")
    add(
        f"<p>Merged regions of {e(name_a)} that share no base with {e(name_b)}: "
        f"{fmt_bp(sum(r.length for r in c.only_a))} bp in total.</p>"
    )
    add(_region_table(c.only_a, "only-a", e))

    add(f"<h2>Regions in B that A misses entirely ({len(c.only_b):,})</h2>")
    add(
        f"<p>Merged regions of {e(name_b)} that share no base with {e(name_a)}: "
        f"{fmt_bp(sum(r.length for r in c.only_b))} bp in total.</p>"
    )
    add(_region_table(c.only_b, "only-b", e))

    # What produced the page, so a saved report still says what it compared.
    add("<h2>Run details</h2><dl>")
    rows = [
        ("Mask A", meta["path_a"]),
        ("Mask B", meta["path_b"]),
        ("Reference", meta.get("ref") or "not given (--ref)"),
        ("Renames", meta.get("rename") or "none"),
        ("Merged regions", f"A {c.regions_a:,} &middot; B {c.regions_b:,}"),
        ("Written", meta["written"]),
        ("refmasker", __version__),
    ]
    for k, v in rows:
        add(f"<dt>{e(k)}</dt><dd>{e(str(v))}</dd>")
    add("</dl>")

    add(f"<script>{JS}{COMPARE_JS}</script></main></body></html>")
    with open_replacing(path) as fh:
        fh.write("".join(parts))
