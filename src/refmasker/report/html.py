"""Self-contained HTML report (works offline; figures are inlined SVG)."""

from __future__ import annotations

import re
from html import escape
from pathlib import Path

from refmasker.report.palette import colour
from refmasker.report.text import SNIPPY_USAGE, WHY, fmt_bp
from refmasker.safeio import open_replacing

# Page styles. Dark mode uses the prefers-color-scheme media query; figures keep a white
# background there because matplotlib draws them with dark text.
CSS = """
:root { --bg:#ffffff; --fg:#1d1d1f; --muted:#666; --line:#e3e3e3; --panel:#f6f7f8;
        --warn-bg:#fff6e0; --warn-line:#e0b040; --accent:#0072B2; }
@media (prefers-color-scheme: dark) {
  :root { --bg:#16181b; --fg:#e8e8e8; --muted:#a0a0a0; --line:#33363b; --panel:#1f2226;
          --warn-bg:#3a3020; --warn-line:#a07a20; --accent:#5aa9e6; }
  .figure svg { background:#ffffff; border-radius:6px; }
}
* { box-sizing:border-box; }
body { margin:0; background:var(--bg); color:var(--fg);
       font:14px/1.5 -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; }
main { max-width:1150px; margin:0 auto; padding:24px 16px 64px; }
h1 { font-size:24px; margin:0 0 4px; }
h2 { font-size:18px; margin:36px 0 10px; padding-bottom:4px; border-bottom:1px solid var(--line); }
.sub { color:var(--muted); margin:0 0 20px; }
.tiles { display:grid; grid-template-columns:repeat(auto-fit,minmax(170px,1fr)); gap:12px; }
.tile { background:var(--panel); border-radius:8px; padding:12px 14px; }
.tile .v { font-size:22px; font-weight:600; }
.tile .k { color:var(--muted); font-size:12px; text-transform:uppercase; letter-spacing:.04em; }
.warn { background:var(--warn-bg); border-left:4px solid var(--warn-line); padding:10px 14px;
        border-radius:4px; margin:12px 0; }
table { border-collapse:collapse; width:100%; font-size:13px; }
th, td { text-align:left; padding:5px 8px; border-bottom:1px solid var(--line); vertical-align:top; }
th { background:var(--panel); position:sticky; top:0; cursor:pointer; user-select:none; }
td.n, th.n { text-align:right; font-variant-numeric:tabular-nums; }
tr.chosen td { font-weight:600; background:var(--panel); }
.swatch { display:inline-block; width:11px; height:11px; border-radius:2px; margin-right:6px;
          vertical-align:-1px; }
.figure { overflow-x:auto; margin:8px 0 16px; text-align:center; }
.figure svg { max-width:100%; height:auto; }
.scroll { max-height:560px; overflow:auto; border:1px solid var(--line); border-radius:6px; }
.controls { display:flex; gap:10px; flex-wrap:wrap; margin:6px 0 10px; }
input, select { font:inherit; padding:4px 8px; border:1px solid var(--line); border-radius:4px;
                background:var(--bg); color:var(--fg); }
code, pre { font-family:ui-monospace, Menlo, Consolas, monospace; font-size:12.5px; }
pre { background:var(--panel); padding:10px 12px; border-radius:6px; overflow-x:auto; }
.genes { color:var(--muted); max-width:380px; }
dl { display:grid; grid-template-columns:max-content 1fr; gap:4px 16px; margin:0; }
dt { color:var(--muted); }
dd { margin:0; word-break:break-all; }
"""

# Scripts for the region table: click a column heading to sort, and filter rows by text or
# category. Numeric cells sort by their data-v attribute (the raw number), not the shown text.
JS = """
function sortTable(th) {
  const table = th.closest('table'), tbody = table.tBodies[0];
  const idx = Array.from(th.parentNode.children).indexOf(th);
  const numeric = th.classList.contains('n');
  const asc = th.dataset.dir !== 'asc';
  th.parentNode.querySelectorAll('th').forEach(h => delete h.dataset.dir);
  th.dataset.dir = asc ? 'asc' : 'desc';
  const rows = Array.from(tbody.rows);
  rows.sort((a, b) => {
    let x = a.cells[idx].dataset.v ?? a.cells[idx].textContent;
    let y = b.cells[idx].dataset.v ?? b.cells[idx].textContent;
    if (numeric) { x = parseFloat(x); y = parseFloat(y); }
    return (x > y ? 1 : x < y ? -1 : 0) * (asc ? 1 : -1);
  });
  rows.forEach(r => tbody.appendChild(r));
}
function filterRegions() {
  const q = document.getElementById('q').value.toLowerCase();
  const c = document.getElementById('cat').value;
  let shown = 0;
  document.querySelectorAll('#regions tbody tr').forEach(r => {
    const ok = (!c || r.dataset.cat === c) && (!q || r.textContent.toLowerCase().includes(q));
    r.style.display = ok ? '' : 'none';
    if (ok) shown++;
  });
  document.getElementById('shown').textContent = shown;
}
"""


def strip_svg_header(text: str) -> str:
    """Remove an SVG's XML and DOCTYPE headers so it can be inlined in HTML."""
    text = re.sub(r"<\?xml[^>]*\?>", "", text)
    text = re.sub(r"<!DOCTYPE[^>]*>", "", text, flags=re.S)
    return text.strip()


def _svg(path: str) -> str:
    """Read an SVG file and strip its headers so it can be inlined in HTML."""
    return strip_svg_header(Path(path).read_text())


def _sw(cat: str) -> str:
    """A small square in the category's colour (the same colour as in the figures)."""
    return f'<span class="swatch" style="background:{colour(cat)}"></span>'


def write_html(data, path: Path, svgs: list[str], labels: dict[str, str]) -> None:
    """Write report.html, a single self-contained page, section by section.

    Every value taken from the data is HTML-escaped, so names cannot break the page.
    """
    e = escape
    total = next(r for r in data.summary if r.record == "ALL" and r.category == "total")
    ref_name = data.reference.name
    # The page is built as a list of strings and joined once at the end.
    parts: list[str] = []
    add = parts.append

    # Header and summary tiles
    add("<!doctype html><html lang='en'><head><meta charset='utf-8'>")
    add("<meta name='viewport' content='width=device-width, initial-scale=1'>")
    add(f"<title>refmasker: {e(ref_name)}</title><style>{CSS}</style></head><body><main>")
    add(f"<h1>Reference mask for {e(ref_name)}</h1>")
    add(
        f"<p class='sub'>refmasker {e(data.versions.get('refmasker', ''))} · {e(data.started)} · "
        f"built for {data.params['read_length']} bp reads</p>"
    )

    add("<div class='tiles'>")
    for k, v in [
        ("Genome size", f"{fmt_bp(total.length)} bp"),
        ("Masked", f"{fmt_bp(total.bp)} bp"),
        ("Fraction masked", f"{total.pct:.2f}%"),
        ("Masked regions", f"{total.regions:,}"),
        ("Records", f"{len(data.records)}"),
    ]:
        add(f"<div class='tile'><div class='k'>{k}</div><div class='v'>{v}</div></div>")
    add("</div>")

    # Warnings (e.g. draft assembly, N bases, mismatched annotation names)
    for w in data.warnings:
        add(f"<div class='warn'><strong>Warning.</strong> {e(w)}</div>")

    # Figures, inlined as SVG
    add("<h2>Where the masked regions are</h2>")
    add(
        "<p>The top track shows every base in <code>mask.bed</code>; the tracks below show each "
        "category before merging. Regions can belong to more than one category.</p>"
    )
    for s in svgs:
        add(f"<div class='figure'>{_svg(s)}</div>")

    # Amount masked per category, and per record when there is more than one
    add("<h2>How much is masked</h2>")
    add(
        "<p>Each category's figure counts bases once even where its own regions overlap. "
        "The total counts each base once across all categories, so it is less than the sum of "
        "the rows when categories overlap.</p>"
    )
    add(_summary_table(data, labels, "ALL"))
    if len(data.records) > 1:
        add("<h3>By record</h3>")
        add(_record_table(data))

    # Reason and detecting tool(s) for each category, plus any notes from its detector
    add("<h2>Why each category is masked</h2><table><thead><tr><th>Category</th>")
    add("<th>Detected by</th><th>Reason</th></tr></thead><tbody>")
    for cat in data.categories:
        sources = sorted({r.source for r in data.regions if r.category == cat}) or ["(none found)"]
        note = " ".join(data.notes.get(cat, []))
        add(
            f"<tr><td>{_sw(cat)}{e(labels.get(cat, cat))}</td><td>{e(', '.join(sources))}</td>"
            f"<td>{e(WHY.get(cat, ''))}{(' <em>' + e(note) + '</em>') if note else ''}</td></tr>"
        )
    add("</tbody></table>")

    # Repeat threshold sensitivity grid (only when repeats were masked)
    if data.sensitivity:
        add("<h2>Repeat threshold sensitivity</h2>")
        add(
            "<p>Fraction of the genome masked as repeats at other minimum length and identity "
            f"settings. The setting used (bold) is {data.params['repeat_min_len']} bp at "
            f"{data.params['repeat_min_id']:g}% identity.</p>"
        )
        add(
            "<table><thead><tr><th class='n'>Min length (bp)</th><th class='n'>Min identity (%)</th>"
            "<th class='n'>Regions</th><th class='n'>bp masked</th><th class='n'>% genome</th>"
            "</tr></thead><tbody>"
        )
        for s in data.sensitivity:
            cls = " class='chosen'" if s.chosen else ""
            add(
                f"<tr{cls}><td class='n'>{s.min_len}</td><td class='n'>{s.min_id:g}</td>"
                f"<td class='n'>{s.regions:,}</td><td class='n'>{s.bp:,}</td>"
                f"<td class='n'>{s.pct:.2f}</td></tr>"
            )
        add("</tbody></table>")

    # Every region, with a search box and a category filter
    add("<h2>All regions</h2>")
    add(
        "<p>Every region before merging, after padding. Click a column heading to sort. "
        f"Showing <span id='shown'>{len(data.regions)}</span> of {len(data.regions)}.</p>"
    )
    add(
        "<div class='controls'><input id='q' placeholder='Search names or genes' oninput='filterRegions()'>"
    )
    add("<select id='cat' onchange='filterRegions()'><option value=''>All categories</option>")
    for cat in data.categories:
        add(f"<option value='{e(cat)}'>{e(labels.get(cat, cat))}</option>")
    add("</select></div>")
    add("<div class='scroll'><table id='regions'><thead><tr>")
    for h, n in [
        ("Record", False),
        ("Start", True),
        ("End", True),
        ("Length", True),
        ("Category", False),
        ("Source", False),
        ("Name", False),
        ("Genes overlapped", False),
    ]:
        add(f"<th{' class=n' if n else ''} onclick='sortTable(this)'>{h}</th>")
    add("</tr></thead><tbody>")
    # Starts are shown 1-based (as in genome browsers); data-v holds the values to sort by.
    for r in data.regions:
        add(
            f"<tr data-cat='{e(r.category)}'><td>{e(r.chrom)}</td>"
            f"<td class='n' data-v='{r.start + 1}'>{r.start + 1:,}</td>"
            f"<td class='n' data-v='{r.end}'>{r.end:,}</td>"
            f"<td class='n' data-v='{r.length}'>{r.length:,}</td>"
            f"<td>{_sw(r.category)}{e(labels.get(r.category, r.category))}</td>"
            f"<td>{e(r.source)}</td><td>{e(r.name)}</td>"
            f"<td class='genes'>{e(', '.join(r.genes))}</td></tr>"
        )
    add("</tbody></table></div>")
    add("<p class='sub'>Start is 1-based here; the BED files are 0-based.</p>")

    # How to use the mask with snippy
    add("<h2>Using the mask</h2>")
    add(f"<pre>{e(SNIPPY_USAGE)}</pre>")
    add(
        "<p>The sequence names in <code>mask.bed</code> are the FASTA IDs of this reference, so "
        "use the same reference FASTA for snippy.</p>"
    )

    # Run details for provenance: inputs, checksums, command, parameters, versions
    add("<h2>Run details</h2><dl>")
    details = [
        ("Reference", str(data.reference)),
        ("Reference md5", data.reference_md5),
        ("Reference sha256", data.reference_sha256),
        ("Annotation", str(data.annotation) if data.annotation else "none"),
        ("Records", "; ".join(f"{r.id} ({r.length:,} bp)" for r in data.records)),
        ("Plasmid records", ", ".join(data.plasmids) or "none"),
        ("Command", data.command),
        ("Started", data.started),
        ("Runtime", f"{data.runtime_s:.0f} s"),
    ]
    details += [
        (k.replace("_", " ").capitalize(), str(v)) for k, v in data.params.items() if v is not None
    ]
    details += [(f"Version: {k}", v) for k, v in data.versions.items()]
    for k, v in details:
        add(f"<dt>{e(k)}</dt><dd>{e(v)}</dd>")
    add("</dl>")
    add(f"<script>{JS}</script></main></body></html>")
    with open_replacing(path) as fh:
        fh.write("\n".join(parts))


def _summary_table(data, labels, record: str) -> str:
    """HTML table of bp and regions masked per category for one record ('ALL' = genome)."""
    rows = [r for r in data.summary if r.record == record]
    out = [
        "<table><thead><tr><th>Category</th><th class='n'>Regions</th>"
        "<th class='n'>bp masked</th><th class='n'>% of genome</th></tr></thead><tbody>"
    ]
    for r in rows:
        if r.category == "total":
            name = "<strong>Total (each base counted once)</strong>"
        else:
            name = _sw(r.category) + escape(labels.get(r.category, r.category))
        out.append(
            f"<tr><td>{name}</td><td class='n'>{r.regions:,}</td>"
            f"<td class='n'>{r.bp:,}</td><td class='n'>{r.pct:.2f}</td></tr>"
        )
    out.append("</tbody></table>")
    return "".join(out)


def _record_table(data) -> str:
    """HTML table of the total masked on each record, marking which records are plasmids."""
    rows = [r for r in data.summary if r.record != "ALL" and r.category == "total"]
    plasmids = set(data.plasmids)
    out = [
        "<table><thead><tr><th>Record</th><th>Type</th><th class='n'>Length</th>"
        "<th class='n'>bp masked</th><th class='n'>% masked</th></tr></thead><tbody>"
    ]
    for r in rows:
        kind = "plasmid" if r.record in plasmids else "chromosome/contig"
        out.append(
            f"<tr><td>{escape(r.record)}</td><td>{kind}</td><td class='n'>{r.length:,}</td>"
            f"<td class='n'>{r.bp:,}</td><td class='n'>{r.pct:.2f}</td></tr>"
        )
    out.append("</tbody></table>")
    return "".join(out)
