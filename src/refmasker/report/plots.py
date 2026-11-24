"""Linear and circular genome plots of the masked regions."""

from __future__ import annotations

import io
import math

import matplotlib

# Non-interactive backend: draw straight to files, no display needed (e.g. on a cluster node).
# It must be chosen before pyplot is imported, hence the imports below it.
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

from refmasker.safeio import open_replacing  # noqa: E402

# svg.fonttype "none" keeps text as real text in SVGs (searchable, smaller files).
# text.parse_math off: a '$' in a record ID is drawn as-is instead of being parsed as maths,
# which garbles the name or, for IDs like 'ctg$\x$1', crashes the report.
plt.rcParams.update(
    {
        "svg.fonttype": "none",
        "text.parse_math": False,
        "font.size": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)

# Okabe-Ito colour-blind-safe palette.
COLOURS = {
    "repeats": "#0072B2",
    "is": "#D55E00",
    "prophage": "#009E73",
    "rrna": "#CC79A7",
    "lowcomplexity": "#E69F00",
    "crispr": "#56B4E9",
    "genes": "#B07AA1",
    "plasmid": "#8C8C3A",
    "user": "#6B6B6B",
    "all": "#222222",
}

MAX_PANELS = 6  # records drawn separately; beyond this, contigs are concatenated


def colour(cat: str) -> str:
    """Colour of a category (grey for unknown ones)."""
    return COLOURS.get(cat, "#999999")


def _scale(length: int) -> tuple[float, str]:
    """Divisor and unit for axis labels: Mb for 1 Mb or more, else kb."""
    return (1e6, "Mb") if length >= 1_000_000 else (1e3, "kb")


def _layout(data):
    """Return a list of panels: (title, length, [(chrom, offset)]).

    Up to MAX_PANELS records each get their own panel, largest first. More records (a draft
    assembly) are drawn end to end in one panel, in FASTA order; `offset` is where each
    record starts in that panel.
    """
    recs = sorted(data.records, key=lambda r: -r.length)
    if len(recs) <= MAX_PANELS:
        return [(r.id, r.length, [(r.id, 0)]) for r in recs]
    offsets, pos = [], 0
    for r in data.records:
        offsets.append((r.id, pos))
        pos += r.length
    return [(f"{len(recs)} contigs, concatenated", pos, offsets)]


def _tracks(data):
    """Track order: all masked bases first, then each category."""
    return ["all"] + list(data.categories)


def _segments(data, cat: str, offsets: dict[str, int]):
    """(start, width) bars for one track, shifted by each record's offset in the panel.

    The 'all' track uses the merged mask; category tracks use their unmerged regions.
    """
    if cat == "all":
        items = [(m.chrom, m.start, m.end) for m in data.merged]
    else:
        items = [(r.chrom, r.start, r.end) for r in data.regions if r.category == cat]
    return [(offsets[c] + s, e - s) for c, s, e in items if c in offsets]


def label(cat: str, labels: dict[str, str]) -> str:
    """Track name shown on the figures."""
    return "All masked" if cat == "all" else labels.get(cat, cat)


def linear_figure(data, labels: dict[str, str]):
    """One horizontal panel per record, with a row of bars per track ('All masked' on top)."""
    panels = _layout(data)
    tracks = _tracks(data)
    row_h = 0.3  # inches per track; each panel's height grows with the number of tracks
    heights = [len(tracks) * row_h + 0.7 for _ in panels]
    fig, axes = plt.subplots(
        len(panels),
        1,
        figsize=(11, sum(heights) + 0.6),
        gridspec_kw={"height_ratios": heights},
        squeeze=False,
    )
    for ax, (title, length, offsets) in zip(axes[:, 0], panels, strict=True):
        off = dict(offsets)
        div, unit = _scale(length)
        # Rows are drawn from the bottom up, so reverse the tracks to put 'All masked' on top.
        for i, cat in enumerate(reversed(tracks)):
            segs = _segments(data, cat, off)
            ax.broken_barh(
                segs,
                (i + 0.15, 0.7),
                facecolors=colour(cat),
                edgecolors=colour(cat),
                linewidth=0.4,
            )
        ax.set_yticks([i + 0.5 for i in range(len(tracks))])
        ax.set_yticklabels([label(c, labels) for c in reversed(tracks)])
        ax.set_ylim(0, len(tracks))
        ax.set_xlim(0, length)
        # Keep only ticks inside the record and label them in kb or Mb.
        ticks = ax.get_xticks()
        ax.set_xticks([t for t in ticks if 0 <= t <= length])
        ax.set_xticklabels([f"{t / div:g}" for t in ax.get_xticks()])
        ax.set_xlabel(f"Position ({unit})")
        ax.set_title(f"{title}  ({length:,} bp)", loc="left", fontsize=9)
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
        # In a concatenated panel, thin grey lines mark where each contig starts.
        if len(offsets) > 1:
            for _c, o in offsets[1:]:
                ax.axvline(o, color="#cccccc", linewidth=0.3, zorder=0)
    fig.tight_layout()
    return fig


def circular_figure(data, labels: dict[str, str]):
    """A circular plot per record (at most the 4 largest), one ring per track.

    'All masked' is the outermost ring. Position 0 is at the top and the genome runs clockwise.
    """
    panels = _layout(data)[:4]
    tracks = _tracks(data)
    n = len(panels)
    fig = plt.figure(figsize=(5.2 * n if n > 1 else 7.5, 6.2 if n > 1 else 7.5))
    # Each track gets an equal-width ring between radius `inner` and `outer`.
    inner, outer = 0.35, 1.0
    ring = (outer - inner) / len(tracks)
    for k, (title, length, offsets) in enumerate(panels):
        ax = fig.add_subplot(1, n, k + 1, projection="polar")
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        off = dict(offsets)
        min_w = 2 * math.pi * 0.0003  # keep tiny regions visible without exaggerating them
        for i, cat in enumerate(tracks):
            bottom = outer - (i + 1) * ring
            segs = _segments(data, cat, off)
            if segs:
                # Convert bp to angles: the record's full length is one turn (2 pi radians).
                theta = [2 * math.pi * s / length for s, _w in segs]
                width = [max(2 * math.pi * w / length, min_w) for _s, w in segs]
                ax.bar(
                    theta,
                    ring * 0.8,
                    width=width,
                    bottom=bottom + ring * 0.1,
                    align="edge",
                    color=colour(cat),
                    edgecolor=colour(cat),
                    linewidth=0.3,
                )
            # A light grey full ring behind each track, so empty tracks are still visible.
            ax.bar(
                0,
                ring * 0.8,
                width=2 * math.pi,
                bottom=bottom + ring * 0.1,
                align="edge",
                color="#f2f2f2",
                zorder=0,
            )
        div, unit = _scale(length)
        # Position labels around the circle at a round spacing.
        step = _nice_step(length / div)
        ticks = [t * step for t in range(int((length / div) // step) + 1)]
        ax.set_xticks([2 * math.pi * t * div / length for t in ticks])
        ax.set_xticklabels([f"{t:g} {unit}" if t else "0" for t in ticks], fontsize=7)
        ax.set_ylim(0, outer + 0.02)
        ax.set_yticks([])
        ax.grid(False)
        ax.spines["polar"].set_visible(False)
        ax.set_title(f"{title}\n{length:,} bp", fontsize=9, pad=14)
    handles = [Patch(color=colour(c), label=label(c, labels)) for c in tracks]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=min(len(handles), 5),
        frameon=False,
        fontsize=8,
        title="Rings, outside to inside",
        title_fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    return fig


def _nice_step(span: float) -> float:
    """Tick spacing giving about 8 ticks, rounded to 1, 2, 2.5 or 5 times a power of ten."""
    raw = span / 8
    mag = 10 ** math.floor(math.log10(raw)) if raw > 0 else 1
    for m in (1, 2, 2.5, 5, 10):
        if raw <= m * mag:
            return m * mag
    return 10 * mag


# Colours for the comparison figures: the two masks and the bases they agree on. From the
# same Okabe-Ito palette as the category colours above.
COMPARE_COLOURS = {"both": "#009E73", "only_a": "#0072B2", "only_b": "#D55E00"}

# The per-category figure cannot use the green above for "shared": the category colours are
# from the same palette, and prophage *is* that green, so a fully shared bar and an all-
# prophage bar would look identical. Shared is neutral grey there, category colour is the
# difference.
CATEGORY_SHARED_COLOUR = "#C9CCD1"


def _compare_layout(rows):
    """Panels for the comparison figures: (title, span, {chrom: offset}).

    Up to MAX_PANELS sequences get a panel each, longest first; more are drawn end to end in
    one panel. `span` is the record length when a reference was given, otherwise how far the
    masks reach (see RecordComparison.span).
    """
    rows = [r for r in rows if r.span > 0]
    if len(rows) <= MAX_PANELS:
        return [(r.record, r.span, {r.record: 0}) for r in rows]
    offsets, pos = {}, 0
    for r in rows:
        offsets[r.record] = pos
        pos += r.span
    return [(f"{len(rows)} sequences, end to end", pos, offsets)]


def _compare_segs(regions, offsets):
    """(start, width) bars for one set of regions, shifted by each record's panel offset."""
    return [(offsets[r.chrom] + r.start, r.length) for r in regions if r.chrom in offsets]


def _trim(name: str, width: int = 26) -> str:
    """Shorten a file name for a figure label, keeping the end (where the name usually is)."""
    return name if len(name) <= width else "..." + name[-(width - 3) :]


def compare_bar_figure(c, name_a: str, name_b: str):
    """How much each mask covers, split into bases the other mask shares and bases it does not."""
    fig, ax = plt.subplots(figsize=(11, 2.0))
    div, unit = _scale(max(c.bp_a, c.bp_b))
    # Drawn bottom-up, so B is listed first to put A on top.
    for i, (name, unique, key) in enumerate(
        [(name_b, c.bp_only_b, "only_b"), (name_a, c.bp_only_a, "only_a")]
    ):
        ax.barh(i, c.bp_both / div, height=0.5, color=COMPARE_COLOURS["both"])
        ax.barh(i, unique / div, left=c.bp_both / div, height=0.5, color=COMPARE_COLOURS[key])
        total = c.bp_both + unique
        ax.text(
            total / div,
            i,
            f"  {total:,} bp",
            va="center",
            ha="left",
            fontsize=8,
            color="#444444",
        )
        ax.text(0, i + 0.42, _trim(name), va="bottom", ha="left", fontsize=8.5)
    ax.set_yticks([])
    ax.set_ylim(-0.5, 1.9)
    # Leave room on the right for the bp labels.
    ax.set_xlim(0, max(c.bp_a, c.bp_b) / div * 1.22 or 1)
    ax.set_xlabel(f"Bases masked ({unit})")
    ax.legend(
        handles=[
            Patch(facecolor=COMPARE_COLOURS["both"], label="masked by both"),
            Patch(facecolor=COMPARE_COLOURS["only_a"], label="only A"),
            Patch(facecolor=COMPARE_COLOURS["only_b"], label="only B"),
        ],
        loc="lower right",
        ncol=3,
        frameon=False,
        fontsize=8,
    )
    fig.tight_layout()
    return fig


def compare_tracks_figure(c, name_a: str, name_b: str):
    """Where the two masks agree and differ along each sequence: one row per mask."""
    panels = _compare_layout(c.per_record)
    heights = [1.5 for _ in panels]
    fig, axes = plt.subplots(
        len(panels),
        1,
        figsize=(11, sum(heights) + 0.8),
        gridspec_kw={"height_ratios": heights},
        squeeze=False,
    )
    known = bool(c.lengths)
    for ax, (title, span, offsets) in zip(axes[:, 0], panels, strict=True):
        div, unit = _scale(span)
        # Drawn bottom-up, so B is listed first to put A on top. Each row holds the shared
        # bases plus that mask's own, so the two colours in a row add up to the whole mask.
        for i, (diff, own) in enumerate([(c.diff_b, "only_b"), (c.diff_a, "only_a")]):
            ax.broken_barh(
                _compare_segs(c.shared, offsets),
                (i + 0.2, 0.6),
                facecolors=COMPARE_COLOURS["both"],
            )
            ax.broken_barh(
                _compare_segs(diff, offsets), (i + 0.2, 0.6), facecolors=COMPARE_COLOURS[own]
            )
        ax.set_yticks([0.5, 1.5])
        ax.set_yticklabels([f"B: {_trim(name_b, 18)}", f"A: {_trim(name_a, 18)}"])
        ax.set_ylim(0, 2)
        ax.set_xlim(0, span)
        ticks = [t for t in ax.get_xticks() if 0 <= t <= span]
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{t / div:g}" for t in ticks])
        ax.set_xlabel(
            f"Position ({unit})"
            if known
            else f"Position ({unit}) - axis ends at the last masked base, not the record end"
        )
        ax.set_title(f"{title}  ({span:,} bp)", loc="left", fontsize=9)
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
        if len(offsets) > 1:
            for off in sorted(offsets.values())[1:]:
                ax.axvline(off, color="#cccccc", linewidth=0.3, zorder=0)
    # Legend above the first panel, so the figure also reads on its own as a saved PNG.
    axes[0, 0].legend(
        handles=[
            Patch(facecolor=COMPARE_COLOURS["both"], label="both"),
            Patch(facecolor=COMPARE_COLOURS["only_a"], label="only A"),
            Patch(facecolor=COMPARE_COLOURS["only_b"], label="only B"),
        ],
        loc="lower right",
        bbox_to_anchor=(1, 1.02),
        ncol=3,
        frameon=False,
        fontsize=8,
    )
    fig.tight_layout()
    return fig


def compare_category_figure(c, name_a: str, name_b: str, labels: dict[str, str]):
    """Per mask category, how much of it the other mask also covers.

    One panel per mask that came with an annotated BED. Each bar is one category: the part
    the other mask shares, then the part only this mask has, in that category's colour from
    the run report.
    """
    panels = [
        (name_a, name_b, c.categories_a),
        (name_b, name_a, c.categories_b),
    ]
    panels = [p for p in panels if p[2]]
    heights = [0.38 * len(rows) + 1.0 for _, _, rows in panels]
    fig, axes = plt.subplots(
        len(panels),
        1,
        figsize=(11, sum(heights)),
        gridspec_kw={"height_ratios": heights},
        squeeze=False,
    )
    for ax, (owner, other, rows) in zip(axes[:, 0], panels, strict=True):
        div, unit = _scale(max(r.bp for r in rows))
        # Rows are drawn from the bottom up, so reverse to put the biggest difference on top.
        shown = list(reversed(rows))
        for i, r in enumerate(shown):
            ax.barh(i, r.bp_shared / div, height=0.6, color=CATEGORY_SHARED_COLOUR)
            ax.barh(
                i, r.bp_only / div, left=r.bp_shared / div, height=0.6, color=colour(r.category)
            )
            ax.text(
                r.bp / div,
                i,
                f"  {r.bp_only:,} bp ({r.pct_missed:.0f}%)",
                va="center",
                ha="left",
                fontsize=8,
                color="#444444",
            )
        ax.set_yticks(range(len(shown)))
        ax.set_yticklabels([labels.get(r.category, r.category) for r in shown])
        ax.set_ylim(-0.6, len(shown) - 0.4)
        ax.set_xlim(0, max(r.bp for r in rows) / div * 1.3)
        ax.set_xlabel(f"Bases masked ({unit})")
        # The title carries the key, so the figure still reads when saved on its own.
        ax.set_title(
            f"{_trim(owner)} by category: grey is what {_trim(other)} also masks, "
            "colour is what it misses",
            loc="left",
            fontsize=9,
        )
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
        ax.legend(
            handles=[Patch(facecolor=CATEGORY_SHARED_COLOUR, label=f"also in {_trim(other, 16)}")],
            loc="lower right",
            bbox_to_anchor=(1, 1.0),
            frameon=False,
            fontsize=8,
        )
    fig.tight_layout()
    return fig


def svg_string(fig) -> str:
    """Render a figure to an SVG string, for inlining without writing a file first."""
    buf = io.BytesIO()
    fig.savefig(buf, format="svg", bbox_inches="tight")
    return buf.getvalue().decode("utf-8")


def save(fig, base, formats: set[str]) -> dict[str, str]:
    """Save `fig` as SVG always (embedded in the HTML) and PNG when asked."""
    paths = {}
    with open_replacing(f"{base}.svg", "wb") as fh:
        fig.savefig(fh, format="svg", bbox_inches="tight")
    paths["svg"] = f"{base}.svg"
    if "png" in formats:
        with open_replacing(f"{base}.png", "wb") as fh:
            fig.savefig(fh, format="png", dpi=200, bbox_inches="tight")
        paths["png"] = f"{base}.png"
    return paths


def save_png(fig, path) -> str:
    """Save `fig` as a PNG only, for a report that inlines its figures instead of linking them."""
    with open_replacing(path, "wb") as fh:
        fig.savefig(fh, format="png", dpi=200, bbox_inches="tight")
    return str(path)


def close(fig) -> None:
    """Free a figure's memory."""
    plt.close(fig)
