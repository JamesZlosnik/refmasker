"""Linear and circular genome plots of the masked regions."""

from __future__ import annotations

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


def close(fig) -> None:
    """Free a figure's memory."""
    plt.close(fig)
