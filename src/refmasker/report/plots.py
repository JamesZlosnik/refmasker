"""Linear and circular genome plots of the masked regions."""

from __future__ import annotations

import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

plt.rcParams.update(
    {
        "svg.fonttype": "none",
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
    return COLOURS.get(cat, "#999999")


def _scale(length: int) -> tuple[float, str]:
    return (1e6, "Mb") if length >= 1_000_000 else (1e3, "kb")


def _layout(data):
    """Return a list of panels: (title, length, [(chrom, offset)])."""
    recs = sorted(data.records, key=lambda r: -r.length)
    if len(recs) <= MAX_PANELS:
        return [(r.id, r.length, [(r.id, 0)]) for r in recs]
    offsets, pos = [], 0
    for r in data.records:
        offsets.append((r.id, pos))
        pos += r.length
    return [(f"{len(recs)} contigs, concatenated", pos, offsets)]


def _tracks(data):
    return ["all"] + list(data.categories)


def _segments(data, cat: str, offsets: dict[str, int]):
    if cat == "all":
        items = [(m.chrom, m.start, m.end) for m in data.merged]
    else:
        items = [(r.chrom, r.start, r.end) for r in data.regions if r.category == cat]
    return [(offsets[c] + s, e - s) for c, s, e in items if c in offsets]


def label(cat: str, labels: dict[str, str]) -> str:
    return "All masked" if cat == "all" else labels.get(cat, cat)


def linear_figure(data, labels: dict[str, str]):
    panels = _layout(data)
    tracks = _tracks(data)
    row_h = 0.3
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
        ticks = ax.get_xticks()
        ax.set_xticks([t for t in ticks if 0 <= t <= length])
        ax.set_xticklabels([f"{t / div:g}" for t in ax.get_xticks()])
        ax.set_xlabel(f"Position ({unit})")
        ax.set_title(f"{title}  ({length:,} bp)", loc="left", fontsize=9)
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
        if len(offsets) > 1:
            for _c, o in offsets[1:]:
                ax.axvline(o, color="#cccccc", linewidth=0.3, zorder=0)
    fig.tight_layout()
    return fig


def circular_figure(data, labels: dict[str, str]):
    panels = _layout(data)[:4]
    tracks = _tracks(data)
    n = len(panels)
    fig = plt.figure(figsize=(5.2 * n if n > 1 else 7.5, 6.2 if n > 1 else 7.5))
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
    raw = span / 8
    mag = 10 ** math.floor(math.log10(raw)) if raw > 0 else 1
    for m in (1, 2, 2.5, 5, 10):
        if raw <= m * mag:
            return m * mag
    return 10 * mag


def save(fig, base, formats: set[str]) -> dict[str, str]:
    """Save `fig` as SVG always (embedded in the HTML) and PNG when asked."""
    paths = {}
    fig.savefig(f"{base}.svg", bbox_inches="tight")
    paths["svg"] = f"{base}.svg"
    if "png" in formats:
        fig.savefig(f"{base}.png", dpi=200, bbox_inches="tight")
        paths["png"] = f"{base}.png"
    return paths


def close(fig) -> None:
    plt.close(fig)
