"""Colours shared by the figures and the HTML reports.

Separate from plots.py, which imports matplotlib: the tables and swatches in a report need
the same colours, but nothing else from the plotting stack, so a report can still be written
on a machine where matplotlib is missing or broken.
"""

from __future__ import annotations

# Okabe-Ito colour-blind-safe palette, one colour per mask category.
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

# The comparison figures: the two masks and the bases they agree on.
COMPARE_COLOURS = {"both": "#009E73", "only_a": "#0072B2", "only_b": "#D55E00"}

# The per-category comparison figure cannot use the green above for "shared": the category
# colours are from the same palette, and prophage *is* that green, so a fully shared bar and
# an all-prophage bar would look identical. Shared is neutral grey there.
CATEGORY_SHARED_COLOUR = "#C9CCD1"


def colour(cat: str) -> str:
    """Colour of a category (grey for unknown ones)."""
    return COLOURS.get(cat, "#999999")
