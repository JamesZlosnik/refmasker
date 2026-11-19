"""Detectors: one per mask category.

Each detector takes a Context and returns a Result holding Regions tagged with
its category. Parsers for tool output are separate pure functions so they can be
unit-tested without the tools installed.
"""

from __future__ import annotations

import logging
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from refmasker.annotation import Feature
from refmasker.regions import Region
from refmasker.seqio import Record

log = logging.getLogger(__name__)


@dataclass
class Context:
    fasta: Path  # normalised FASTA (ID-only headers)
    records: list[Record]
    workdir: Path
    threads: int
    options: Any  # argparse.Namespace from the CLI
    features: list[Feature] | None = None
    annotation_path: Path | None = None

    @property
    def lengths(self) -> dict[str, int]:
        return {r.id: r.length for r in self.records}


@dataclass
class Result:
    regions: list[Region] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    extra: dict[str, Any] = field(default_factory=dict)


@dataclass
class Category:
    name: str
    label: str
    default: bool
    description: str
    tools: Callable[[Any], list[str]]
    run: Callable[[Context], Result]


def _registry() -> dict[str, Category]:
    from refmasker.detectors import (
        crispr,
        genes,
        insertion_sequences,
        lowcomplexity,
        plasmids,
        prophages,
        repeats,
        rrna,
        user,
    )

    cats = [
        Category(
            "repeats",
            "Repeats",
            True,
            "Sequence present in more than one copy (self-alignment with nucmer)",
            lambda o: ["nucmer", "show-coords"],
            repeats.detect,
        ),
        Category(
            "is",
            "IS / transposons",
            True,
            "Insertion sequences and transposases (ISEScan, plus annotation)",
            lambda o: ["isescan.py"],
            insertion_sequences.detect,
        ),
        Category(
            "prophage",
            "Prophages",
            True,
            "Integrated phage (geNomad or PhiSpy)",
            lambda o: ["genomad"] if o.phage_tool == "genomad" else ["PhiSpy.py"],
            prophages.detect,
        ),
        Category(
            "rrna",
            "rRNA operons",
            True,
            "Ribosomal RNA operons (barrnap)",
            lambda o: ["barrnap"],
            rrna.detect,
        ),
        Category(
            "lowcomplexity",
            "Low complexity",
            False,
            "Low-complexity sequence (dustmasker) and homopolymer runs",
            lambda o: ["dustmasker"],
            lowcomplexity.detect,
        ),
        Category(
            "crispr",
            "CRISPR arrays",
            False,
            "CRISPR repeat-spacer arrays (MinCED)",
            lambda o: ["minced"],
            crispr.detect,
        ),
        Category(
            "genes",
            "Named genes",
            False,
            "Annotated genes matching --mask-genes-matching",
            lambda o: [],
            genes.detect,
        ),
        Category(
            "plasmid",
            "Plasmids",
            False,
            "Whole plasmid records",
            lambda o: [],
            plasmids.detect,
        ),
        Category(
            "user",
            "User regions",
            False,
            "Regions supplied with --extra-bed",
            lambda o: [],
            user.detect,
        ),
    ]
    return {c.name: c for c in cats}


CATEGORIES: dict[str, Category] = {}


def categories() -> dict[str, Category]:
    if not CATEGORIES:
        CATEGORIES.update(_registry())
    return CATEGORIES


def default_categories() -> list[str]:
    return [c.name for c in categories().values() if c.default]
