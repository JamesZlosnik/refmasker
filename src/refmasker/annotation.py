"""Read GFF3 or GenBank annotations and look up the genes a region overlaps."""

from __future__ import annotations

import re
from bisect import bisect_left
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote

from refmasker.seqio import open_text

# Feature types reported as "genes" in the region table.
GENE_TYPES = {"CDS", "rRNA", "tRNA", "ncRNA", "tmRNA", "pseudogene"}

# Feature types that are themselves mobile elements.
MOBILE_TYPES = {
    "mobile_genetic_element",
    "mobile_element",
    "insertion_sequence",
    "transposable_element",
}

# Product names that mark transposase / IS genes (PGAP, Prokka and Bakta wording).
TRANSPOSASE_RE = re.compile(
    r"transposase|insertion sequence|insertion element|\bIS[0-9]+[A-Za-z0-9]*\b", re.IGNORECASE
)


@dataclass(frozen=True)
class Feature:
    seqid: str
    start: int  # 0-based
    end: int  # exclusive
    strand: str
    ftype: str
    locus_tag: str = ""
    gene: str = ""
    product: str = ""

    @property
    def label(self) -> str:
        tag = self.locus_tag or self.gene or self.ftype
        if self.gene and self.locus_tag:
            tag = f"{self.locus_tag}({self.gene})"
        return tag

    @property
    def is_mobile(self) -> bool:
        if self.ftype in MOBILE_TYPES:
            return True
        return self.ftype in {"CDS", "pseudogene"} and bool(TRANSPOSASE_RE.search(self.product))


def detect_format(path: str | Path) -> str:
    with open_text(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("LOCUS"):
                return "genbank"
            if line.startswith("##gff-version") or len(line.split("\t")) == 9:
                return "gff3"
            if line.startswith("#"):
                continue
            break
    raise ValueError(f"{path}: could not tell whether this is GFF3 or GenBank")


def load_annotation(path: str | Path) -> list[Feature]:
    fmt = detect_format(path)
    return parse_gff3(path) if fmt == "gff3" else parse_genbank(path)


def _gff_attrs(text: str) -> dict[str, str]:
    attrs = {}
    for item in text.strip().strip(";").split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            attrs[k.strip()] = unquote(v.strip())
    return attrs


def parse_gff3(path: str | Path) -> list[Feature]:
    """Parse a GFF3 file from NCBI (PGAP), Prokka, Bakta or Ensembl.

    CDS/RNA lines often carry only a product; their locus tag and gene name live on the
    parent gene line, so names are filled in by walking up the Parent chain.
    """
    rows = []
    by_id: dict[str, dict[str, str]] = {}
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                continue
            a = _gff_attrs(cols[8])
            a["_type"] = cols[2]
            if "ID" in a:
                by_id[a["ID"]] = a
            rows.append((cols, a))

    def lineage(a: dict[str, str]) -> list[dict[str, str]]:
        chain = [a]
        seen = set()
        while "Parent" in chain[-1] and len(chain) < 4:
            pid = chain[-1]["Parent"].split(",")[0]
            if pid in seen or pid not in by_id:
                break
            seen.add(pid)
            chain.append(by_id[pid])
        return chain

    def first(chain: list[dict[str, str]], *keys: str) -> str:
        for a in chain:
            for k in keys:
                if a.get(k):
                    return a[k]
        return ""

    features = []
    for cols, a in rows:
        seqid, _src, ftype, start, end, _score, strand, _phase, _attr = cols
        if ftype == "gene" and a.get("pseudo") == "true":
            ftype = "pseudogene"
        if ftype not in GENE_TYPES and ftype not in MOBILE_TYPES:
            continue
        chain = lineage(a)
        genes = [x for x in chain if x["_type"] in {"gene", "pseudogene"}]
        gene_name = first(genes, "gene", "Name") if genes else a.get("gene", "")
        locus = first(chain, "locus_tag") or first(genes, "gene_id")
        if gene_name == locus:
            gene_name = ""
        features.append(
            Feature(
                seqid=seqid,
                start=int(start) - 1,
                end=int(end),
                strand=strand,
                ftype=ftype,
                locus_tag=locus,
                gene=gene_name,
                product=first(chain, "product", "mobile_element_type", "description"),
            )
        )
    return features


def parse_genbank(path: str | Path) -> list[Feature]:
    from Bio import SeqIO

    features = []
    with open_text(path) as fh:
        for rec in SeqIO.parse(fh, "genbank"):
            for f in rec.features:
                ftype = f.type
                q = f.qualifiers
                if ftype == "gene" and "pseudo" in q:
                    ftype = "pseudogene"
                if ftype not in GENE_TYPES and ftype not in MOBILE_TYPES:
                    continue
                product = q.get("product", q.get("mobile_element_type", [""]))[0]
                features.append(
                    Feature(
                        seqid=rec.id,
                        start=int(f.location.start),
                        end=int(f.location.end),
                        strand={1: "+", -1: "-"}.get(f.location.strand, "."),
                        ftype=ftype,
                        locus_tag=q.get("locus_tag", [""])[0],
                        gene=q.get("gene", [""])[0],
                        product=product,
                    )
                )
    return features


class FeatureIndex:
    """Fast overlap lookup of gene features per sequence."""

    def __init__(self, features: list[Feature]):
        self._by_seq: dict[str, list[Feature]] = defaultdict(list)
        for f in features:
            if f.ftype in GENE_TYPES:
                self._by_seq[f.seqid].append(f)
        self._starts: dict[str, list[int]] = {}
        self._maxlen: dict[str, int] = {}
        for seqid, fs in self._by_seq.items():
            fs.sort(key=lambda f: (f.start, f.end))
            self._starts[seqid] = [f.start for f in fs]
            self._maxlen[seqid] = max(f.end - f.start for f in fs)

    @property
    def seqids(self) -> set[str]:
        return set(self._by_seq)

    def overlapping(self, seqid: str, start: int, end: int) -> list[Feature]:
        fs = self._by_seq.get(seqid)
        if not fs:
            return []
        starts = self._starts[seqid]
        lo = bisect_left(starts, start - self._maxlen[seqid])
        hi = bisect_left(starts, end)
        return [f for f in fs[lo:hi] if f.end > start]
