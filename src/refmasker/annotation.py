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
    """One annotated feature, with BED-style coordinates (0-based start, exclusive end)."""

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
        """Name shown in reports: 'locus_tag(gene)', else whichever exists, else the type."""
        tag = self.locus_tag or self.gene or self.ftype
        if self.gene and self.locus_tag:
            tag = f"{self.locus_tag}({self.gene})"
        return tag

    @property
    def is_mobile(self) -> bool:
        """True for mobile elements, and for CDS/pseudogenes with a transposase/IS product."""
        if self.ftype in MOBILE_TYPES:
            return True
        return self.ftype in {"CDS", "pseudogene"} and bool(TRANSPOSASE_RE.search(self.product))


def detect_format(path: str | Path) -> str:
    """Return 'genbank' or 'gff3' by looking at the first line that is not a comment.

    GenBank files start with 'LOCUS'; GFF3 has a '##gff-version' header or 9 tab-separated
    columns.
    """
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
    """Read a GFF3 or GenBank file (format detected automatically) into Features."""
    fmt = detect_format(path)
    return parse_gff3(path) if fmt == "gff3" else parse_genbank(path)


def _gff_attrs(text: str) -> dict[str, str]:
    """Parse a GFF3 attribute column ('ID=x;Name=y') into a dict, decoding %-escapes."""
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
    # First pass: keep the 9-column feature lines and index each by its ID, so a feature's
    # parent gene can be found whatever order the lines are in.
    rows = []
    by_id: dict[str, dict[str, str]] = {}
    with open_text(path) as fh:
        for line in fh:
            # Some GFF3 files (e.g. Prokka's) end with the sequence after ##FASTA; stop there.
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
        """The feature, its parent, grandparent, ... (at most 4; stops at unknown IDs or loops)."""
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
        """First non-empty value of any of `keys`, searching from the feature upwards."""
        for a in chain:
            for k in keys:
                if a.get(k):
                    return a[k]
        return ""

    # Second pass: build Features for gene-like and mobile-element lines only.
    features = []
    for cols, a in rows:
        seqid, _src, ftype, start, end, _score, strand, _phase, _attr = cols
        # NCBI marks pseudogenes as 'gene' lines with pseudo=true.
        if ftype == "gene" and a.get("pseudo") == "true":
            ftype = "pseudogene"
        if ftype not in GENE_TYPES and ftype not in MOBILE_TYPES:
            continue
        # The gene name comes from the parent gene line; the locus tag from anywhere in the
        # chain (Ensembl-style files use gene_id instead).
        chain = lineage(a)
        genes = [x for x in chain if x["_type"] in {"gene", "pseudogene"}]
        gene_name = first(genes, "gene", "Name") if genes else a.get("gene", "")
        locus = first(chain, "locus_tag") or first(genes, "gene_id")
        # Some annotators repeat the locus tag as the gene name; drop it to avoid 'X(X)' labels.
        if gene_name == locus:
            gene_name = ""
        features.append(
            Feature(
                seqid=seqid,
                start=int(start) - 1,  # GFF3 is 1-based inclusive; convert to 0-based
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
    """Parse a GenBank file with Biopython (its locations are already 0-based, end-exclusive)."""
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
                # Biopython stores each qualifier as a list; take the first value.
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
    """Fast lookup of the gene features overlapping an interval, per sequence.

    Features are sorted by start. A feature overlapping [start, end) must start before `end`,
    and cannot start more than the longest feature's length before `start`, so a binary search
    narrows the candidates before each one is checked.
    """

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
        """Sequence names that have at least one gene feature."""
        return set(self._by_seq)

    def overlapping(self, seqid: str, start: int, end: int) -> list[Feature]:
        """Gene features on `seqid` overlapping the half-open interval [start, end)."""
        fs = self._by_seq.get(seqid)
        if not fs:
            return []
        starts = self._starts[seqid]
        lo = bisect_left(starts, start - self._maxlen[seqid])
        hi = bisect_left(starts, end)
        return [f for f in fs[lo:hi] if f.end > start]
