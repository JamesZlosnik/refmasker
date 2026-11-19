"""Minimal FASTA handling (plain or gzipped)."""

from __future__ import annotations

import gzip
import hashlib
from dataclasses import dataclass
from pathlib import Path


def open_text(path: str | Path):
    """Open a text file, transparently decompressing gzip."""
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path)


@dataclass
class Record:
    id: str
    description: str
    seq: str

    @property
    def length(self) -> int:
        return len(self.seq)


def read_fasta(path: str | Path) -> list[Record]:
    records: list[Record] = []
    header = None
    chunks: list[str] = []
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    records.append(_make_record(header, chunks))
                header = line[1:].strip()
                chunks = []
            elif header is not None:
                chunks.append(line.strip())
            elif line.strip():
                raise ValueError(f"{path}: does not look like FASTA (no '>' header line)")
    if header is not None:
        records.append(_make_record(header, chunks))
    if not records:
        raise ValueError(f"{path}: no FASTA records found")
    seen = set()
    for r in records:
        if not r.id:
            raise ValueError(f"{path}: a FASTA record has an empty ID")
        if r.id in seen:
            raise ValueError(f"{path}: duplicate FASTA record ID '{r.id}'")
        if r.length == 0:
            raise ValueError(f"{path}: FASTA record '{r.id}' has no sequence")
        seen.add(r.id)
    return records


def _make_record(header: str, chunks: list[str]) -> Record:
    parts = header.split(None, 1)
    rid = parts[0] if parts else ""
    desc = parts[1] if len(parts) > 1 else ""
    return Record(rid, desc, "".join(chunks))


def write_fasta(path: str | Path, records: list[Record], width: int = 80) -> None:
    """Write records with ID-only headers so every external tool sees the same names."""
    with open(path, "w") as fh:
        for r in records:
            fh.write(f">{r.id}\n")
            for i in range(0, len(r.seq), width):
                fh.write(r.seq[i : i + width] + "\n")


def checksums(path: str | Path) -> dict[str, str]:
    """md5 (to match NCBI's md5checksums.txt) and sha256 of a file, read once."""
    md5, sha256 = hashlib.md5(), hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            md5.update(block)
            sha256.update(block)
    return {"md5": md5.hexdigest(), "sha256": sha256.hexdigest()}
