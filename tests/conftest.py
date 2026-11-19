import gzip
import random
from pathlib import Path

import pytest


def random_seq(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


@pytest.fixture
def genome(tmp_path: Path) -> Path:
    """A 30 kb chromosome with a 1 kb repeat planted twice, plus a 5 kb plasmid."""
    rep = random_seq(1000, 1)
    chrom = random_seq(5000, 2) + rep + random_seq(14000, 3) + rep + random_seq(5000, 4) + "A" * 12
    plasmid = random_seq(5000, 5)
    path = tmp_path / "ref.fna"
    path.write_text(
        f">chr1 Test organism chromosome\n{chrom}\n>pA Test organism plasmid pA\n{plasmid}\n"
    )
    return path


@pytest.fixture
def genome_gz(tmp_path: Path, genome: Path) -> Path:
    out = tmp_path / "ref.fna.gz"
    with gzip.open(out, "wt") as fh:
        fh.write(genome.read_text())
    return out
