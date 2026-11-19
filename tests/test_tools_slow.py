"""Integration tests that run the real external tools. Run with: pytest -m slow"""

import shutil

import pytest

from refmasker.cli import main

pytestmark = pytest.mark.slow


def need(*tools):
    missing = [t for t in tools if shutil.which(t) is None]
    return pytest.mark.skipif(bool(missing), reason=f"not installed: {', '.join(missing)}")


@need("nucmer", "show-coords")
def test_planted_repeat_is_found(tmp_path, genome):
    out = tmp_path / "out"
    assert (
        main(["run", str(genome), "-o", str(out), "--categories", "repeats", "--formats", "none"])
        == 0
    )
    rows = [line.split("\t") for line in (out / "mask.bed").read_text().splitlines()]
    # Random flanking bases can extend the alignment by a few bp by chance.
    assert len(rows) == 2
    for (chrom, s, e), (want_s, want_e) in zip(rows, [(5000, 6000), (20000, 21000)], strict=True):
        assert chrom == "chr1"
        assert abs(int(s) - want_s) <= 5 and abs(int(e) - want_e) <= 5


@need("barrnap", "dustmasker", "minced")
def test_other_tools_run(tmp_path, genome):
    out = tmp_path / "out"
    rc = main(
        [
            "run",
            str(genome),
            "-o",
            str(out),
            "--categories",
            "rrna,lowcomplexity,crispr",
            "--formats",
            "html",
        ]
    )
    assert rc == 0
    assert "chr1\t26000\t26012" in (out / "mask.bed").read_text()  # the poly-A tail
