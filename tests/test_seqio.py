import pytest

from refmasker.seqio import read_fasta, write_fasta


def test_read_fasta(genome):
    recs = read_fasta(genome)
    assert [r.id for r in recs] == ["chr1", "pA"]
    assert recs[1].description == "Test organism plasmid pA"
    assert recs[0].length == 26012


def test_read_gzipped(genome_gz):
    assert len(read_fasta(genome_gz)) == 2


@pytest.mark.parametrize(
    "text",
    [">a\nACGT\n>a\nACGT\n", ">a\n\n", "ACGT\n", ""],
    ids=["duplicate", "empty-seq", "no-header", "empty-file"],
)
def test_read_fasta_errors(tmp_path, text):
    p = tmp_path / "x.fa"
    p.write_text(text)
    with pytest.raises(ValueError):
        read_fasta(p)


def test_write_fasta_strips_descriptions(tmp_path, genome):
    out = tmp_path / "o.fa"
    write_fasta(out, read_fasta(genome))
    headers = [line for line in out.read_text().splitlines() if line.startswith(">")]
    assert headers == [">chr1", ">pA"]
