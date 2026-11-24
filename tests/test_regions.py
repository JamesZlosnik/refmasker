import pytest

from refmasker.regions import (
    Region,
    intersect,
    intersect_bp,
    merge,
    pad_and_clip,
    read_bed,
    subtract,
    union_bp,
    write_annotated_bed,
    write_mask_bed,
)


def R(c, s, e, cat="x"):
    return Region(c, s, e, cat)


def test_merge_overlapping_and_bookended():
    m = merge([R("a", 0, 10, "x"), R("a", 5, 20, "y"), R("a", 20, 30, "z"), R("a", 40, 50)])
    assert [(x.start, x.end) for x in m] == [(0, 30), (40, 50)]
    assert m[0].categories == ("x", "y", "z")


def test_merge_keeps_chromosomes_apart():
    m = merge([R("a", 0, 10), R("b", 5, 15)])
    assert {(x.chrom, x.start, x.end) for x in m} == {("a", 0, 10), ("b", 5, 15)}


def test_merge_gap():
    assert len(merge([R("a", 0, 10), R("a", 15, 20)], gap=5)) == 1
    assert len(merge([R("a", 0, 10), R("a", 16, 20)], gap=5)) == 2


def test_union_counts_overlaps_once():
    assert union_bp([R("a", 0, 10), R("a", 5, 15), R("b", 0, 5)]) == 20


def test_pad_and_clip():
    out = pad_and_clip([R("a", 2, 10), R("a", 95, 100)], 5, {"a": 100})
    assert [(r.start, r.end) for r in out] == [(0, 15), (90, 100)]


def test_intersect_bp():
    a = [R("a", 0, 100), R("a", 200, 300)]
    b = [R("a", 50, 250), R("b", 0, 100)]
    assert intersect_bp(a, b) == 50 + 50


def test_read_bed(tmp_path):
    p = tmp_path / "x.bed"
    p.write_text("# comment\ntrack name=x\nchr1\t10\t20\tfoo\nchr1\t30\t40\n\n")
    rs = read_bed(p)
    assert [(r.chrom, r.start, r.end, r.name) for r in rs] == [
        ("chr1", 10, 20, "foo"),
        ("chr1", 30, 40, ""),
    ]


@pytest.mark.parametrize("line", ["chr1\t10", "chr1\tx\t20", "chr1\t30\t20"])
def test_read_bed_rejects_bad_lines(tmp_path, line):
    p = tmp_path / "x.bed"
    p.write_text(line + "\n")
    with pytest.raises(ValueError):
        read_bed(p)


def test_mask_bed_has_no_header(tmp_path):
    """snippy-core treats every line as data, so mask.bed must not have a header."""
    p = tmp_path / "mask.bed"
    write_mask_bed(p, merge([R("chr1", 0, 10), R("chr1", 5, 20)]))
    assert p.read_text() == "chr1\t0\t20\n"


def test_annotated_bed(tmp_path):
    p = tmp_path / "a.bed"
    write_annotated_bed(p, [Region("c", 0, 5, "repeats", "nucmer", "n1", ("g1", "g2"))])
    lines = p.read_text().splitlines()
    assert lines[0].startswith("#chrom")
    assert lines[1].split("\t") == ["c", "0", "5", "repeats", "nucmer", "5", "n1", "g1,g2"]


def test_overlaps_and_intersect_regions():
    a = [Region("c1", 0, 100), Region("c1", 200, 300)]
    b = [Region("c1", 50, 250), Region("c2", 0, 10)]
    assert intersect(a, b) == [Region("c1", 50, 100), Region("c1", 200, 250)]
    assert intersect_bp(a, b) == 100
    # Chromosomes present in only one set contribute nothing.
    assert intersect([Region("c3", 0, 10)], b) == []


def test_subtract_is_per_base():
    a = [Region("c1", 0, 100), Region("c1", 500, 600)]
    b = [Region("c1", 40, 60)]
    # The covered middle is cut out, leaving the two ends and the untouched region.
    assert subtract(a, b) == [
        Region("c1", 0, 40),
        Region("c1", 60, 100),
        Region("c1", 500, 600),
    ]
    assert subtract(a, []) == [Region("c1", 0, 100), Region("c1", 500, 600)]
    assert subtract(a, a) == []


def test_subtract_handles_regions_spanning_several_holes():
    a = [Region("c1", 0, 1000)]
    b = [Region("c1", 100, 200), Region("c1", 300, 400), Region("c1", 900, 2000)]
    assert subtract(a, b) == [
        Region("c1", 0, 100),
        Region("c1", 200, 300),
        Region("c1", 400, 900),
    ]
