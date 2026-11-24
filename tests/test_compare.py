"""`refmasker compare`: the comparison maths, its warnings, and the HTML report."""

from pathlib import Path

import pytest

from refmasker.cli import compare_output_base, main
from refmasker.compare import compare, compare_files
from refmasker.regions import Region

# Two masks over one chromosome: a partial overlap, an identical region, and one region
# unique to each. A: 2000 + 2000 + 500 bp. B: 2000 + 2000 + 1000 bp. Shared: 1500 + 2000.
MASK_A = "chr1\t1000\t3000\nchr1\t10000\t12000\nchr1\t50000\t50500\n"
MASK_B = "chr1\t1500\t3500\nchr1\t10000\t12000\nchr1\t80000\t81000\n"


@pytest.fixture
def masks(tmp_path):
    a, b = tmp_path / "a.bed", tmp_path / "published.bed"
    a.write_text(MASK_A)
    b.write_text(MASK_B)
    return a, b


@pytest.fixture
def reference(tmp_path):
    """A 200 kb chromosome, long enough to hold every region in the two masks."""
    path = tmp_path / "ref.fna"
    path.write_text(">chr1 test chromosome\n" + "ACGT" * 50_000 + "\n")
    return path


def test_bp_accounting(masks):
    c = compare_files(*masks)
    assert (c.bp_a, c.bp_b, c.bp_both) == (4500, 5000, 3500)
    assert c.bp_only_a == 1000 and c.bp_only_b == 1500
    assert c.jaccard == pytest.approx(3500 / 6000)


def test_per_base_difference_and_whole_missing_regions_are_not_the_same(masks):
    c = compare_files(*masks)
    # Per-base: the 500 bp overhang of the first region, plus the region B has not got.
    assert sorted((r.start, r.end) for r in c.diff_a) == [(1000, 1500), (50000, 50500)]
    # Whole regions B does not touch at all: only the last one.
    assert [(r.start, r.end) for r in c.only_a] == [(50000, 50500)]


def test_shared_intervals(masks):
    c = compare_files(*masks)
    assert [(r.start, r.end) for r in c.shared] == [(1500, 3000), (10000, 12000)]


def test_percentages_need_a_reference(masks, reference):
    assert compare_files(*masks).pct(4500) is None
    c = compare_files(*masks, ref=reference)
    assert c.genome_length == 200_000
    assert c.pct(c.bp_a) == pytest.approx(2.25)


def test_per_record_rows(masks, reference):
    c = compare_files(*masks, ref=reference)
    (row,) = c.per_record
    assert row.record == "chr1"
    assert (row.bp_a, row.bp_b, row.bp_both) == (4500, 5000, 3500)
    assert row.length == 200_000 and row.span == 200_000
    assert row.pct(row.bp_a) == pytest.approx(2.25)


def test_without_a_reference_the_span_is_the_masked_extent(masks):
    (row,) = compare_files(*masks).per_record
    assert row.length is None
    assert row.span == row.extent == 81_000  # the furthest masked base, in either mask


def test_warns_when_the_masks_share_no_sequence_names():
    c = compare([Region("chr1", 0, 10)], [Region("NC_000962.3", 0, 10)])
    assert c.bp_both == 0
    assert "share no sequence names" in " ".join(c.warnings)
    assert "--rename" in " ".join(c.warnings)


def test_warns_about_sequences_only_one_mask_names():
    c = compare([Region("chr1", 0, 10), Region("pA", 0, 10)], [Region("chr1", 0, 10)])
    assert "pA" in " ".join(c.warnings)


def test_warns_when_regions_run_past_the_reference(masks, tmp_path):
    short = tmp_path / "short.fna"
    short.write_text(">chr1\n" + "ACGT" * 100 + "\n")  # 400 bp, far shorter than the masks
    c = compare_files(*masks, ref=short)
    assert "run past the end" in " ".join(c.warnings)


def test_rename_makes_the_masks_comparable(tmp_path):
    a, b = tmp_path / "a.bed", tmp_path / "b.bed"
    a.write_text("chr1\t0\t100\n")
    b.write_text("NC_000962.3\t0\t100\n")
    assert compare_files(a, b).bp_both == 0
    assert compare_files(a, b, {"NC_000962.3": "chr1"}).bp_both == 100


@pytest.mark.parametrize(
    "a,b,prefix,expected",
    [
        ("ref_mask/mask.bed", "snippy_tb.bed", None, "compare_mask_vs_snippy_tb"),
        ("a.bed", "b.bed", "vs_snippy", "vs_snippy.compare"),
        # Both inputs called mask.bed: the folder name keeps the two reports apart.
        ("ST8/mask.bed", "ST15/mask.bed", None, "compare_ST8_mask_vs_ST15_mask"),
        # Characters that have no business in a file name are replaced.
        ("we ird$.bed", "b.bed", None, "compare_we_ird_vs_b"),
    ],
)
def test_output_base_naming(a, b, prefix, expected):
    assert compare_output_base(a, b, prefix) == Path(expected)


def test_cli_writes_the_report_by_default(masks, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert main(["compare", str(masks[0]), str(masks[1])]) == 0
    page = tmp_path / "compare_a_vs_published.html"
    assert page.is_file()
    text = page.read_text()
    assert text.startswith("<!doctype html>")
    assert "Mask comparison: a.bed vs published.bed" in text
    assert "<svg" in text  # both figures inlined, so the page needs no other file
    assert "No reference FASTA was given" in text  # the unknown-length caveat


def test_cli_report_uses_the_reference_when_given(masks, reference, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert main(["compare", str(masks[0]), str(masks[1]), "--ref", str(reference)]) == 0
    text = (tmp_path / "compare_a_vs_published.html").read_text()
    assert "No reference FASTA was given" not in text
    assert "2.25%" in text  # 4,500 bp of a 200 kb genome


def test_cli_out_prefix_names_every_output(masks, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert main(["compare", str(masks[0]), str(masks[1]), "--out-prefix", "vs_pub"]) == 0
    assert (tmp_path / "vs_pub.compare.html").is_file()
    assert (tmp_path / "vs_pub.only_a.bed").read_text() == "chr1\t50000\t50500\n"
    assert (tmp_path / "vs_pub.only_b.bed").read_text() == "chr1\t80000\t81000\n"


def test_cli_formats_none_writes_nothing(masks, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert main(["compare", str(masks[0]), str(masks[1]), "--formats", "none"]) == 0
    assert not list(tmp_path.glob("*.html"))


def test_cli_png_writes_figures_but_no_stray_svg(masks, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    args = ["compare", str(masks[0]), str(masks[1]), "--formats", "html,png"]
    assert main(args) == 0
    assert (tmp_path / "compare_a_vs_published.bar.png").is_file()
    assert (tmp_path / "compare_a_vs_published.tracks.png").is_file()
    # The page inlines its figures, so no SVG file should be left behind.
    assert not list(tmp_path.glob("*.svg"))


def test_cli_rejects_pdf(masks, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert main(["compare", str(masks[0]), str(masks[1]), "--formats", "pdf"]) == 2


def test_hostile_sequence_name_is_escaped(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    a, b = tmp_path / "a.bed", tmp_path / "b.bed"
    a.write_text("<script>alert(1)</script>\t0\t100\n")
    b.write_text("<script>alert(1)</script>\t50\t150\n")
    assert main(["compare", str(a), str(b)]) == 0
    text = (tmp_path / "compare_a_vs_b.html").read_text()
    assert "<script>alert(1)</script>" not in text.replace("<script>function", "")
    assert "&lt;script&gt;alert(1)&lt;/script&gt;" in text


def test_empty_masks_do_not_crash_the_report(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    a, b = tmp_path / "a.bed", tmp_path / "b.bed"
    a.write_text("")
    b.write_text("")
    assert main(["compare", str(a), str(b)]) == 0
    text = (tmp_path / "compare_a_vs_b.html").read_text()
    assert "nothing to draw" in text
