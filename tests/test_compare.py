"""`refmasker compare`: the comparison maths, its warnings, and the HTML report."""

from pathlib import Path

import pytest

from refmasker.cli import compare_output_base, main
from refmasker.compare import annotated_source, compare, compare_files
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


# ---------------------------------------------------------------- mask categories
#
# A plain mask.bed has three columns, so the categories behind each difference come from a
# mask.annotated.bed: the input itself, the file beside it, or --annotated-a/-b.

ANNOTATED_A = (
    "#chrom\tstart\tend\tcategory\tsource\tlength\tname\tgenes\n"
    "chr1\t1000\t3000\trepeats\tnucmer\t2000\tmatches 1 other copy\t.\n"
    "chr1\t10000\t12000\tis\tISEScan\t2000\tIS6110\t.\n"
    "chr1\t50000\t50500\tprophage\tgeNomad\t500\tscore=0.9\t.\n"
)


@pytest.fixture
def annotated_masks(tmp_path):
    """A refmasker output folder (mask.bed with its annotated twin) and a plain published mask."""
    folder = tmp_path / "ST8_mask"
    folder.mkdir()
    (folder / "mask.bed").write_text(MASK_A)
    (folder / "mask.annotated.bed").write_text(ANNOTATED_A)
    published = tmp_path / "published.bed"
    published.write_text(MASK_B)
    return folder / "mask.bed", published


def test_annotated_source_finds_the_file_beside_the_mask(annotated_masks):
    mask, published = annotated_masks
    assert annotated_source(mask, None) == mask.with_suffix(".annotated.bed")
    assert annotated_source(published, None) is None


def test_annotated_source_accepts_an_annotated_bed_as_the_mask(annotated_masks):
    annotated = annotated_masks[0].with_suffix(".annotated.bed")
    assert annotated_source(annotated, None) == annotated


def test_annotated_source_override_and_off(annotated_masks):
    mask, published = annotated_masks
    annotated = mask.with_suffix(".annotated.bed")
    assert annotated_source(published, str(annotated)) == annotated
    assert annotated_source(mask, "none") is None


def test_annotated_source_rejects_a_bad_override(annotated_masks, tmp_path):
    mask, published = annotated_masks
    with pytest.raises(FileNotFoundError):
        annotated_source(mask, str(tmp_path / "nope.bed"))
    with pytest.raises(ValueError, match="not a refmasker mask.annotated.bed"):
        annotated_source(mask, str(published))


def test_categories_report_what_the_other_mask_misses(annotated_masks):
    mask, published = annotated_masks
    c = compare_files(mask, published, annotated_a=mask.with_suffix(".annotated.bed"))
    rows = {r.category: r for r in c.categories_a}
    # B covers 1500-3000 of the repeat, all of the IS, and none of the prophage.
    assert (rows["repeats"].bp, rows["repeats"].bp_shared, rows["repeats"].bp_only) == (
        2000,
        1500,
        500,
    )
    assert rows["is"].bp_only == 0 and rows["is"].pct_missed == 0
    assert rows["prophage"].bp_only == 500 and rows["prophage"].pct_missed == 100
    # Biggest difference first, so the reader sees what matters at the top.
    assert [r.category for r in c.categories_a][0] == "prophage"
    assert c.categories_b == []  # the published mask has no categories to report


def test_unique_regions_carry_their_categories_and_names(annotated_masks):
    mask, published = annotated_masks
    c = compare_files(mask, published, annotated_a=mask.with_suffix(".annotated.bed"))
    (region,) = c.unique_a
    assert (region.start, region.end) == (50000, 50500)
    assert region.categories == ("prophage",)
    assert region.names == ("score=0.9",)


def test_unique_regions_without_categories_are_still_listed(masks):
    c = compare_files(*masks)
    assert [(r.start, r.end) for r in c.unique_a] == [(50000, 50500)]
    assert c.unique_a[0].categories == ()


def test_rename_keeps_the_categories(annotated_masks, tmp_path):
    mask, _ = annotated_masks
    other = tmp_path / "other.bed"
    other.write_text("NC_1\t1000\t3000\n")
    c = compare_files(mask, other, {"NC_1": "chr1"}, annotated_a=mask.with_suffix(".annotated.bed"))
    assert c.bp_both == 2000
    assert {r.category for r in c.categories_a} == {"repeats", "is", "prophage"}


def test_cli_finds_the_annotated_file_by_itself(annotated_masks, tmp_path, monkeypatch, caplog):
    monkeypatch.chdir(tmp_path)
    mask, published = annotated_masks
    assert main(["compare", str(mask), str(published)]) == 0
    assert "Categories for A from" in caplog.text
    text = (tmp_path / "compare_mask_vs_published.html").read_text()
    assert "What kind of masking the differences are" in text  # the per-category figure
    assert "by category" in text
    assert "Prophages" in text  # the category label, from the detector registry
    assert "score=0.9" in text  # the region's name, in the unique-regions table


def test_cli_can_switch_the_categories_off(annotated_masks, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    mask, published = annotated_masks
    assert main(["compare", str(mask), str(published), "--annotated-a", "none"]) == 0
    text = (tmp_path / "compare_mask_vs_published.html").read_text()
    assert "Neither mask came with categories" in text
    assert "What kind of masking the differences are" not in text


def test_cli_rejects_a_plain_bed_as_annotated(annotated_masks, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    mask, published = annotated_masks
    assert main(["compare", str(mask), str(published), "--annotated-a", str(published)]) == 2
