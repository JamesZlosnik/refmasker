"""Parsers for each external tool's output, tested on captured output snippets."""

from refmasker.detectors.crispr import parse_minced_gff
from refmasker.detectors.insertion_sequences import parse_isescan_tsv
from refmasker.detectors.lowcomplexity import find_homopolymers, parse_dust_intervals
from refmasker.detectors.plasmids import plasmid_ids
from refmasker.detectors.prophages import (
    genomad_db_version,
    parse_genomad_virus_summary,
    parse_phispy_coordinates,
)
from refmasker.detectors.repeats import hits_to_regions, parse_show_coords
from refmasker.detectors.rrna import group_operons, parse_barrnap_gff
from refmasker.seqio import Record

# show-coords -T -H -r -c -l (nucmer 4.0.1 on H37Rv)
SHOW_COORDS = "\n".join(
    [
        "1\t4411532\t1\t4411532\t4411532\t4411532\t100.00\t4411532\t4411532\t100.00\t100.00\tNC\tNC",
        "80185\t80373\t3082466\t3082654\t189\t189\t100.00\t4411532\t4411532\t0.00\t0.00\tNC\tNC",
        "80927\t82664\t2882593\t2884330\t1738\t1738\t81.16\t4411532\t4411532\t0.04\t0.04\tNC\tNC",
        "103706\t105076\t3884921\t3883550\t1371\t1372\t97.16\t4411532\t4411532\t0.03\t0.03\tNC\tNC",
        "3884921\t3883550\t103706\t105076\t1372\t1371\t97.16\t4411532\t4411532\t0.03\t0.03\tNC\tNC",
    ]
)


def test_show_coords_drops_self_hit_and_flips_reverse():
    hits = parse_show_coords(SHOW_COORDS)
    assert len(hits) == 4
    rev = hits[2]
    assert (rev.start2, rev.end2) == (3883549, 3884921)


def test_show_coords_from_open_file(tmp_path):
    path = tmp_path / "self.coords.tsv"
    path.write_text(SHOW_COORDS)
    with open(path) as fh:
        assert parse_show_coords(fh) == parse_show_coords(SHOW_COORDS)


def test_repeat_thresholds_and_both_copies():
    hits = parse_show_coords(SHOW_COORDS)
    regions = hits_to_regions(hits, 200, 90)
    assert sorted((r.start, r.end) for r in regions) == [(103705, 105076), (3883549, 3884921)]
    assert len(hits_to_regions(hits, 100, 80)) == 6
    assert regions[0].name == "matches 1 other copy, up to 97.2% identity"


def test_multicopy_repeat_is_one_region_per_locus():
    from refmasker.detectors.repeats import RepeatHit

    hits = [
        RepeatHit("c", 100, 1500, "c", 5000, 6400, 1400, 99.0),
        RepeatHit("c", 102, 1500, "c", 9000, 10398, 1398, 98.0),
        RepeatHit("c", 5000, 6400, "c", 9000, 10400, 1400, 99.5),
    ]
    regions = hits_to_regions(hits, 200, 90)
    assert [(r.start, r.end) for r in regions] == [(100, 1500), (5000, 6400), (9000, 10400)]
    assert regions[0].name == "matches 2 other copies, up to 99.0% identity"


ISESCAN = (
    "seqID\tfamily\tcluster\tisBegin\tisEnd\tisLen\tncopy4is\tstart1\tend1\tstart2\tend2\tscore"
    "\tirId\tirLen\tnGaps\torfBegin\torfEnd\tstrand\torfLen\tE-value\tE-value4copy\ttype\tov\ttir\n"
    "NC\tIS21\tIS21_259\t1\t1889\t1889\t1\t1\t14\t1876\t1889\t20\t12\t14\t0\t1\t1889\t+\t1889"
    "\t1e-10\t1e-10\tp\t1\tTTG:TTG\n"
    "NC\tIS6\tIS6_1\t889000\t890360\t1361\t16\t889000\t889027\t890333\t890360\t52\t28\t28\t0"
    "\t889100\t890300\t+\t1200\t1e-80\t1e-80\tc\t1\tGGT:GGT\n"
)


def test_isescan_complete_only_by_default(tmp_path):
    p = tmp_path / "x.tsv"
    p.write_text(ISESCAN)
    (r,) = parse_isescan_tsv(p)
    assert (r.chrom, r.start, r.end, r.name) == ("NC", 888999, 890360, "IS6/IS6_1 (complete)")
    assert len(parse_isescan_tsv(p, include_partial=True)) == 2


GENOMAD = (
    "seq_name\tlength\ttopology\tcoordinates\tn_genes\tgenetic_code\tvirus_score\tfdr"
    "\tn_hallmarks\tmarker_enrichment\ttaxonomy\n"
    "chr1|provirus_1001_41000\t40000\tProvirus\t1001-41000\t50\t11\t0.97\tNA\t8\t20.1"
    "\tViruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;;\n"
    "pA\t5000\tNo terminal repeats\tNA\t6\t11\t0.81\tNA\t1\t3.0\tViruses;;;;;;\n"
)


def test_genomad_proviruses_and_whole_sequences(tmp_path):
    p = tmp_path / "v.tsv"
    p.write_text(GENOMAD)
    rs = parse_genomad_virus_summary(p, {"chr1": 100000, "pA": 5000})
    assert [(r.chrom, r.start, r.end) for r in rs] == [("chr1", 1000, 41000), ("pA", 0, 5000)]
    assert "Caudoviricetes" in rs[0].name


def test_genomad_db_version(tmp_path):
    assert genomad_db_version(None) == "not set"
    assert genomad_db_version(tmp_path).startswith("unknown")
    (tmp_path / "version.txt").write_text("1.9\n")
    assert genomad_db_version(tmp_path) == "1.9"


def test_phispy(tmp_path):
    p = tmp_path / "prophage_coordinates.tsv"
    p.write_text("pp1\tchr1\t1000\t20000\t990\t1003\t19990\t20003\tAAA\tAAA\tatt\n")
    (r,) = parse_phispy_coordinates(p)
    assert (r.chrom, r.start, r.end, r.name) == ("chr1", 999, 20000, "pp1")


BARRNAP = """##gff-version 3
NC\tbarrnap:0.9\trRNA\t1471852\t1473380\t0\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA
NC\tbarrnap:0.9\trRNA\t1473661\t1476789\t0\t+\t.\tName=23S_rRNA;product=23S ribosomal RNA
NC\tbarrnap:0.9\trRNA\t1476903\t1476978\t1e-09\t+\t.\tName=5S_rRNA;product=5S ribosomal RNA (partial)
NC\tbarrnap:0.9\trRNA\t2000000\t2001528\t0\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA
"""


def test_barrnap_and_operons():
    genes = parse_barrnap_gff(BARRNAP)
    assert [g.name for g in genes] == ["16S", "23S", "5S(partial)", "16S"]
    ops = group_operons(genes, 1000)
    assert [(o.start, o.end) for o in ops] == [(1471851, 1476978), (1999999, 2001528)]
    assert ops[0].name == "operon 16S,23S,5S(partial)"
    assert len(group_operons(genes, 0)) == 4


def test_dustmasker_intervals_are_inclusive():
    rs = parse_dust_intervals(">NC\n1552 - 1558\n17896 - 17903\n>pA\n0 - 9\n")
    assert [(r.chrom, r.start, r.end) for r in rs] == [
        ("NC", 1552, 1559),
        ("NC", 17896, 17904),
        ("pA", 0, 10),
    ]


def test_homopolymers():
    recs = [Record("c", "", "ACGT" + "A" * 9 + "CG" + "T" * 7 + "GG")]
    rs = find_homopolymers(recs, 8)
    assert [(r.start, r.end, r.name) for r in rs] == [(4, 13, "polyA x9")]


def test_minced():
    gff = (
        "##gff-version 3\n"
        "NC\tminced:0.4.2\trepeat_region\t3119185\t3120539\t19\t.\t.\tID=CRISPR1;rpt_type=direct\n"
    )
    (r,) = parse_minced_gff(gff)
    assert (r.start, r.end, r.name) == (3119184, 3120539, "CRISPR1 (19 repeats)")


def test_plasmid_ids():
    recs = [
        Record("chr", "chromosome", "A"),
        Record("p1", "plasmid pX", "A"),
        Record("c3", "", "A"),
    ]
    assert plasmid_ids(recs) == ["p1"]
    assert plasmid_ids(recs, ["c3"]) == ["p1", "c3"]
