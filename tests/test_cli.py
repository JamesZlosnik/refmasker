"""End-to-end runs using only categories that need no external tools."""

import json
import os

import pytest

from refmasker.cli import main
from refmasker.pipeline import resolve_categories


def test_resolve_categories():
    assert resolve_categories("default", False, False) == ["repeats", "is", "prophage", "rrna"]
    assert resolve_categories("repeats", True, True) == ["repeats", "plasmid", "user"]
    with pytest.raises(ValueError):
        resolve_categories("nonsense", False, False)
    with pytest.raises(ValueError):
        resolve_categories("user", False, False)


def test_run_end_to_end(tmp_path, genome):
    bed = tmp_path / "extra.bed"
    bed.write_text("chr1\t100\t200\tmy_region\nchr1\t150\t300\n")
    out = tmp_path / "out"
    rc = main(
        [
            "run",
            str(genome),
            "-o",
            str(out),
            "--categories",
            "plasmid",
            "--extra-bed",
            str(bed),
            "--pad",
            "10",
            "--formats",
            "html,pdf,png",
        ]
    )
    assert rc == 0
    for name in [
        "mask.bed",
        "mask.annotated.bed",
        "summary.tsv",
        "report.html",
        "report.pdf",
        "run_info.json",
        "refmasker.log",
        "figures/linear.png",
        "figures/linear.svg",
        "figures/circular.png",
    ]:
        assert (out / name).exists(), name
    assert (out / "mask.bed").read_text() == "chr1\t90\t310\npA\t0\t5000\n"
    info = json.loads((out / "run_info.json").read_text())
    assert info["masked_bp"] == 220 + 5000
    assert info["plasmids"] == ["pA"]
    assert not (out / "work").exists()
    html = (out / "report.html").read_text()
    assert "my_region" in html and "<svg" in html


def test_run_refuses_nonempty_outdir(tmp_path, genome):
    out = tmp_path / "out"
    out.mkdir()
    (out / "x").write_text("x")
    assert main(["run", str(genome), "-o", str(out), "--categories", "plasmid"]) == 2


def test_force_refuses_foreign_folder(tmp_path, genome):
    """--force must not clear a folder refmasker did not create (e.g. a Nextflow work/)."""
    out = tmp_path / "project"
    (out / "work").mkdir(parents=True)
    (out / "work" / "keep.txt").write_text("x")
    rc = main(["run", str(genome), "-o", str(out), "--categories", "plasmid", "--force"])
    assert rc == 2
    assert (out / "work" / "keep.txt").exists()


def test_force_reruns_into_previous_output(tmp_path, genome):
    out = tmp_path / "out"
    args = ["run", str(genome), "-o", str(out), "--categories", "plasmid", "--formats", "none"]
    assert main(args) == 0
    assert main(args + ["--force"]) == 0


def test_redact_hides_paths_and_descriptions(tmp_path, genome):
    bed = tmp_path / "extra.bed"
    bed.write_text("chr1\t100\t99999\n")  # runs past the end, so a note names the file
    out = tmp_path / "out"
    rc = main(
        ["run", str(genome), "-o", str(out), "--categories", "plasmid", "--extra-bed", str(bed)]
        + ["--formats", "html", "--redact"]
    )
    assert rc == 0
    folder = str(tmp_path)
    info_text = (out / "run_info.json").read_text()
    info = json.loads(info_text)
    assert info["redacted"] is True
    assert info["reference"] == "ref.fna"
    assert all("description" not in r for r in info["records"])
    assert info["tool_paths"] is None
    assert "Test organism" not in info_text
    html = (out / "report.html").read_text()
    for text in (info_text, html):
        assert folder not in text
    assert "extra.bed" in html  # the clipping note survives, with the file name only
    assert folder in (out / "refmasker.log").read_text()


def test_unredacted_keeps_descriptions(tmp_path, genome):
    out = tmp_path / "out"
    assert (
        main(["run", str(genome), "-o", str(out), "--categories", "plasmid", "--formats", "none"])
        == 0
    )
    info = json.loads((out / "run_info.json").read_text())
    assert info["redacted"] is False
    assert info["records"][0]["description"] == "Test organism chromosome"
    assert info["tool_paths"] == {}  # the plasmid category needs no external tools


def test_run_reports_missing_tools(tmp_path, genome, monkeypatch):
    monkeypatch.setenv("PATH", str(tmp_path))
    rc = main(["run", str(genome), "-o", str(tmp_path / "o"), "--categories", "repeats"])
    assert rc == 2


def test_formats_none(tmp_path, genome):
    out = tmp_path / "out"
    rc = main(["run", str(genome), "-o", str(out), "--categories", "plasmid", "--formats", "none"])
    assert rc == 0
    assert (out / "mask.bed").exists()
    assert not (out / "report.html").exists()
    assert not (out / "figures").exists()


def test_force_does_not_write_through_planted_symlink(tmp_path, genome):
    out = tmp_path / "out"
    args = ["run", str(genome), "-o", str(out), "--categories", "plasmid", "--formats", "none"]
    assert main(args) == 0
    victim = tmp_path / "victim.txt"
    victim.write_text("important\n")
    (out / "mask.bed").unlink()
    (out / "mask.bed").symlink_to(victim)
    assert main(args + ["--force"]) == 2
    assert victim.read_text() == "important\n"


def test_new_output_folder_is_closed_to_other_users(tmp_path, genome):
    out = tmp_path / "out"
    old = os.umask(0o022)
    try:
        assert main(["run", str(genome), "-o", str(out), "--categories", "plasmid"]) == 0
    finally:
        os.umask(old)
    assert out.stat().st_mode & 0o007 == 0


def test_failed_run_still_removes_work(tmp_path, genome):
    """work/ holds a copy of the reference, so it must not be left behind by a failed run."""
    bed = tmp_path / "extra.bed"
    bed.write_text("not_in_fasta\t1\t10\n")
    out = tmp_path / "out"
    rc = main(
        ["run", str(genome), "-o", str(out), "--categories", "plasmid", "--extra-bed", str(bed)]
    )
    assert rc == 2
    assert not (out / "work").exists()
    assert (out / "refmasker.log").exists()


def test_hostile_record_id_is_escaped_and_does_not_crash_reports(tmp_path):
    ref = tmp_path / "ref.fna"
    ref.write_text(">ctg$\\foo$<img_src=x>1 plasmid\n" + "ACGT" * 500 + "\n")
    out = tmp_path / "out"
    rc = main(["run", str(ref), "-o", str(out), "--categories", "plasmid", "--formats", "html,pdf"])
    assert rc == 0
    html = (out / "report.html").read_text()
    assert "<img" not in html
    assert "&lt;img" in html


def test_compare(tmp_path, capsys):
    a = tmp_path / "a.bed"
    b = tmp_path / "b.bed"
    a.write_text("c\t0\t100\nc\t500\t600\n")
    b.write_text("c\t50\t150\n")
    assert main(["compare", str(a), str(b), "--out-prefix", str(tmp_path / "cmp")]) == 0
    out = capsys.readouterr().out
    assert "Jaccard (bp):               0.200" in out
    assert (tmp_path / "cmp.only_a.bed").read_text() == "c\t500\t600\n"
